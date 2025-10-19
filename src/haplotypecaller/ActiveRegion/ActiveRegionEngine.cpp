
#include "ActiveRegionEngine.h"

#include <math.h>

#include <iostream>
#include <limits>
#include <numeric>

const double ROOT_TWO_PI = 2.5066282746310002;

const int wgs_ring_buffer_sz = 0x400;
const int wes_ring_buffer_sz = 0x400;
ActiveRegionEngine::ActiveRegionEngine(int min_region_size, int max_region_size, int filer_size, double sigma, BedLoader *bed_loader,
                                       contig_info_t *fasta_info, bool force_non_active)
    : m_bed_loader(bed_loader)
    , m_fasta_info(fasta_info)
    , m_ring_buffer(bed_loader ? wes_ring_buffer_sz : wgs_ring_buffer_sz)
    , m_min_region_size(min_region_size)
    , m_max_region_size(max_region_size)
    , m_force_non_active(force_non_active)
{
    int kernel_size = 2 * filer_size + 1;
    m_gaussian_kernel.reserve(kernel_size);
    for (int i = 0; i < kernel_size; i++) {
        double value = exp(-(i - filer_size) * (i - filer_size) / (2 * sigma * sigma)) / (ROOT_TWO_PI * sigma);
        m_gaussian_kernel.push_back(value);
    }
    double sum_gaussian = std::accumulate(m_gaussian_kernel.begin(), m_gaussian_kernel.end(), 0.0);
    std::transform(m_gaussian_kernel.begin(), m_gaussian_kernel.end(), m_gaussian_kernel.begin(),
                   [&](double x) { return x / sum_gaussian; });
}

bool ActiveRegionEngine::append(std::shared_ptr<ActiveBase> base)
{
    if (base->get_work_id() < record.current_work_id) {
        return false;
    }
    /**1.记录一下当前最大id号;
     * 2.将当前base插入;
     */
    if (base->get_work_id() == 0) {
        record.current_work_id = 0;
    }
    record.max_work_id = std::max(record.max_work_id, base->get_work_id());
    m_active_bases.insert({base->get_work_id(), base});
    if (record.current_work_id == -1) return true;
    // 迭代有效的base:
    //  1.调用incorporate_status，更新RingBuffer状态，如果RingBuffer almost full，则需要调用 pop_active_region，将结果更新region。详细见
    //  incorporate_status的注释
    //  2.调用 pop_active_region , 将输出的结果放在regions里面
    // 迭代完后将id追加到release_ids
    for (; record.current_work_id <= record.max_work_id; record.current_work_id++) {
        if (!m_active_bases.count(record.current_work_id)) break;
        incorporate_status();
        pop_active_region();
        release_ids.push_back(record.current_work_id);
    }

    if (record.current_work_id > record.last_work_id && record.last_work_id != -1) {
        record.extension_pos = m_ring_buffer.get_end();
        pop_active_region();
        if (!m_bed_loader)
            pop_non_active_value(m_fasta_info->idict[m_ring_buffer.get_tid()] - 1);
        else
            pop_non_active_value_target(m_fasta_info->idict[m_ring_buffer.get_tid()] - 1);
        m_is_finish = true;
    }

    return true;
}

#include <iomanip>

void ActiveRegionEngine::incorporate_status()
{
    // 迭代base里面的符合条件的位点（在target区间里面，且active Value > 0 )
    // 向前、后扩展c_filter_size，去累计对应ActiveValue*高斯值（根据距离计算）
    // 1、判断是否向前扩展的坐标已经让AcitveRegionRingBuffer满了，
    // 2、则需要调向AcitveRegionRingBuffer补全0；
    // 3、再调用pop_active_region,更新AcitveRegionRingBuffer。再跳到第一步
    // 4、完成扩展；
    std::shared_ptr<ActiveBase> base = m_active_bases[record.current_work_id];
    int max_chrom_lengh = m_fasta_info->idict[base->get_tid()];
    hts_pos_t actual_start = base->get_acutal_start();
    bool need_refresh = true;
    for (int offset = 0; offset < base->get_actual_size(); offset++) {
        if (!base->test_in_actual_target(offset)) continue;
        refresh_ring_buffer(base, offset + actual_start);
        double active_value = base->get_active_actual_result(offset).acitve_value;
        need_refresh = false;
        if (active_value == 0.0) continue;
        int repeat = base->get_extend_actual_length(offset);

        if (!m_bed_loader) {
            for (int filter_offset = 0; filter_offset < 2 * c_filter_size + 1; filter_offset++) {
                if (base->get_acutal_start() + offset + filter_offset - c_filter_size >= max_chrom_lengh) continue;
                m_ring_buffer.increase_likelihood(base->get_tid(), base->get_acutal_start() + offset + filter_offset - c_filter_size,
                                                  repeat * active_value * m_gaussian_kernel[filter_offset]);
            }
        }
        else {
            for (int filter_offset = c_filter_size, temp_offset = offset; filter_offset >= 0; filter_offset--, temp_offset--) {
                if (temp_offset >= 0 && !base->test_in_actual_target(temp_offset)) break;
                m_ring_buffer.increase_likelihood(base->get_tid(), base->get_acutal_start() + offset + filter_offset - c_filter_size,
                                                  repeat * active_value * m_gaussian_kernel[filter_offset]);
            }
            for (int filter_offset = c_filter_size + 1, temp_offset = offset; filter_offset < 2 * c_filter_size + 1;
                 filter_offset++, temp_offset++) {
                if (temp_offset < base->get_actual_size() && !base->test_in_actual_target(temp_offset)) break;
                if (base->get_acutal_start() + offset + filter_offset - c_filter_size >= max_chrom_lengh) continue;
                m_ring_buffer.increase_likelihood(base->get_tid(), base->get_acutal_start() + offset + filter_offset - c_filter_size,
                                                  repeat * active_value * m_gaussian_kernel[filter_offset]);
            }
        }
        record.extension_current_tid = base->get_tid();
        record.extension_pos = base->get_acutal_start() + offset - c_filter_size;
    }
    if (need_refresh) {
        refresh_ring_buffer(base, base->get_stop() - 1);
    }
}

hts_pos_t ActiveRegionEngine::find_best_suite(int tid, hts_pos_t pos)
{
    double minP = std::numeric_limits<double>::max();
    int minI = m_max_region_size - 1;
    for (int i = minI, j = 0; i >= m_min_region_size - 1; i--, j++) {
        double curr = m_ring_buffer.get_likelihood(tid, pos - j);
        if (curr < minP && pos - j + 1 <= m_ring_buffer.get_end() && curr <= m_ring_buffer.get_likelihood(tid, pos - j + 1) &&
            curr < m_ring_buffer.get_likelihood(tid, pos - j - 1)) {
            minI = i;
            minP = curr;
        }
    }
    return minI;
}

/**p_hc_region_active_storage C native struct*/
p_hc_region_active_storage new_region_active(int tid, hts_pos_t start, hts_pos_t end, bool active)
{
    p_hc_region_active_storage item = (p_hc_region_active_storage)malloc(sizeof(hc_region_active_storage));
    item->tid = tid;
    item->start_index = start;
    item->end_index = end;
    item->active = active;
    return item;
}

void ActiveRegionEngine::refresh_ring_buffer(std::shared_ptr<ActiveBase> base, hts_pos_t pos)
{
    while (base->get_tid() != m_ring_buffer.get_tid()) {
        // 清空数据
        change_target(m_ring_buffer.get_tid());
        if ((current_region && m_ring_buffer.get_end() != current_region->end_index) || (!current_region && !m_ring_buffer.empty())) {
            m_ring_buffer.fill_all();
            record.extension_pos = m_ring_buffer.get_end() > (hts_pos_t)m_fasta_info->idict[m_ring_buffer.get_tid()] - 1
                                       ? m_fasta_info->idict[m_ring_buffer.get_tid()] - 1
                                       : m_ring_buffer.get_end();
            pop_active_region();
        }
        if (!m_bed_loader)
            pop_non_active_value(m_fasta_info->idict[m_ring_buffer.get_tid()] - 1);
        else
            pop_non_active_value_target(m_fasta_info->idict[m_ring_buffer.get_tid()] - 1);
        if (current_region) regions.push_back(current_region);
        current_region = nullptr;
        m_ring_buffer.reset(m_ring_buffer.get_tid() + 1, 0);
        record.extension_pos = -1;
    }
    if (pos - c_filter_size - 1 < m_ring_buffer.get_start()) return;

    if (m_ring_buffer.check_site_is_overflow(m_ring_buffer.get_tid(), pos - c_filter_size - 1)) {
        change_target(m_ring_buffer.get_tid());
        if ((current_region && m_ring_buffer.get_end() != current_region->end_index) || (!current_region && !m_ring_buffer.empty())) {
            m_ring_buffer.fill_all();
            record.extension_pos = m_ring_buffer.get_end();
            pop_active_region();
        }
        if (!m_bed_loader)
            pop_non_active_value(pos - c_filter_size - 1);
        else
            pop_non_active_value_target(pos - c_filter_size - 1);
        record.extension_pos = pos - c_filter_size - 1;

        m_ring_buffer.reset(m_ring_buffer.get_tid(), pos - c_filter_size);
    }
    if (m_ring_buffer.check_site_is_overflow(base->get_tid(), pos + c_filter_size)) {
        record.extension_pos = pos - c_filter_size - 1;
        m_ring_buffer.fill_all();
        pop_active_region();
    }
}

void ActiveRegionEngine::pop_active_region()
{
    // 迭代新算好的active value的base对象
    // 1、计算激活状态；
    // 2, 判断激活状态是否与区间相同:
    // 2-1 激活状态相同时，就延伸current_region的end，如果超过max_region_size,则调用裁剪findBestCutSite逻辑;
    // 2-2 激活状态不同时，则追加到regions,将current_region置于空
    int tid = m_ring_buffer.get_tid();
    for (hts_pos_t pos = current_region == nullptr ? m_ring_buffer.get_start() : current_region->end_index + 1; pos <= record.extension_pos;
         pos++) {
        bool is_active = m_ring_buffer.get_likelihood(tid, pos) > active_prob_threshold;
        // if (is_active) {
        //     std::cerr << pos + 1 << " " << m_ring_buffer.get_likelihood(tid, pos) << std::endl;
        // }
        if (!check_in_target(tid, pos)) {
            if (current_region) {
                regions.push_back(current_region);
                current_region = nullptr;
            }
            m_ring_buffer.pop_n_element(tid, pos);
        }
        else if (!current_region) {
            if (m_force_non_active || is_active) {
                current_region = new_region_active(tid, pos, pos, is_active);
            }
            else {
                if (pos > m_ring_buffer.get_start()) m_ring_buffer.pop_n_element(tid, pos);
            }
        }
        else if (current_region->active != is_active || !(current_region->tid == tid && pos == current_region->end_index + 1)) {
            regions.push_back(current_region);
            m_ring_buffer.pop_n_element(tid, pos - 1);
            if (m_force_non_active) {
                current_region = new_region_active(tid, pos, pos, is_active);
            }
            else {
                current_region = nullptr;
            }
        }
        else {
            current_region->end_index++;
            if (current_region->end_index - current_region->start_index + 1 >= m_max_region_size) {
                if (current_region->active) {
                    hts_pos_t split_site = find_best_suite(tid, pos);
                    split_site += current_region->start_index;
                    if (split_site < current_region->end_index) {
                        current_region->end_index = split_site;
                        regions.push_back(current_region);
                        m_ring_buffer.pop_n_element(tid, split_site);
                        current_region = new_region_active(tid, split_site + 1, pos, is_active);
                    }
                    else {
                        regions.push_back(current_region);
                        m_ring_buffer.pop_n_element(tid, pos);
                        current_region = nullptr;
                    }
                }
                else {
                    regions.push_back(current_region);
                    m_ring_buffer.pop_n_element(tid, pos);
                    current_region = nullptr;
                }
            }
        }
    }
}

void ActiveRegionEngine::pop_non_active_value(hts_pos_t pos)
{
    hts_pos_t start = m_ring_buffer.get_start();
    if (current_region) {
        if (current_region->active) {
            regions.push_back(current_region);
            m_ring_buffer.pop_n_element(current_region->tid, current_region->end_index);
            start = current_region->end_index + 1;
            current_region = nullptr;
        }
        else {
            m_ring_buffer.pop_n_element(current_region->tid, record.extension_pos);
            if (pos <= record.extension_pos) {
                regions.push_back(current_region);
                current_region = nullptr;
                return;
            }
            else {
                if (current_region->start_index + m_max_region_size < pos) {
                    current_region->end_index = current_region->start_index + m_max_region_size - 1;
                    start = current_region->start_index + m_max_region_size;
                }
                else {
                    current_region->end_index = pos;
                    start = pos + 1;
                }
                regions.push_back(current_region);
                current_region = nullptr;
            }
        }
    }

    int tid = m_ring_buffer.get_tid();
    if (m_force_non_active) {
        for (hts_pos_t iter_pos = start; iter_pos <= pos; iter_pos += m_max_region_size) {
            if (iter_pos + m_max_region_size - 1 > pos) {
                current_region = new_region_active(tid, iter_pos, pos, 0);
                m_ring_buffer.pop_no_reset_n_element(current_region->tid, pos);
            }
            else {
                current_region = new_region_active(tid, iter_pos, iter_pos + m_max_region_size - 1, 0);
                m_ring_buffer.pop_no_reset_n_element(current_region->tid, current_region->end_index);
                regions.push_back(current_region);
                current_region = nullptr;
            }
        }
    }
    else {
        m_ring_buffer.pop_no_reset_n_element(tid, pos);
    }
}

void ActiveRegionEngine::pop_non_active_value_target(hts_pos_t pos)
{
    int tid = m_ring_buffer.get_tid();
    if (!current_interval) return;
    if (current_region) {
        if (current_region->active) {
            regions.push_back(current_region);
            m_ring_buffer.pop_n_element(current_region->tid, current_region->end_index);
            hts_pos_t next_index = current_region->end_index + 1;
            current_region = nullptr;
            if (m_force_non_active) {
                if (next_index > current_interval->start[interval_offset]) {
                    if (pos >= current_interval->end[interval_offset]) {
                        current_region = new_region_active(tid, next_index, current_interval->end[interval_offset] - 1, 0);
                        interval_offset++;
                    }
                    else
                        current_region = new_region_active(tid, next_index, pos, 0);
                    regions.push_back(current_region);
                    current_region = nullptr;
                }
            }
        }
        else {
            if (record.extension_pos != -1) m_ring_buffer.pop_n_element(current_region->tid, record.extension_pos);
            if (pos <= record.extension_pos) {
                regions.push_back(current_region);
                current_region = nullptr;
                return;
            }
            else {
                if (!current_interval) return;
                while (interval_offset < current_interval->n && current_region->end_index >= current_interval->end[interval_offset])
                    interval_offset++;
                if (interval_offset < current_interval->n && current_region->end_index >= current_interval->start[interval_offset]) {
                    if (pos >= current_interval->end[interval_offset]) {
                        current_region->end_index = current_interval->end[interval_offset] - 1;
                        interval_offset++;
                        regions.push_back(current_region);
                        current_region = nullptr;
                    }
                    else {
                        current_region->end_index = pos - 1;
                        regions.push_back(current_region);
                        current_region = new_region_active(tid, pos, pos, 0);
                        return;
                    }
                }
                else {
                    return;
                }
            }
        }
    }
    while (interval_offset < current_interval->n && pos >= current_interval->start[interval_offset]) {
        if (m_ring_buffer.get_start() > current_interval->end[interval_offset] - 1) {
            interval_offset++;
        }
        else if (pos >= current_interval->end[interval_offset]) {
            if (m_force_non_active) {
                if (record.extension_pos + 1 >= current_interval->start[interval_offset] &&
                    record.extension_pos + 1 <= current_interval->end[interval_offset] - 1)
                    current_region = new_region_active(tid, record.extension_pos + 1, current_interval->end[interval_offset] - 1, 0);
                else
                    current_region =
                        new_region_active(tid, current_interval->start[interval_offset], current_interval->end[interval_offset] - 1, 0);
                regions.push_back(current_region);
                current_region = nullptr;
            }
            interval_offset++;
        }
        else {
            if (m_force_non_active) {
                if (m_ring_buffer.get_start() <= current_interval->start[interval_offset])
                    current_region = new_region_active(tid, current_interval->start[interval_offset], pos - 1, 0);
                else
                    current_region = new_region_active(tid, m_ring_buffer.get_start(), pos - 1, 0);
                regions.push_back(current_region);
                current_region = new_region_active(tid, pos, pos, 0);
            }
            break;
        }
    }
}

void ActiveRegionEngine::flush_region_WES(int tid)
{
    change_target(tid);
    if (!current_interval) return;
    while (interval_offset < current_interval->n) {
        current_region = new_region_active(tid, current_interval->start[interval_offset], current_interval->end[interval_offset] - 1, 0);
        regions.push_back(current_region);
        current_region = nullptr;
        interval_offset++;
    }
}
p_hc_region_active_storage ActiveRegionEngine::poll()
{
    if (regions.empty()) {
        if (m_is_finish) {
            p_hc_region_active_storage ret = current_region;
            current_region = nullptr;
            return ret;
        }
        else
            return nullptr;
    }
    p_hc_region_active_storage region = regions.front();
    regions.pop_front();
    return region;
}