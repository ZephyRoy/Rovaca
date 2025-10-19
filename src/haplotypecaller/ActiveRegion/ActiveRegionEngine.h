#ifndef HC_ACTIVE_REGION_ENGINE_H
#define HC_ACTIVE_REGION_ENGINE_H
#include <iostream>
#include <list>
#include <map>
#include <memory>
#include <vector>

#include "ActiveBase.h"
#include "assemble_interface.h"
#include "bed_loader.h"
#include "fasta_loader.h"
#include "htslib/sam.h"

struct ActiveRegionIndexRecord
{
    int max_work_id = -1;
    /**  当前正在处理extend的work ID */
    int current_work_id = -1;

    int extension_current_tid = -1;

    hts_pos_t extension_pos = -1;

    int last_work_id = -1;
};
class AcitveRegionRingBuffer
{
public:
    int length = 0x400;
    int mask = length - 1;

    AcitveRegionRingBuffer()
    {
        m_ring_buffer.resize(length);
        std::fill(m_ring_buffer.begin(), m_ring_buffer.end(), 0.0);
    }

    AcitveRegionRingBuffer(int buffer_sz)
        : length(buffer_sz)
        , mask(buffer_sz - 1)
    {
        m_ring_buffer.resize(length);
        std::fill(m_ring_buffer.begin(), m_ring_buffer.end(), 0.0);
    }

    bool check_site_is_overflow(const int tid, const hts_pos_t pos)
    {
        if (tid != m_tid) return true;
        assert(pos >= m_start_coor);
        return pos - m_start_coor + 1 > mask;
    }

    void increase_likelihood(const int tid, const hts_pos_t pos, const double value)
    {
        if (tid != m_tid) return;
        if (pos < m_start_coor) return;
        if ((pos - m_start_coor + 1) > ((m_end_index - m_start_index) & mask)) {
            m_end_index = (m_start_index + (pos - m_start_coor + 1)) & (mask);
        }
        m_ring_buffer[(pos - m_start_coor + m_start_index) & mask] += value;
    }

    double get_likelihood(const int tid, const hts_pos_t pos)
    {
        if (tid != m_tid) return 0.0;
        assert(pos >= m_start_coor && pos <= get_end());
        return m_ring_buffer[(pos - m_start_coor + m_start_index) & mask];
    }

    hts_pos_t get_start() { return m_start_coor; }

    bool empty() { return m_end_index == m_start_index; }

    hts_pos_t get_end()
    {
        if (m_end_index == m_start_index)
            return -1;
        else
            return m_start_coor + ((m_end_index - m_start_index - 1) & mask);
    }

    void pop_n_element(const int tid, const hts_pos_t pos)
    {
        if (m_tid != tid) return;
        int new_idx = (pos - m_start_coor + m_start_index + 1) & mask;
        for (; m_start_index != new_idx; m_start_index = (m_start_index + 1) & mask) {
            m_ring_buffer[m_start_index] = 0;
            m_start_coor++;
        }
    }
    void pop_no_reset_n_element(const int tid, const hts_pos_t pos)
    {
        if (m_tid != tid) return;
        m_start_index = (pos - m_start_coor + m_start_index + 1) & mask;
        m_start_coor = pos + 1;
    }
    int get_tid() { return m_tid; }

    bool pop_front(int &tid, hts_pos_t &pos, double &value)
    {
        if (m_start_index == m_end_index) return false;
        tid = m_tid;
        pos = m_start_coor;
        m_start_coor++;
        value = m_ring_buffer[m_start_index];
        m_ring_buffer[m_start_index] = 0;
        m_start_index = (m_start_index + 1) & mask;
        return true;
    }

    bool reset(const int tid, const hts_pos_t pos)
    {
        if (tid < m_tid) return false;
        m_tid = tid;
        m_start_coor = pos;
        m_start_index = m_end_index = 0;
        return true;
    }

    void fill_all() { m_end_index = (m_start_index - 1) & mask; }

private:
    std::vector<double> m_ring_buffer;
    int m_tid = 0;
    hts_pos_t m_start_coor = 0;
    int m_start_index = 0;
    int m_end_index = 0;
};
class ActiveRegionEngine
{
public:
    static constexpr int c_filter_size = 50;

    static constexpr double gussian_sigma = 17;

    static constexpr double active_prob_threshold = 0.002000;

    ActiveRegionEngine(int min_region_size, int max_region_size, BedLoader *bed_loader, contig_info_t *fasta_info, bool force_non_active)
        : ActiveRegionEngine(min_region_size, max_region_size, c_filter_size, gussian_sigma, bed_loader, fasta_info, force_non_active)
    {}

    ActiveRegionEngine(int min_region_size, int max_region_size, BedLoader *bed_loader, contig_info_t *fasta_info)
        : ActiveRegionEngine(min_region_size, max_region_size, c_filter_size, gussian_sigma, bed_loader, fasta_info, true)
    {}

    ActiveRegionEngine(int min_region_size, int max_region_size, int filer_size, double sigma, BedLoader *bed_loader,
                       contig_info_t *fasta_info)
        : ActiveRegionEngine(min_region_size, max_region_size, filer_size, sigma, bed_loader, fasta_info, true)
    {}

    /* 仅用于测试*/
    ActiveRegionIndexRecord &get_record_stauts() { return record; }

    ActiveRegionEngine(int min_region_size, int max_region_size, int filer_size, double sigma, BedLoader *bed_loader,
                       contig_info_t *fasta_info, bool force_non_active);
    /**
     * @brief 将ActiveBase追加到队列里，追加后通过已有缓存高斯核扩展每个位点的ActiveValue
     *
     * @param base 待Extnension的ActiveBase的类
     * @return int
     */
    bool append(std::shared_ptr<ActiveBase> base);

    /**
     * @brief 返回内部的缓存Region
     *
     * @param force 强制清空所有的缓存Region
     * @return RawRegion*  若无缓存区间，则返回NULL
     */
    p_hc_region_active_storage poll();

    hts_pos_t find_best_suite(int tid, hts_pos_t pos);

    void update_finish(int last_work_id) { record.last_work_id = last_work_id; }

    bool is_finish() { return m_is_finish; }

    std::vector<double> &get_gaussian_kernel() { return m_gaussian_kernel; }

    std::vector<int> &get_release_id() { return release_ids; }

    void release_expired_work_id()
    {
        for (auto i : release_ids) {
            if (m_active_bases.count(i)) m_active_bases.erase(i);
        }
        release_ids.clear();
    }

    AcitveRegionRingBuffer &get_ring_buffer() { return m_ring_buffer; }

    void pop_active_region();

    void incorporate_status();

    bool isWES() { return m_bed_loader != nullptr; }

    void dummy_append(std::shared_ptr<ActiveBase> base) { m_active_bases[base->get_work_id()] = base; }

    void dummy_set_finish() { m_is_finish = true; }

    /**
     * @brief  对当前缓存的region进行切分，切分成两段，第一段的长度至少min_region_size，第二段开始点是最小的极小值。
     *
     * @return hts_pos_t
     */

    void change_target(int tid)
    {
        if (tid != interval_tid && m_bed_loader) {
            std::string &contig_name = m_fasta_info->key[tid];
            const std::map<std::string, p_bed_intervals> &all_interval = m_bed_loader->get_bed_intervals();
            auto it = all_interval.find(contig_name);
            if (it != all_interval.end())
                current_interval = it->second;
            else
                current_interval = nullptr;
            interval_tid = tid;
            interval_offset = 0;
        }
    }

    bool check_in_target(int tid, hts_pos_t pos)
    {
        if (!m_bed_loader) return true;

        if (tid != interval_tid) {
            std::string &contig_name = m_fasta_info->key[tid];
            const std::map<std::string, p_bed_intervals> &all_interval = m_bed_loader->get_bed_intervals();
            auto it = all_interval.find(contig_name);
            if (it != all_interval.end())
                current_interval = it->second;
            else
                current_interval = nullptr;
            interval_tid = tid;
            interval_offset = 0;
        }
        if (!current_interval) return false;
        while (interval_offset < current_interval->n && pos >= current_interval->end[interval_offset]) interval_offset++;

        return interval_offset < current_interval->n && pos >= current_interval->start[interval_offset];
    }

    void flush_region_WES(int tid);

    bool region_is_empty() { return regions.empty(); }

    bool check_first_is_same_chrom(int tid)
    {
        if (regions.empty()) return true;
        return tid == regions.front()->tid;
    }

    int get_interval_tid()
    {
        if (regions.empty()) return m_ring_buffer.get_tid();
        return -1;
    }

private:
    void pop_non_active_value_target(hts_pos_t pos);

    void pop_non_active_value(hts_pos_t pos);

    void refresh_ring_buffer(std::shared_ptr<ActiveBase> base, hts_pos_t pos);

private:
    BedLoader *m_bed_loader;

    contig_info_t *m_fasta_info;
    /** 缓存未处理的ActiveBase**/
    std::map<int, std::shared_ptr<ActiveBase>> m_active_bases{};

    AcitveRegionRingBuffer m_ring_buffer;

    int m_min_region_size;

    int m_max_region_size;

    bool m_force_non_active;
    /** max_work_id */
    ActiveRegionIndexRecord record;

    p_hc_region_active_storage current_region = nullptr;

    std::list<p_hc_region_active_storage> regions;

    bool m_is_finish = false;

    std::vector<double> m_gaussian_kernel;

    std::vector<int> release_ids;

    p_bed_intervals current_interval{nullptr};

    int interval_tid{-1};

    int interval_offset{0};
};
#endif