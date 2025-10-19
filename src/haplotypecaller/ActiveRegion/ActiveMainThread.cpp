#include "ActiveMainThread.h"

#include <algorithm>
#include <memory_resource>

#include "HcActiveBase.h"
#include "rovaca_logger.h"

void ActiveMainThreadDispatchTasks::run()
{
    bool finish = false;
    int current_tid = -1;
    hts_pos_t start = HTS_POS_MIN;
    hts_pos_t end = HTS_POS_MIN;
    hts_pos_t actual_start = HTS_POS_MIN;
    hts_pos_t actual_end = HTS_POS_MIN;
    int work_id = 0;
    int ref_len = 0;
    int max_size = 0;
    while (!finish) {
        std::shared_ptr<BamSource> bam_source = m_bam_resource->pop();
        bam_source->new_mempool();
        std::pmr::vector<bam1_t *> cover_reads(bam_source->pool);
        // std::cerr << "dispatch id" << work_id << std::endl;
        finish = fetch_block_reads();

        bool ret = split_reads(current_tid, start, end, actual_start, actual_end, cover_reads);
        finish &= (!m_block_resource->m_reads_buffer.empty() &&
                   (m_block_resource->m_reads_buffer.back().tid < current_tid ||
                    (m_block_resource->m_reads_buffer.back().tid == current_tid && m_block_resource->m_reads_buffer.back().end < end)));

        if (!ret) break;
        std::shared_ptr<char> ref_base = m_fasta_loader->get(current_tid, ref_len);
        std::shared_ptr<ActiveBaseResource> resource = m_base_resource->pop();
        get_target_biset(current_tid, actual_start, actual_end, resource->target, resource->m_target);
        // std::cerr << "start pool work id " << work_id << "tid:" << current_tid << "real:" << start << "-" << end
        //           << "actual:" << actual_start << "-" << actual_end << "m_bit:" << actual_end - actual_start << std::endl;
        max_size = std::max(max_size, (int)(actual_end - actual_start));
        boost::asio::post(*m_run_pool, [=] {
            pthread_setname_np(pthread_self(), "ActiveRegionBase");
            resource->new_mempool();
            ActiveBase *base = new HcActiveBase(resource->pool, 2, current_tid, start, end, actual_start, actual_end,
                                                m_bed_loader ? resource->target : nullptr, work_id, ref_base.get(), ref_len);

            for (auto reads : cover_reads) {
                base->process_bam_to_slot(reads, 0, 0);
            }
            base->compulte_all_likelihood();
            BaseSource basesource({resource, base, finish});
            m_reduce_queue->push(basesource);
            m_bam_resource->push(bam_source);
        });
        work_id++;
    }
    RovacaLogger::info("region over. dispatch finish max size: {} last tid {}", max_size, current_tid);
}

typedef void (*BamFunction)(bam1_t *);

bool ActiveMainThreadDispatchTasks::fetch_block_reads()
{
    std::lock_guard<std::mutex> lock(m_block_resource->m_mutex);
    int reads_num = 0;
    bool pos_changes = false;
    int file_index;

    auto cleanup_bam = [this](bam1_t *bam) {
        if (m_stream->get_type() == stream_type::transformed) {
            m_stream->deallocate_bqsr_mem(bam);
        }
        else {
            m_stream->mem_recovery(bam);
        }
    };

    do {
        if (!m_stream->has_next()) {
            m_block_resource->m_reads_buffer.back().is_out_chrom = true;
            m_block_resource->m_reads_buffer.back().is_out_interval = true;
            return true;
        }
        bam1_t *bam = m_stream->next(file_index);
        if (!bam) continue;
        pos_changes = m_block_resource->is_next_adjacent(bam);
        if (reads_num + 1 > batch_reads && pos_changes) {
            m_block_resource->finalize();
            m_block_resource->insert(bam);
            cleanup_bam(bam);
            break;
        }
        if (m_block_resource->insert(bam)) {
            cleanup_bam(bam);
            break;
        };
        cleanup_bam(bam);
        reads_num++;

    } while (1);
    return false;
}

void ActiveMainThreadDispatchTasks::split_block_resource(BamBlockList &block, int &tid, hts_pos_t &end, hts_pos_t &actual_start,
                                                         hts_pos_t &actual_end, std::pmr::vector<bam1_t *> &cover_reads,
                                                         bool cover_upstream)
{
    if (!cover_upstream) actual_start = block.start;
    for (auto read : block.bams) {
        cover_reads.push_back(read.bam);
    }
    if (block.is_out_interval) {
        actual_end = block.max_end;
        if (block.is_out_chrom) {
            end = m_fasta_info->idict[tid];
        }
        else {
            end = actual_end;
        }
    }
    else {
        actual_end = block.end + 1;
        end = block.end + 1;
    }
}

enum BlockStatus { INIT = 1, PRELOAD = 2, ADJACENT = 4 };

bool ActiveMainThreadDispatchTasks::split_reads(int &tid, hts_pos_t &start, hts_pos_t &end, hts_pos_t &actual_start, hts_pos_t &actual_end,
                                                std::pmr::vector<bam1_t *> &cover_reads)
{
    std::lock_guard<std::mutex> lock(m_block_resource->m_mutex);
    if (m_block_resource->m_reads_buffer.empty()) return false;

    if (tid == -1 || start == -1 || end == -1) {
        if (m_block_resource->m_reads_buffer.front().tid != 0) {
            tid = 0;
            start = actual_start = actual_end = 0;
            end = (hts_pos_t)m_fasta_info->idict[0];
        }
        else {
            // 将第一个元素插入cover reads
            tid = m_block_resource->m_reads_buffer.front().tid;
            start = 0;
            split_block_resource(m_block_resource->m_reads_buffer.front(), tid, end, actual_start, actual_end, cover_reads, false);
            return true;
        }
    }
    else if (end == (hts_pos_t)m_fasta_info->idict[tid]) {
        // blockList的第Current的元素里的所有reads插入到cover_reads
        tid++;
        start = 0;
        if ((int)m_fasta_info->key.size() == tid) return false;
        auto it = std::find_if(m_block_resource->m_reads_buffer.begin(), m_block_resource->m_reads_buffer.end(),
                               [&tid](BamBlockList i) { return i.tid == tid; });

        if (it != m_block_resource->m_reads_buffer.end()) {
            split_block_resource(*it, tid, end, actual_start, actual_end, cover_reads, false);
        }
        else {
            actual_start = actual_end = 0;
            end = m_fasta_info->idict[tid];
        }
    }
    else {
        // blockList的第Current的元素里的所有reads插入到cover_reads
        BlockStatus status = BlockStatus::INIT;
        start = end;
        actual_start = end;
        for (auto &block_reads : m_block_resource->m_reads_buffer) {
            if (block_reads.tid < tid || (block_reads.max_end < start && block_reads.tid == tid)) {
                continue;
            }
            else if (block_reads.tid == tid) {
                if (status & (BlockStatus::INIT | BlockStatus::ADJACENT) && block_reads.max_end > start && block_reads.end < start) {
                    for (auto it = block_reads.bams.rbegin(); it != block_reads.bams.rend(); it++) {
                        if ((*it).ref_end_pos > end && (*it).bam->core.tid == tid) {
                            cover_reads.push_back((*it).bam);
                            actual_end = std::max(actual_end, (*it).ref_end_pos);
                        }
                        else if ((*it).bam->core.pos < end - skip_bp_length) {
                            break;
                        }
                    }
                    status = BlockStatus::ADJACENT;
                }
                else if (status == BlockStatus::ADJACENT && block_reads.start >= start) {
                    split_block_resource(block_reads, tid, end, actual_start, actual_end, cover_reads, true);
                    break;
                }
                else if (status == BlockStatus::INIT && block_reads.start > start) {
                    split_block_resource(block_reads, tid, end, actual_start, actual_end, cover_reads, false);
                    status = BlockStatus::PRELOAD;
                    break;
                }
            }
            else {
                RovacaLogger::error("tid={} start={} end={} block_reads.tid {},size={} start={} end={}", tid, start, end, block_reads.tid,
                                  block_reads.bams.size(), block_reads.start, block_reads.end);

                assert(0);
                break;
            }
        }
        if (start == end) {
            actual_start = actual_end = 0;
            end = m_fasta_info->idict[tid];
        }
    }
    return true;
}

bool ActiveMainThreadDispatchTasks::get_target_biset(int tid, hts_pos_t start, hts_pos_t end, boost::dynamic_bitset<> *target, int &m)
{
    if (tid == -1 || start > end) {
        RovacaLogger::error("invalid region {}{}-{}", tid, start, end);
        return false;
    }
    if ((size_t)end - start > target->size()) {
        m = (end - start) * 2;
        target->resize(m, true);
    }
    if (m_bed_loader) {
        if (tid != interval_tid) {
            const std::map<std::string, p_bed_intervals> &all_interval = m_bed_loader->get_bed_intervals();
            std::string &contig_name = m_fasta_info->key[tid];
            auto it = all_interval.find(contig_name);
            if (it != all_interval.end())
                current_interval = it->second;
            else
                current_interval = nullptr;
            interval_tid = tid;
            interval_offset = 0;
        }
        if (!current_interval) {
            for (hts_pos_t pos = start; pos < end; pos++) {
                target->set(pos - start, 0);
            }
            return true;
        }
        for (hts_pos_t pos = start; pos < end; pos++) {
            while (interval_offset < current_interval->n && pos >= current_interval->end[interval_offset]) {
                interval_offset++;
            }

            if (interval_offset < current_interval->n && pos >= current_interval->start[interval_offset]) {
                target->set(pos - start, 1);
            }
            else {
                target->set(pos - start, 0);
            }
        }
    }
    else {
        target->set();
    }
    return true;
}

void ActiveMainThreadReduce::region_handle(int &region_id)
{
    p_hc_region_active_storage region = nullptr;

    do {
        std::shared_ptr<RegionResource> region_resouce = m_region_resource->pop();
        std::vector<RegionResult> result;
        // force_non_active_ true means GVCF
        uint32_t buffer_threshold = force_non_active_ ? region_resouce->buffer_size_ * 1 / 4 : region_resouce->buffer_size_ * 1 / 2;
        // TODO: 如果reigon tid  得分批次发
        while (region_resouce->pool_->get_used_data_size() < buffer_threshold) {
            if (!result.empty() && !engine->check_first_is_same_chrom(result.front().region->tid)) {
                break;
            }

            if ((region = engine->poll()) == nullptr) {
                break;
            }
            while (last_tid < region->tid) {
                m_fasta_loader->pop();
                last_tid++;
            }
            if (!force_non_active_ && !region->active) {
                delete region;
                continue;
            }

            std::vector<bam1_t *> reads;
            bool record_index = false;
            for (auto it = m_block_resource->m_reads_buffer.begin(); it != m_block_resource->m_reads_buffer.end();) {
                if ((it->max_end < region->start_index - read_padded_span && it->tid == region->tid) || it->tid < region->tid) {
                    std::lock_guard<std::mutex> lock(m_block_resource->m_mutex);
                    m_block_resource->m_reads_idle.splice(m_block_resource->m_reads_idle.begin(), it->bams);
                    it = m_block_resource->m_reads_buffer.erase(it);
                }
                else if (it->tid > region->tid || it->start > region->end_index + read_padded_span) {
                    break;
                }
                else {
                    for (auto read_it = block_itertor == it ? bam_itertor : it->bams.begin(); read_it != it->bams.end(); read_it++) {
                        if (__glibc_unlikely(!record_index && read_it->ref_end_pos > region->end_index - read_padded_span)) {
                            block_itertor = it;
                            bam_itertor = read_it;
                            record_index = true;
                        }
                        if (read_it->bam->core.pos > region->end_index + read_padded_span) {
                            break;
                        }
                        else if (read_it->ref_end_pos > region->start_index - read_padded_span) {
                            bam1_t *bam = region_resouce->pool_->alloc_bam_struct_pool();
                            if (!bam) {
                                bam = bam_init1();
                            }
                            if (!bam_copy1(bam, read_it->bam)) {
                                return;
                            }
                            region_resouce->pool_->peek_bam(bam);
                            reads.push_back(bam);
                        }
                    }
                    it++;
                }
            }
            result.emplace_back(reads, region, region_id);
            region_id++;
        }
        if (!result.empty()) {
            int target_len = 0;
            std::shared_ptr<char> ref_base = m_fasta_loader->get(last_tid, target_len);
            std::shared_ptr<RegionSource> source =
                std::make_shared<RegionSource>(ref_base, last_tid, target_len, result, 0, false, region_resouce);
            m_region_queue->push(source);
        }
        else {
            m_region_resource->push(region_resouce);
        }
        // TODO: Regions超过1k才能输出到最后强制输出，需要engine记录一下激活的region数量和非激活region的数量
    } while (region);
}

void ActiveMainThreadReduce::region_finish()
{
    std::vector<RegionResult> result;
    std::shared_ptr<RegionSource> source = std::make_shared<RegionSource>(nullptr, -1, -1, result, -1, true, nullptr);
    m_region_queue->push(source);
}

void ActiveMainThreadReduce::run()
{
    int work_id = 0;
    while (!engine->is_finish()) {
        // 从队列里获取计算完似然值的activeBase
        BaseSource base_source = m_dispatch_queue->pop();
        {
            std::shared_ptr<ActiveBase> base(base_source.base);
            work_id = base->get_work_id();
            if (base_source.finish) {
                engine->update_finish(work_id);
            }
            // record  pair of work resouce and work id
            m_cache.insert({work_id, base_source.base_source});
            // std::cerr << "get base id" << work_id << "start:" << base->get_start() << "end" << base->get_stop() << std::endl;
            engine->append(base);
            // Release base resource
        }
        std::vector<int> release_ids = engine->get_release_id();
        engine->release_expired_work_id();
        for (auto id : release_ids) {
            m_base_resource->push(m_cache[id]);
            m_cache.erase(id);
        }
        while (!engine->region_is_empty()) region_manager->poll(false);
        region_manager->flush(true);
    }
    region_manager->poll(true);
    region_manager->update_finish_sig();
}

bool ActiveMainThreadReduce::get_region_reads_block(rovaca::BamDataPool *pool, p_hc_region_active_storage region,
                                                    std::pmr::vector<bam1_t *> &reads)
{
    std::scoped_lock(m_block_resource->m_mutex);
    for (BamBlockList &block : m_block_resource->m_reads_buffer) {
        if (block.tid < region->tid || (block.tid == region->tid && block.max_end < region->start_index - read_padded_span)) {
            // 直接删除逻辑
            continue;
        }
        else if (block.tid > region->tid || block.start > region->end_index + read_padded_span) {
            break;
        }
        else {
            for (auto &read : block.bams) {
                if (read.bam->core.pos > region->end_index + read_padded_span) {
                    break;
                }
                else if (__glibc_unlikely(read.ref_end_pos > region->start_index - read_padded_span)) {
                    bam1_t *bam = pool->alloc_bam_struct_pool();
                    if (!bam) {
                        bam = bam_init1();
                    }
                    if (!bam_copy1(bam, read.bam)) {
                        return false;
                    }
                    pool->peek_bam(bam);
                    reads.push_back(bam);
                }
            }
        }
    }
    return true;
}
//
bool ActiveMainThreadReduce::remove_non_overlapping_reads_region_block(p_hc_region_active_storage region)
{
    std::lock_guard<std::mutex> lock(m_block_resource->m_mutex);
    for (auto it = m_block_resource->m_reads_buffer.begin(); it != m_block_resource->m_reads_buffer.end();) {
        if ((it->max_end < region->start_index - read_padded_span && it->tid == region->tid) || it->tid < region->tid) {
            m_block_resource->m_reads_idle.splice(m_block_resource->m_reads_idle.begin(), it->bams);
            it = m_block_resource->m_reads_buffer.erase(it);
        }
        else {
            break;
        }
    }
    return true;
}

void FakeConsumer::run()
{
    int last_tid = -1;
    hts_pos_t last_pos = -1;
    while (1) {
        std::shared_ptr<RegionSource> source = m_region_queue->pop();
        if (source->finish_) break;
        // 必须source不用再把resource放回
        std::shared_ptr<RegionResource> resource = source->resource_;

        int read_num = 0;
        for (const auto &result : source->regions_) {
            if (out.is_open()) {
                out << "chrom=" << hdr->target_name[result.region->tid] << " Start=" << result.region->start_index + 1
                    << " End=" << result.region->end_index + 1 << " Active=" << result.region->active
                    << " number=" << result.padding_reads.size() << std::endl;
            }
            read_num += result.padding_reads.size();

            if (result.region->tid < last_tid || (result.region->tid == last_tid && result.region->start_index < last_pos)) {
                RovacaLogger::error("out of order: prev tid={}, last pos={}, chrom={}, start={}, end={}, active={}", last_tid, last_pos,
                                  hdr->target_name[result.region->tid], result.region->start_index + 1, result.region->end_index + 1,
                                  result.region->active);
            }
            for (auto read : result.padding_reads) {
                if (read->core.n_cigar > 10) {
                    RovacaLogger::warn("incorrect cigar {} chrom={}, start={}, end={}, active={}", read->core.n_cigar,
                                     hdr->target_name[result.region->tid], result.region->start_index + 1, result.region->end_index + 1,
                                     result.region->active);
                }
            }
            last_tid = result.region->tid;
            last_pos = result.region->end_index;
        }
        if (out.is_open()) {
            out << "Reads Number=" << read_num << " Region Number=" << source->regions_.size() << std::endl;
        }

        source->release();
        resource->finalize();
        m_region_source->push(resource);
    }
    m_reference_manager->assign_finish();
}