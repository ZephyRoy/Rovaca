#include "RegionManager.h"

#include "rovaca_logger.h"

const int LITTLE_CHROM_SIZE = 500000;
const double LITTLE_CHROM_FACTOR = 1 / 32;
const double GVCF_BUFFER_THRES_FACTOR = 1 / 4;
const double VCF_BUFFER_THRES_FACTOR = 1 / 2;
const int EXTRA_FACTOR = 16;

void RegionManager::poll(bool force)
{
    p_hc_region_active_storage region = nullptr;

    bool is_ref_force = false;

    if (m_resource == nullptr) {
        m_resource = m_region_resource->pop();
        buffer_threshold =
            force_non_active_ ? m_resource->buffer_size_ * GVCF_BUFFER_THRES_FACTOR : m_resource->buffer_size_ * VCF_BUFFER_THRES_FACTOR;
    }

    while (m_resource->pool_->get_unsed_data_size() > max_used_region * EXTRA_FACTOR) {
        if (!result.empty() && !engine->check_first_is_same_chrom(result.front().region->tid)) {
            break;
        }
        if (m_fasta_loader->get_contig().idict[last_tid] < LITTLE_CHROM_SIZE &&
            m_resource->pool_->get_used_data_size() > (buffer_threshold * LITTLE_CHROM_FACTOR)) {
            is_ref_force = true;
            break;
        }
        if ((region = engine->poll()) == nullptr) {
            int extension_current_tid = engine->get_record_stauts().extension_current_tid;
            while (result.size() == 0 && last_tid < extension_current_tid) {
                RovacaLogger::info("region tid {} finished", last_tid);
                m_fasta_loader->pop();
                last_tid++;
            }
            break;
        }
        while (last_tid < region->tid) {
            RovacaLogger::info("region tid {} finished", last_tid);
            m_fasta_loader->pop();
            last_tid++;
        }
        if (!force_non_active_ && !region->active) {
            delete region;
            continue;
        }
        bool record_index = false;
        uint32_t current_used_size = m_resource->pool_->get_used_data_size();
        hts_pos_t padding_region_end = region->end_index + read_padded_span;
        hts_pos_t padding_region_start = region->start_index - read_padded_span;
        hts_pos_t maybe_next_region_start = region->end_index - read_padded_span;
        for (auto it = m_block_resource->m_reads_buffer.begin(); it != m_block_resource->m_reads_buffer.end();) {
            if ((it->max_end < padding_region_start && it->tid == region->tid) || it->tid < region->tid) {
                std::lock_guard<std::mutex> lock(m_block_resource->m_mutex);
                m_block_resource->m_reads_idle.splice(m_block_resource->m_reads_idle.begin(), it->bams);
                it = m_block_resource->m_reads_buffer.erase(it);
            }
            else if (it->tid > region->tid || it->start > padding_region_end) {
                break;
            }
            else {
                for (auto read_it = block_itertor == it ? bam_itertor : it->bams.begin(); read_it != it->bams.end(); read_it++) {
                    if (__glibc_unlikely(!record_index && read_it->ref_end_pos > maybe_next_region_start)) {
                        block_itertor = it;
                        bam_itertor = read_it;
                        record_index = true;
                    }
                    if (read_it->bam->core.pos > padding_region_end) {
                        break;
                    }
                    else if (read_it->ref_end_pos > padding_region_start) {
                        bam1_t* bam = m_resource->pool_->alloc_bam_struct_pool();
                        if (!bam_copy1(bam, read_it->bam)) {
                            return;
                        }
                        m_resource->pool_->peek_bam(bam);
                        reads.push_back(bam);
                    }
                }
                it++;
            }
        }
        max_used_region = std::max(m_resource->pool_->get_used_data_size() - current_used_size, max_used_region);
        result.emplace_back(reads, region, region_id);
        reads.clear();
        region_id++;
    }

    if ((force && result.size() > 0) || is_ref_force || m_resource->pool_->get_unsed_data_size() <= max_used_region * EXTRA_FACTOR ||
        (!result.empty() && !engine->check_first_is_same_chrom(result.front().region->tid))) {  // 符合输出条件
        int target_len = 0;
        std::shared_ptr<char> ref_base = m_fasta_loader->get(result.front().region->tid, target_len);
        std::shared_ptr<RegionSource> source =
            std::make_shared<RegionSource>(ref_base, result.front().region->tid, target_len, result, source_id, false, m_resource);

        if (db_) {
            calculate_db(source);
        }

        m_region_queue->push(source);
        result.clear();
        m_resource = nullptr;
        source_id++;
        max_used_region = 0;
    }
}

void RegionManager::flush(bool force)
{
    if (force) {
        int interval_tid = engine->get_interval_tid();

        if (interval_tid == -1) return;
        // RovacaLogger::info("here {}-{}", last_tid, interval_tid);
        if (last_tid != interval_tid) poll(true);
        for (; last_tid < interval_tid; last_tid++) {
            RovacaLogger::info("flush region tid {} finished", last_tid);
            m_fasta_loader->pop();
        }
    }
}

void RegionManager::update_finish_sig()
{
    p_hc_region_active_storage region = nullptr;
    last_tid++;
    while (last_tid < static_cast<int32_t>(m_fasta_loader->get_contig().key.size()) && force_non_active_) {
        int target_len = 0;
        std::shared_ptr<char> ref_base = m_fasta_loader->get(last_tid, target_len);
        if (!m_resource) m_resource = m_region_resource->pop();

        if (!engine->isWES()) {
            region = (p_hc_region_active_storage)malloc(sizeof(hc_region_active_storage));
            region->tid = last_tid;
            region->active = 0;
            region->start_index = 0;
            region->end_index = target_len;
            result.emplace_back(reads, region, region_id);
            std::shared_ptr<RegionSource> source =
                std::make_shared<RegionSource>(ref_base, last_tid, target_len, result, source_id, false, m_resource);

            if (db_) {
                calculate_db(source);
            }

            m_region_queue->push(source);
            region_id++;
            source_id++;
            m_resource = nullptr;
        }
        else {
            engine->flush_region_WES(last_tid);
            if (!engine->region_is_empty()) {
                while ((region = engine->poll()) != nullptr) {
                    /* code */
                    result.emplace_back(reads, region, region_id);
                    region_id++;
                }
                std::shared_ptr<RegionSource> source =
                    std::make_shared<RegionSource>(ref_base, last_tid, target_len, result, source_id, false, m_resource);

                if (db_) {
                    calculate_db(source);
                }

                m_region_queue->push(source);
                source_id++;
                m_resource = nullptr;
            }
        }
        m_fasta_loader->pop();
        last_tid++;
        result.clear();
    }
    std::shared_ptr<RegionSource> source = std::make_shared<RegionSource>(nullptr, -1, -1, result, source_id, true, nullptr);
    m_region_queue->push(source);

    if (db_) {
        db_->update_finish();
    }
}

void RegionManager::calculate_db(std::shared_ptr<RegionSource>& source)
{
    if (__glibc_unlikely(db_tid_ != source->tid_)) {
        db_offset_ = 0;
        db_->notify(db_tid_);
        db_tid_ = source->tid_;
    }

    // 此函数可能阻塞
    source->db_data_ = db_->get(source->tid_);

    int32_t tid_data_len = static_cast<int32_t>(source->db_data_->size());
    for (; db_offset_ < tid_data_len; ++db_offset_) {
        if (source->db_data_->at(db_offset_)->pos >= result.front().region->start_index) {
            source->db_offset_ = db_offset_;
            break;
        }
    }
}
