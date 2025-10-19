#include "bam_loader.h"

#include <algorithm>
#include <iostream>
#include <tuple>

#include "../rovaca_logger/rovaca_logger.h"

BamLoader::BamLoader(const std::vector<const char*>& bam_list, const htsThreadPool* pool, bool force_bed)
    : nfiles(bam_list.size())
    , multi_reg(force_bed)
    , all_finished(false)
    , go_on_pop(false)
    , specified_target(false)
    , bam_fp(nfiles, nullptr)
    , bam_itr(nfiles, nullptr)
    , bam_header(nfiles, nullptr)
    , bam_index(nfiles, nullptr)
    , bam_paths(nfiles, nullptr)
    , reads_cache({})
    , thread_pool(pool)
    , read_pool{BAM1POOL_READCOUNT, BAM1POOL_READSIZE}
{
    for (int i = 0; i < nfiles; i++) {
        bam_fp[i] = sam_open(bam_list[i], "r");
        if (bam_fp[i] == nullptr) {
            RovacaLogger::error("invalid bam file ptr");
            exit(EXIT_FAILURE);
        }
        bam_header[i] = sam_hdr_read(bam_fp[i]);
        if (bam_header[i] == nullptr) {
            RovacaLogger::error("invalid bam header ptr");
            exit(EXIT_FAILURE);
        }
        bam_index[i] = sam_index_load(bam_fp[i], bam_list[i]);
        if (bam_index[i] == nullptr) {
            RovacaLogger::error("invalid bam index ptr");
            exit(EXIT_FAILURE);
        }
        bam_paths[i] = bam_list[i];
    }
    hts_set_opt(bam_fp[0], HTS_OPT_THREAD_POOL, thread_pool);
    hts_set_opt(bam_fp[0], HTS_OPT_CACHE_SIZE, BAM_HTS_CACHE_SIZE);
}

BamLoader::BamLoader(const std::vector<const char*>& bam_list, bool force_bed)
    : BamLoader(bam_list, nullptr, force_bed)
{}

BamLoader::BamLoader(const char* bam_path, bool force_bed)
    : BamLoader(std::vector<const char*>{bam_path}, nullptr, force_bed)
{}

BamLoader::BamLoader(const char* bam_path, const htsThreadPool* pool, bool force_bed)
    : BamLoader(std::vector<const char*>{bam_path}, pool, force_bed)
{}

BamLoader::~BamLoader()
{
    for (int i = 0; i < nfiles; i++) {
        if (bam_index[i]) hts_idx_destroy(bam_index[i]);
        if (bam_header[i]) bam_hdr_destroy(bam_header[i]);
        if (bam_itr[i]) hts_itr_destroy(bam_itr[i]);
        if (bam_fp[i]) sam_close(bam_fp[i]);
    }
}

void BamLoader::set_target(const p_bed_intervals interval)
{
    if (interval->regarry) {
        char** regarray = const_cast<char**>(interval->regarry);
        all_finished = false;
        multi_reg = true;
        for (int i = 0; i < nfiles; i++) {
            if (bam_itr[i]) hts_itr_destroy(bam_itr[i]);
            bam_itr[i] = sam_itr_regarray(bam_index[i], bam_header[i], regarray, interval->n);
        }
    }
}

void BamLoader::set_target(const char* targets)
{
    all_finished = false;
    specified_target = true;
    for (int i = 0; i < nfiles; i++) {
        if (bam_itr[i]) hts_itr_destroy(bam_itr[i]);
        bam_itr[i] = sam_itr_querys(bam_index[i], bam_header[i], targets);
    }
}

bool BamLoader::has_next() { return !(all_finished && reads_cache.empty()); }

bam1_t* BamLoader::get_next_read(int& min_index)
{
    if (!(all_finished || go_on_pop) || reads_cache.empty()) {
        process_load();
        if (nfiles > 1)
            std::sort(reads_cache.begin(), reads_cache.end(), [](const auto& a, const auto& b) {
                uint64_t pos_a = (uint64_t)((a.first->core.pos + 1) << 1) | bam_is_rev(a.first);
                uint64_t pos_b = (uint64_t)((b.first->core.pos + 1) << 1) | bam_is_rev(b.first);
                return std::tie(a.first->core.tid, pos_a, a.second) < std::tie(b.first->core.tid, pos_b, b.second);
            });
    }
    if (!reads_cache.empty()) {
        go_on_pop = false;
        auto ret = reads_cache.front().first;
        min_index = reads_cache.front().second;
        auto next_itr = std::next(reads_cache.begin());
        if (next_itr != reads_cache.end() && next_itr->first->core.pos == ret->core.pos && next_itr->second == min_index) {
            go_on_pop = true;
        }
        reads_cache.erase(reads_cache.begin());
        return ret;
    }
    return nullptr;
}

void BamLoader::process_load()
{
    for (int i = 0; i < nfiles; i++) {
        all_finished = true;
        if (multi_reg && (bam_itr[i] == NULL || bam_itr[i]->finished)) {
            continue;
        }
        auto read = read_pool.allocate();
        if (multi_reg) {
            all_finished = (hts_itr_multi_next(bam_fp[i], bam_itr[i], read) < 0);
        }
        else if (specified_target) {
            all_finished = (sam_itr_next(bam_fp[i], bam_itr[i], read) < 0);
        }
        else {
            all_finished = (sam_read1(bam_fp[i], bam_header[i], read) < 0);
        }
        if (read != nullptr && !all_finished) {
            reads_cache.emplace_back(std::make_pair(read, i));
        }
        else {
            read_pool.deallocate(read);
        }
    }
}
const std::vector<bam_hdr_t*>& BamLoader::get_sam_hdr() { return bam_header; }
void BamLoader::read_recovery(bam1_t* it) { read_pool.deallocate(it); }