#ifndef BAM_LOADER_H
#define BAM_LOADER_H

#include <vector>

#include "bam1_pool.h"
#include "bed_loader.h"
#include "htslib/sam.h"
#include "htslib/thread_pool.h"
class BamLoader
{
private:
#define BAM1POOL_READCOUNT (50)
#define BAM1POOL_READSIZE  (4098)
#define BAM_HTS_CACHE_SIZE (64 * 1024 * 1024)
private:
    int nfiles;
    bool multi_reg;
    bool all_finished;
    bool go_on_pop;
    bool specified_target;
    std::vector<samFile*> bam_fp;
    std::vector<hts_itr_t*> bam_itr;
    std::vector<bam_hdr_t*> bam_header;
    std::vector<hts_idx_t*> bam_index;
    std::vector<const char*> bam_paths;
    std::vector<std::pair<bam1_t*, int>> reads_cache;
    const htsThreadPool* thread_pool;
    Bam1Pool read_pool;

public:
    BamLoader(const char* bam_path, bool force_bed);
    BamLoader(const char* bam_path, const htsThreadPool* pool, bool force_bed);
    BamLoader(const std::vector<const char*>& bam_list, bool force_bed);
    BamLoader(const std::vector<const char*>& bam_list, const htsThreadPool* pool, bool force_bed);
    bool has_next();
    bam1_t* get_next_read(int& min_index);
    ~BamLoader();

    // for wes, load read by regions.
    void set_target(const p_bed_intervals interval);

    // for specified target region.
    void set_target(const char* targets);

    const std::vector<bam_hdr_t*>& get_sam_hdr();
    void read_recovery(bam1_t*);

private:
    void process_load();
};

#endif  // BAM_LOADER_H