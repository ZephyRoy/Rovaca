#ifndef __SHARED_HC_ACTIVE_H__
#define __SHARED_HC_ACTIVE_H__

#include <stdbool.h>

#include "assemble_interface.h"
#include "hc_assemble_gatk_sw.h"
#include "hc_assemble_graph.h"
#include "hc_assemble_simple_thread_mem.h"
#include "mem_pool_auto_enlarge.h"
#include "rbtree_shared.h"
#include "smithwaterman_common.h"
#include "uthash.h"

#define HC_APPLY_MAX_READ_BUFFER   (1024)
#define HC_APPLY_MAX_READ_LEN      (512)
#define HC_APPLY_READ_PADDED_SPAN  (100)
#define HC_APPLY_REFERENCE_PADDING (500)
#define RING_MAX_BAM1_INIT_LEN     (4096)
#define ASSEMBLE_CIGAR_STRING_MEM  (128)  // 须与AssembleEngine中k_assemble_cigar_string_mem定义一致
#define ASSEMBLE_CIGAR_CACHE_MEM   (128)  // 须与AssembleEngine中k_assemble_cigar_cache_mem定义一致
#define ASSEMBLE_READ_CACHE_MEM    (ASSEMBLE_CIGAR_STRING_MEM + 2 * ASSEMBLE_CIGAR_CACHE_MEM * sizeof(uint32_t))

typedef enum { SOFTCLIP = BAM_CSOFT_CLIP, HARDCLIP = BAM_CHARD_CLIP } clipping_type;

typedef struct AssembleReadsBuffer
{
    uint8_t* buffer_;
    uint32_t used_;
    uint32_t capacity_;
} AssembleReadsBuffer, *pAssembleReadsBuffer;

typedef struct hc_apply_t
{
    p_hc_apply_one_read reads;
    struct rb_root read_sorted;
    uint32_t reads_count;
    p_hc_region_active_storage region;  // active span
    const uint8_t* ref;                 // padded with 500
    const uint8_t* ref_haplotype;       // within padding
    uint32_t span_start;                // padded span(100)
    uint32_t span_end;
    uint32_t ref_start;
    uint32_t ref_end;

    uint32_t ref_chr_len;
    uint32_t ref_len;
    uint32_t ref_hap_len;
    p_assemble_graph read_thread_graph;  // reference path.
    hc_assemble_read_graph read_graph;
    gatk_sw_storage gatk_sw;
    p_lib_sw_avx sw_avx;
    pAssembleReadsBuffer reads_buffer_mem;

} hc_apply, *p_hc_apply;

// 在上层初始化，获取组装模块的指定参数值
typedef struct ApplyArgument_t
{
    struct Kmer
    {
        size_t nkmer;
        uint32_t kmer[0];
    }* kmers;

    bool recover_dangling_branches;
    bool recover_all_dangling_branches;
    bool allow_non_unique_kmer_in_ref;
    bool soft_clip_low_quality_ends;
    bool removePathsNotConnectedToRef;
    bool generateSeqGraph;
    bool pruneBeforeCycleCounting;
    uint8_t minBaseQualityToUseInAssembly;
    uint32_t min_dangling_branch_length;
    uint32_t min_pruning;
    uint32_t max_unpruned_variants;
    uint32_t num_pruning_samples;
    uint32_t max_num_haplotypes_inpopulation;
    uint32_t kmer_length_for_read_error_correction;
    uint32_t min_observations_for_kmer_tobesolid;
    int min_matchingbases_to_dangling_end_recovery;
    double initial_error_rate_for_pruning;
    double pruning_logodds_threshold;
    double pruning_seeding_logodds_threshold;
    double pileup_error_correction_logodds;
    struct Debug
    {
        bool debug_print;
        bool debugGraphTransformations;
    } debug;
} ApplyArgument, *pApplyArgument;

extern pApplyArgument apply_arguments;

#endif  // !__SHARED_HC_ACTIVE_H__