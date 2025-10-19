#ifndef __ASSEMBLE_READS_H
#define __ASSEMBLE_READS_H

#include "htslib/sam.h"
#include "rbtree_shared.h"

typedef struct hc_apply_one_read_t
{
    bam1_t read;
    unsigned char* read_data;
    char* cigar_str;
    hts_pos_t pos_start;
    hts_pos_t pos_end;
    uint32_t cigar_len;
    uint32_t read_len;
    uint32_t ref_len;
    hts_pos_t insert_size;

    uint32_t* cigar;
    uint32_t* cigar_back;
    uint8_t* seq;
    uint8_t* qual;

    // cache for modifying reads.
    uint32_t* cigar_cache;
    uint32_t* cigar_cache_back;
    uint8_t* seq_cache;
    struct hc_apply_one_read_t* prev;
    struct hc_apply_one_read_t* next;
    // rb tree sort
    struct rb_node node;
    // Mem resource flag of cigar_str/cigar_cache/cigar_cache_back
    // 0 for malloc, 1 for user defined pool.
    uint32_t mempolicy;
} hc_apply_one_read, *p_hc_apply_one_read;

typedef struct hc_region_active_storage_t
{
    int tid;
    uint32_t active;

    hts_pos_t start_index;
    hts_pos_t end_index;

    struct ActiveSpan
    {
        hts_pos_t start;
        hts_pos_t end;
    } activeSpan;

    struct PaddedSpan
    {
        hts_pos_t start;
        hts_pos_t end;
    } paddedSpan;

} hc_region_active_storage, *p_hc_region_active_storage;

#endif
