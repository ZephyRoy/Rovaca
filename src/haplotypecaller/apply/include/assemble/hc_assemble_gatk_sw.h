#ifndef __SHARED_HC_ASSEMBLE_GATK_SW_H__
#define __SHARED_HC_ASSEMBLE_GATK_SW_H__

#include<stdint.h>
#define GATK_SW_CIGAR_MAX_LEN (1024)
#define GATK_SW_SEQ_MAX_LEN   (1024)
#define GATK_SW_DATA_LEN      (GATK_SW_SEQ_MAX_LEN * GATK_SW_SEQ_MAX_LEN * 2 + GATK_SW_SEQ_MAX_LEN * 4)
// match=1, mismatch = -1/3, gap=-(1+k/3)
#define GATK_SW_ORIGINAL_DEFAULT {3, -1, -4, -3};
#define GATK_SW_STANDARD_NGS     {25, -50, -110, -6};

enum GATK_SW_OVERHANG_STRATEGY {
    /*
     * Add softclips for the overhangs
     */
    GATK_SW_OVERHANG_STRATEGY_SOFTCLIP,

    /*
     * Treat the overhangs as proper insertions/deletions
     */
    GATK_SW_OVERHANG_STRATEGY_INDEL,

    /*
     * Treat the overhangs as proper insertions/deletions for leading (but not trailing) overhangs.
     * This is useful e.g. when we want to merge dangling tails in an assembly graph: because we don't
     * expect the dangling tail to reach the end of the reference path we are okay ignoring trailing
     * deletions - but leading indels are still very much relevant.
     */
    GATK_SW_OVERHANG_STRATEGY_LEADING_INDEL,

    /*
     * Just ignore the overhangs
     */
    GATK_SW_OVERHANG_STRATEGY_IGNORE
};

// the set of weights to use to configure the alignment
typedef struct gatk_sw_parameters_t
{
    int32_t matchValue;
    int32_t mismatchPenalty;
    int32_t gapOpenPenalty;
    int32_t gapExtendPenalty;
} gatk_sw_parameters, *p_gatk_sw_parameters;

typedef struct gatk_sw_storage_t
{
    // input
    struct
    {
        uint8_t* ref;      // ref sequence
        uint8_t* alt;      // alt sequence
        uint32_t ref_len;  // ref sequence length
        uint32_t alt_len;  // alt sequence length

        uint32_t overhang_strategy;  // the strategy to use for dealing with overhangs
    } input;

    // paramates
    gatk_sw_parameters paramates;  // the set of weights to use to configure the alignment

    struct
    {
        int32_t* sw;      // the Smith-Waterman matrix to populate
        int32_t* btrack;  // the back track matrix to populate

        // access is pricey if done enough times so we extract those out
        int32_t* best_gap_v;
        int32_t* gap_size_v;
        int32_t* best_gap_h;
        int32_t* gap_size_h;

        uint32_t x_len;
        uint32_t y_len;

        int32_t* current_data;
        int32_t data[GATK_SW_DATA_LEN];
    } cacl_cache;

    // output
    struct
    {
        uint32_t cigar[GATK_SW_CIGAR_MAX_LEN];
        uint32_t cigar_len;

        int alignment_offset;
    } output;
} gatk_sw_storage, *p_gatk_sw_storage;

#endif