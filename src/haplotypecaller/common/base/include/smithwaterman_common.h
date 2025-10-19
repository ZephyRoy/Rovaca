/**
 * The MIT License (MIT)
 *
 * Copyright (c) 2016-2021 Intel Corporation
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in all
 * copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE.
 */
#ifndef SMITHWATERMAN_COMMON_H
#define SMITHWATERMAN_COMMON_H

#ifdef __cplusplus
extern "C"
{
#endif

#include <immintrin.h>
#include <stdint.h>
#include <string.h>
#include <x86intrin.h>  // SIMD intrinsics for GCC

#define SW_AVX_MAX_SEQ_SUPPORTED (65536)
#define SW_AVX_MAX_DATA_LEN      (65536 * 5)

    enum LIB_SW_AVX_OVERHANG_STRATEGY_STATUS {
        LIB_SW_AVX_OVERHANG_STRATEGY_MATCH = 0,
        LIB_SW_AVX_OVERHANG_STRATEGY_INSERT = 1,
        LIB_SW_AVX_OVERHANG_STRATEGY_DELETE = 2,
        LIB_SW_AVX_OVERHANG_STRATEGY_INSERT_EXT = 4,
        LIB_SW_AVX_OVERHANG_STRATEGY_DELETE_EXT = 8,
        LIB_SW_AVX_OVERHANG_STRATEGY_SOFT_CLIP = 9,
        LIB_SW_AVX_OVERHANG_STRATEGY_INDEL = 10,
        LIB_SW_AVX_OVERHANG_STRATEGY_LEADING_INDEL = 11,
        LIB_SW_AVX_OVERHANG_STRATEGY_IGNORE = 12
    };

    typedef struct sw_parameters_t
    {
        int32_t match;
        int32_t mismatch;
        int32_t gap_open;
        int32_t gap_extend;
    } sw_parameters;

    typedef struct lib_sw_avx_st
    {
        uint8_t* seq1;  // ref
        uint8_t* seq2;  // seq
        int16_t len1;   // ref len
        int16_t len2;   // seq len
        int32_t score;
        int16_t max_i;
        int16_t max_j;
        int16_t* btrack;
        uint32_t* cigar;
        uint32_t cigar_count;
        uint32_t cigar_max_len;
        int16_t alignment_offset;

        int32_t max_seq_len;
        int32_t seq_len_changed;

        struct
        {
            int32_t* e;
            int16_t* back_track;
            int16_t* cigar_buf;
        } mm_aloc;

        sw_parameters para;

        uint32_t data_len;
        uint32_t max_data_len;

        struct
        {
            uint32_t* seq1_rev;
            uint32_t* seq2;
        } cacl_cache;
        void (*avx_function)(struct lib_sw_avx_st* p);

        uint8_t data[];
    } lib_sw_avx, *p_lib_sw_avx;

    p_lib_sw_avx sw_avx_init(uint8_t is_assemble);
    void sw_avx_finit(p_lib_sw_avx ret);
    void sw_avx512_run(p_lib_sw_avx p);
    void sw_avx2_run(p_lib_sw_avx p);

#ifdef __cplusplus
}
#endif
#endif