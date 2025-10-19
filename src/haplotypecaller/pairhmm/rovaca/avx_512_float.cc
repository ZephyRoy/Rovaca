#include "avx_512_float.h"

#include <x86intrin.h>

#include <algorithm>
#include <cassert>
#include <cmath>
#include <memory_resource>

#include "context.h"
#include "internal.h"

namespace rovaca
{

using MainType = float;
constexpr uint32_t k_bit_count_per_byte = 8;
constexpr uint32_t k_main_type_size = k_bit_count_per_byte * sizeof(MainType);
constexpr uint32_t k_avx_length = k_avx512_bit_count / k_main_type_size;

static inline __mmask16 generate_read_len_mask(uint32_t read_idx, const uint32_t* lens)
{
    __m512i idx = _mm512_set1_epi32(read_idx);
    __m512i len = _mm512_set_epi32(lens[15], lens[14], lens[13], lens[12], lens[11], lens[10], lens[9], lens[8], lens[7], lens[6], lens[5],
                                   lens[4], lens[3], lens[2], lens[1], lens[0]);
    return _mm512_cmp_epu32_mask(len, idx, _MM_CMPINT_NLE);
}

void compute_full_prob_avx512_float_1xn(const TestCase& tc, float* result)
{
    Context<float> ctx{};
    float zero = ctx._(0.0);
    float init_d = rovaca::Context<float>::INITIAL_CONSTANT / tc.max_hap_len;
    __m512 mm[tc.max_hap_len + 1], ii[tc.max_hap_len + 1], dd[tc.max_hap_len + 1];
    mm[0] = _mm512_set1_ps(zero);
    ii[0] = _mm512_set1_ps(zero);
    dd[0] = _mm512_set1_ps(init_d);
    for (uint32_t i = 1; i <= tc.max_hap_len; ++i) {
        mm[i] = mm[0];
        ii[i] = ii[0];
        dd[i] = dd[0];
    }

    __m512i r, h;
    __m512 M_t, I_t, D_t;
    __mmask16 mask, len_mask;
    __m512 distm, _1_distm, distm_chosen;
    __m512 p_gapm, p_mm, p_mx, p_xx, p_my, p_yy;
    __m512 M, M_i1, M_j1, M_i1j1, I, I_i1, I_j1, I_i1j1, D, D_i1, D_j1, D_i1j1;

    uint32_t hap_base_idx, read_base_idx, read_arr_offset;

    // Process reads with common length
    for (read_base_idx = 0; read_base_idx < tc.min_read_len; ++read_base_idx) {
        read_arr_offset = read_base_idx * k_avx_length;
        r = _mm512_load_epi32(tc.reads + read_arr_offset);
        p_mm = _mm512_load_ps(tc.mm + read_arr_offset);
        p_mx = _mm512_load_ps(tc.mi + read_arr_offset);
        p_xx = _mm512_load_ps(tc.ii + read_arr_offset);
        p_my = _mm512_load_ps(tc.md + read_arr_offset);
        p_yy = _mm512_load_ps(tc.dd + read_arr_offset);
        p_gapm = _mm512_load_ps(tc.gapm + read_arr_offset);
        distm = _mm512_load_ps(tc.distm + read_arr_offset);
        _1_distm = _mm512_load_ps(tc._1_distm + read_arr_offset);

        M_j1 = I_j1 = D_j1 = M_i1j1 = I_i1j1 = _mm512_setzero_ps();
        D_i1j1 = read_base_idx == 0 ? _mm512_set1_ps(init_d) : _mm512_set1_ps(zero);

        M_i1 = mm[0];
        I_i1 = ii[0];
        D_i1 = dd[0];

        for (hap_base_idx = 0; hap_base_idx < tc.min_hap_len; ++hap_base_idx) {
            h = _mm512_set1_epi32(tc.haps[hap_base_idx]);

            // mem pre fetch
            // _mm_prefetch(mm + hap_base_idx + 1, _MM_HINT_T0);
            // _mm_prefetch(ii + hap_base_idx + 1, _MM_HINT_T0);
            // _mm_prefetch(dd + hap_base_idx + 1, _MM_HINT_T0);

            mask = _mm512_test_epi32_mask(r, h);
            distm_chosen = _mm512_mask_blend_ps(mask, distm, _1_distm);

            M = _mm512_mul_ps(
                _mm512_add_ps(_mm512_add_ps(_mm512_mul_ps(M_i1j1, p_mm), _mm512_mul_ps(I_i1j1, p_gapm)), _mm512_mul_ps(D_i1j1, p_gapm)),
                distm_chosen);
            I = _mm512_add_ps(_mm512_mul_ps(M_i1, p_mx), _mm512_mul_ps(I_i1, p_xx));
            D = _mm512_add_ps(_mm512_mul_ps(M_j1, p_my), _mm512_mul_ps(D_j1, p_yy));

            M_i1j1 = M_i1;
            I_i1j1 = I_i1;
            D_i1j1 = D_i1;

            M_j1 = M;
            I_j1 = I;
            D_j1 = D;
            mm[hap_base_idx] = M;
            ii[hap_base_idx] = I;
            dd[hap_base_idx] = D;

            M_i1 = mm[hap_base_idx + 1];
            I_i1 = ii[hap_base_idx + 1];
            D_i1 = dd[hap_base_idx + 1];
        }
    }

    // Use mask to process reads with extra length
    for (; read_base_idx < tc.max_read_len; ++read_base_idx) {
        read_arr_offset = read_base_idx * k_avx_length;
        r = _mm512_load_epi32(tc.reads + read_arr_offset);
        p_mm = _mm512_load_ps(tc.mm + read_arr_offset);
        p_mx = _mm512_load_ps(tc.mi + read_arr_offset);
        p_xx = _mm512_load_ps(tc.ii + read_arr_offset);
        p_my = _mm512_load_ps(tc.md + read_arr_offset);
        p_yy = _mm512_load_ps(tc.dd + read_arr_offset);
        p_gapm = _mm512_load_ps(tc.gapm + read_arr_offset);
        distm = _mm512_load_ps(tc.distm + read_arr_offset);
        _1_distm = _mm512_load_ps(tc._1_distm + read_arr_offset);

        len_mask = generate_read_len_mask(read_base_idx, tc.read_len_arr);

        M_j1 = I_j1 = D_j1 = M_i1j1 = I_i1j1 = _mm512_setzero_ps();
        D_i1j1 = read_base_idx == 0 ? _mm512_set1_ps(init_d) : _mm512_set1_ps(zero);

        M_i1 = mm[0];
        I_i1 = ii[0];
        D_i1 = dd[0];

        for (hap_base_idx = 0; hap_base_idx < tc.min_hap_len; ++hap_base_idx) {
            h = _mm512_set1_epi32(tc.haps[hap_base_idx]);

            // mem pre fetch
            // _mm_prefetch(mm + hap_base_idx + 1, _MM_HINT_T0);
            // _mm_prefetch(ii + hap_base_idx + 1, _MM_HINT_T0);
            // _mm_prefetch(dd + hap_base_idx + 1, _MM_HINT_T0);

            mask = _mm512_test_epi32_mask(r, h);
            distm_chosen = _mm512_mask_blend_ps(mask, distm, _1_distm);
            M_t = _mm512_mul_ps(
                _mm512_add_ps(_mm512_add_ps(_mm512_mul_ps(M_i1j1, p_mm), _mm512_mul_ps(I_i1j1, p_gapm)), _mm512_mul_ps(D_i1j1, p_gapm)),
                distm_chosen);
            I_t = _mm512_add_ps(_mm512_mul_ps(M_i1, p_mx), _mm512_mul_ps(I_i1, p_xx));
            D_t = _mm512_add_ps(_mm512_mul_ps(M_j1, p_my), _mm512_mul_ps(D_j1, p_yy));

            M = _mm512_mask_blend_ps(len_mask, M_i1, M_t);
            I = _mm512_mask_blend_ps(len_mask, I_i1, I_t);
            D = _mm512_mask_blend_ps(len_mask, D_i1, D_t);

            M_i1j1 = M_i1;
            I_i1j1 = I_i1;
            D_i1j1 = D_i1;

            M_j1 = M;
            I_j1 = I;
            D_j1 = D;
            mm[hap_base_idx] = M;
            ii[hap_base_idx] = I;
            dd[hap_base_idx] = D;

            M_i1 = mm[hap_base_idx + 1];
            I_i1 = ii[hap_base_idx + 1];
            D_i1 = dd[hap_base_idx + 1];
        }
    }

    __m512 sum_m = _mm512_setzero_ps();
    __m512 sum_i = _mm512_setzero_ps();
    for (hap_base_idx = 0; hap_base_idx < tc.min_hap_len; ++hap_base_idx) {
        sum_m = _mm512_add_ps(sum_m, mm[hap_base_idx]);
        sum_i = _mm512_add_ps(sum_i, ii[hap_base_idx]);
    }

    float m_result_temp[k_avx_length] __attribute__((aligned(64)));
    float i_result_temp[k_avx_length] __attribute__((aligned(64)));
    _mm512_store_ps(m_result_temp, sum_m);
    _mm512_store_ps(i_result_temp, sum_i);

    for (uint32_t i{0}; i < k_avx_length; ++i) {
        result[i] = m_result_temp[i] + i_result_temp[i];
    }
}

}  // namespace rovaca
