#include "avx_impl.h"
#include "avx_type.h"
#include "common.h"
#include "const.h"
#include "context.h"

using MainType = float;   // 精度类型
using SimdType = __m512;  // 并行向量类型
using MaskType = uint32_t;
using UnionType = mix_F512;
using BitMaskVec = __m512i;
constexpr MaskType s_mask_bit_src = 0x1;
constexpr uint32_t s_simd_data_len = 512;  // 并行向量长度
constexpr uint32_t s_precision_len = sizeof(MainType);
constexpr MaskType s_mask_all_ones = s_mask_all_ones_s;
constexpr uint32_t s_main_type_size = s_precision_len * s_bits_per_byte;
constexpr uint32_t s_avx_length = s_simd_data_len / s_main_type_size;

using ProbalityArrays = ProbalityPtrDataAvx512S;

#define vec_set1_val_i(v) _mm512_set1_epi32(v)
#define vec_set1_val_f(v) _mm512_set1_ps(v)
#define vec_add(v1, v2)   _mm512_add_ps(v1, v2)
#define vec_sub(v1, v2)   _mm512_sub_ps(v1, v2)
#define vec_mul(v1, v2)   _mm512_mul_ps(v1, v2)
#define vec_div(v1, v2)   _mm512_div_ps(v1, v2)
#define vec_set_high(v)   _mm512_set_ps(v, zero, zero, zero, zero, zero, zero, zero, zero, zero, zero, zero, zero, zero, zero, zero)

static void precompute_masks_avx512_s(const TestCase &tc, uint32_t cols, uint32_t num_mask_vecs, MaskType (*mask_arr)[s_num_distinct_chars])
{
    constexpr uint32_t mask_bit_cnt = s_main_type_size;

    for (uint32_t vi{0}; vi < num_mask_vecs; ++vi) {
        for (uint32_t rs{0}; rs < s_num_distinct_chars; ++rs) {
            mask_arr[vi][rs] = 0;
        }
        mask_arr[vi][s_ambig_char] = s_mask_all_ones;
    }

    for (uint32_t col{1}; col < cols; ++col) {
        uint32_t index = (col - 1) / mask_bit_cnt;
        uint32_t offset = (col - 1) % mask_bit_cnt;
        uint8_t hap_char = ConvertChar::get(tc.hap[col - 1]);
        MaskType bit_mask = s_mask_bit_src << (mask_bit_cnt - 1 - offset);
        mask_arr[index][hap_char] |= bit_mask;

        if (unlikely(hap_char == s_ambig_char)) {
            for (uint32_t ci{0}; ci < s_num_distinct_chars; ++ci) {
                mask_arr[index][ci] |= bit_mask;
            }
        }
    }
}

static void initialize_vectors_avx512_s(uint32_t rows, uint32_t cols, MainType *shift_out_m, MainType *shift_out_x, MainType *shift_out_y,
                                        Context<MainType> &ctx, const TestCase &tc, const ProbalityArrays &pa, SimdType *distm1d)
{
    MainType zero = Context<MainType>::_(0.0);
    MainType init_y = Context<MainType>::INITIAL_CONSTANT / static_cast<MainType>(tc.haplen);
    for (uint32_t s{0}; s < rows + cols + s_avx_length; s++) {
        shift_out_m[s] = zero;
        shift_out_x[s] = zero;
        shift_out_y[s] = init_y;
    }

    auto *ptr_p_mm = reinterpret_cast<MainType *>(pa.p_mm_arr);
    auto *ptr_p_gapm = reinterpret_cast<MainType *>(pa.p_gapm_arr);
    auto *ptr_p_mx = reinterpret_cast<MainType *>(pa.p_mx_arr);
    auto *ptr_p_xx = reinterpret_cast<MainType *>(pa.p_xx_arr);
    auto *ptr_p_my = reinterpret_cast<MainType *>(pa.p_my_arr);
    auto *ptr_p_yy = reinterpret_cast<MainType *>(pa.p_yy_arr);
    auto *ptr_distm1d = reinterpret_cast<MainType *>(distm1d);

    *ptr_p_mm = Context<MainType>::_(0.0);
    *ptr_p_xx = Context<MainType>::_(0.0);
    *ptr_p_yy = Context<MainType>::_(0.0);
    *ptr_p_mx = Context<MainType>::_(0.0);
    *ptr_p_my = Context<MainType>::_(0.0);
    *ptr_p_gapm = Context<MainType>::_(0.0);

    uint32_t remaining_rows = (rows - 1) % s_avx_length;
    uint32_t stripe_cnt = (rows - 1) / s_avx_length;

    uint8_t _i, _d, _c, _q;
    uint32_t base_idx, stripe_idx = 0;
    uint32_t head_offset, tail_offset;
    for (; stripe_idx < stripe_cnt; ++stripe_idx) {
        base_idx = stripe_idx * s_avx_length;
        for (head_offset = 0, tail_offset = s_avx_length - 1; head_offset < s_avx_length; ++head_offset, --tail_offset) {
            _i = tc.i[base_idx + head_offset] & 127;
            _d = tc.d[base_idx + head_offset] & 127;
            _c = tc.c[base_idx + head_offset] & 127;
            _q = tc.q[base_idx + head_offset] & 127;

            *(ptr_p_mm + base_idx + tail_offset) = ctx.set_mm_prob(_i, _d);
            *(ptr_p_gapm + base_idx + tail_offset) = Context<MainType>::_(1.0) - Context<MainType>::ph2pr[_c];
            *(ptr_p_mx + base_idx + tail_offset) = Context<MainType>::ph2pr[_i];
            *(ptr_p_xx + base_idx + tail_offset) = Context<MainType>::ph2pr[_c];
            *(ptr_p_my + base_idx + tail_offset) = Context<MainType>::ph2pr[_d];
            *(ptr_p_yy + base_idx + tail_offset) = Context<MainType>::ph2pr[_c];
            *(ptr_distm1d + base_idx + tail_offset) = Context<MainType>::ph2pr[_q];
        }
    }

    // 不足一个 s_avx_length 的 stripe 初始化
    base_idx = stripe_idx * s_avx_length;
    for (head_offset = 0, tail_offset = s_avx_length - 1; head_offset < remaining_rows; ++head_offset, --tail_offset) {
        _i = tc.i[base_idx + head_offset] & 127;
        _d = tc.d[base_idx + head_offset] & 127;
        _c = tc.c[base_idx + head_offset] & 127;
        _q = tc.q[base_idx + head_offset] & 127;

        *(ptr_p_mm + base_idx + tail_offset) = ctx.set_mm_prob(_i, _d);
        *(ptr_p_gapm + base_idx + tail_offset) = Context<MainType>::_(1.0) - Context<MainType>::ph2pr[_c];
        *(ptr_p_mx + base_idx + tail_offset) = Context<MainType>::ph2pr[_i];
        *(ptr_p_xx + base_idx + tail_offset) = Context<MainType>::ph2pr[_c];
        *(ptr_p_my + base_idx + tail_offset) = Context<MainType>::ph2pr[_d];
        *(ptr_p_yy + base_idx + tail_offset) = Context<MainType>::ph2pr[_c];
        *(ptr_distm1d + base_idx + tail_offset) = Context<MainType>::ph2pr[_q];
    }
}

static void stripe_initialization_avx512_s(uint32_t stripe_idx, SimdType &p_gapm, SimdType &p_mm, SimdType &p_mx, SimdType &p_xx,
                                           SimdType &p_my, SimdType &p_yy, SimdType &distm, SimdType &_1_distm, const SimdType *distm1d,
                                           const ProbalityArrays &pa, SimdType &m_t_1, SimdType &m_t_2, SimdType &x_t_1, SimdType &x_t_2,
                                           SimdType &y_t_1, SimdType &y_t_2, SimdType &m_t_1_y, const MainType *shift_out_x,
                                           const MainType *shift_out_m, const TestCase &tc)
{
    p_mm = pa.p_mm_arr[stripe_idx];
    p_gapm = pa.p_gapm_arr[stripe_idx];
    p_mx = pa.p_mx_arr[stripe_idx];
    p_xx = pa.p_xx_arr[stripe_idx];
    p_my = pa.p_my_arr[stripe_idx];
    p_yy = pa.p_yy_arr[stripe_idx];

    const MainType zero = Context<MainType>::_(0.0);
    const MainType init_y = Context<MainType>::INITIAL_CONSTANT / static_cast<MainType>(tc.haplen);
    const SimdType packed1 = vec_set1_val_f(1.0);
    const SimdType packed3 = vec_set1_val_f(3.0);

    distm = distm1d[stripe_idx];
    _1_distm = vec_sub(packed1, distm);
    distm = vec_div(distm, packed3);

    /* initialize m_t_2, m_t_1, x_t_2, x_t_1, y_t_2, y_t_1 */
    m_t_2 = vec_set1_val_f(zero);
    x_t_2 = vec_set1_val_f(zero);

    if (stripe_idx == 0) {
        m_t_1 = vec_set1_val_f(zero);
        x_t_1 = vec_set1_val_f(zero);
        y_t_2 = vec_set_high(init_y);
        y_t_1 = vec_set1_val_f(zero);
    }
    else {
        x_t_1 = vec_set_high(shift_out_x[s_avx_length]);
        m_t_1 = vec_set_high(shift_out_m[s_avx_length]);
        y_t_2 = vec_set1_val_f(zero);
        y_t_1 = vec_set1_val_f(zero);
    }

    m_t_1_y = m_t_1;
}

static inline void init_masks_for_row_avx512_s(const TestCase &tc, uint8_t *rs_arr, BitMaskVec &last_mask_shift_out,
                                               uint32_t begin_row_index, uint32_t num_rows_to_process)
{
    for (uint32_t ri{0}; ri < num_rows_to_process; ++ri) {
        rs_arr[ri] = ConvertChar::get(tc.rs[ri + begin_row_index - 1]);
    }

    last_mask_shift_out = vec_set1_val_i(0);
}

// 每个 cc 文件单独修改
static void update_masks_for_cols_avx512_s(uint32_t mask_index, BitMaskVec &bit_mask_vec, MaskType (*mask_arr)[s_num_distinct_chars],
                                           const uint8_t *rs_arr, BitMaskVec &last_mask_shift_out)
{
    MaskType *arr = mask_arr[mask_index];

    auto val0 = static_cast<int32_t>(arr[rs_arr[s_rbi_0]]);
    auto val1 = static_cast<int32_t>(arr[rs_arr[s_rbi_1]]);
    auto val2 = static_cast<int32_t>(arr[rs_arr[s_rbi_2]]);
    auto val3 = static_cast<int32_t>(arr[rs_arr[s_rbi_3]]);
    auto val4 = static_cast<int32_t>(arr[rs_arr[s_rbi_4]]);
    auto val5 = static_cast<int32_t>(arr[rs_arr[s_rbi_5]]);
    auto val6 = static_cast<int32_t>(arr[rs_arr[s_rbi_6]]);
    auto val7 = static_cast<int32_t>(arr[rs_arr[s_rbi_7]]);
    auto val8 = static_cast<int32_t>(arr[rs_arr[s_rbi_8]]);
    auto val9 = static_cast<int32_t>(arr[rs_arr[s_rbi_9]]);
    auto valA = static_cast<int32_t>(arr[rs_arr[s_rbi_A]]);
    auto valB = static_cast<int32_t>(arr[rs_arr[s_rbi_B]]);
    auto valC = static_cast<int32_t>(arr[rs_arr[s_rbi_C]]);
    auto valD = static_cast<int32_t>(arr[rs_arr[s_rbi_D]]);
    auto valE = static_cast<int32_t>(arr[rs_arr[s_rbi_E]]);
    auto valF = static_cast<int32_t>(arr[rs_arr[s_rbi_F]]);

    BitMaskVec src = _mm512_set_epi32(val0, val1, val2, val3, val4, val5, val6, val7, val8, val9, valA, valB, valC, valD, valE, valF);
    BitMaskVec right_shift = _mm512_set_epi32(s_mc_0, s_mc_1, s_mc_2, s_mc_3, s_mc_4, s_mc_5, s_mc_6, s_mc_7, s_mc_8, s_mc_9, s_mc_A,
                                              s_mc_B, s_mc_C, s_mc_D, s_mc_E, s_mc_F);
    BitMaskVec mask_vec = _mm512_srlv_epi32(src, right_shift);
    bit_mask_vec = _mm512_or_si512(mask_vec, last_mask_shift_out);

    BitMaskVec reserved_mask = _mm512_set_epi32(s_rm_0, s_rm_1, s_rm_2, s_rm_3, s_rm_4, s_rm_5, s_rm_6, s_rm_7, s_rm_8, s_rm_9, s_rm_10,
                                                s_rm_11, s_rm_12, s_rm_13, s_rm_14, s_rm_15);
    BitMaskVec reserved_src = _mm512_and_si512(src, reserved_mask);
    BitMaskVec left_shift = _mm512_set_epi32(s_sc_20, s_sc_1F, s_sc_1E, s_sc_1D, s_sc_1C, s_sc_1B, s_sc_1A, s_sc_19, s_sc_18, s_sc_17,
                                             s_sc_16, s_sc_15, s_sc_14, s_sc_13, s_sc_12, s_sc_11);
    last_mask_shift_out = _mm512_sllv_epi32(reserved_src, left_shift);
}

// // 每个 cc 文件单独修改
static inline void compute_dist_vec_avx512_s(BitMaskVec &bit_mask_vec, SimdType &distm_chosen, const SimdType &distm,
                                             const SimdType &_1_distm)
{
    // mask 最高位
    __m512i val = _mm512_set1_epi32(static_cast<int32_t>(0x80000000));
    __mmask16 mask = _mm512_test_epi32_mask(val, bit_mask_vec);
    distm_chosen = _mm512_mask_blend_ps(mask, distm, _1_distm);
    bit_mask_vec = _mm512_slli_epi32(bit_mask_vec, 0x1);
}

static inline void compute_mxy_avx512_s(SimdType &m_t, SimdType &m_t_y, SimdType &x_t, SimdType &y_t, const SimdType &m_t_2,
                                        const SimdType &x_t_2, const SimdType &y_t_2, const SimdType &m_t_1, const SimdType &x_t_1,
                                        const SimdType &m_t_1_y, const SimdType &y_t_1, const SimdType &p_mm, const SimdType &p_gapm,
                                        const SimdType &p_mx, const SimdType &p_xx, const SimdType &p_my, const SimdType &p_yy,
                                        const SimdType &distm_sel)
{
    /* compute m_t <= distm * (p_mm*m_t_2 + p_gapm*x_t_2 + p_gapm*y_t_2) */
    m_t = vec_mul(vec_add(vec_add(vec_mul(m_t_2, p_mm), vec_mul(x_t_2, p_gapm)), vec_mul(y_t_2, p_gapm)), distm_sel);

    m_t_y = m_t;

    /* compute x_t */
    x_t = vec_add(vec_mul(m_t_1, p_mx), vec_mul(x_t_1, p_xx));

    /* compute y_t */
    y_t = vec_add(vec_mul(m_t_1_y, p_my), vec_mul(y_t_1, p_yy));
}

// 每个 cc 文件单独修改
static inline void vector_shift_avx512_s(SimdType &x, MainType shift_in, MainType &shift_out)
{
    shift_out = _mm512_cvtss_f32(x);

    IF_32 shiftInH;
    shiftInH.f = shift_in;
    x = _mm512_castsi512_ps(_mm512_alignr_epi32(_mm512_set1_epi32(shiftInH.i), _mm512_castps_si512(x), 0x1));
}

// 每个 cc 文件单独修改
static inline void vector_shift_last_avx512_s(SimdType &x, MainType shift_in)
{
    IF_32 shiftInH;
    shiftInH.f = shift_in;
    x = _mm512_castsi512_ps(_mm512_alignr_epi32(_mm512_set1_epi32(shiftInH.i), _mm512_castps_si512(x), 0x1));
}

float compute_fp_avx512_s(const TestCase &tc)
{
    // 此函数应该放在上层，仅调用一次
    ConvertChar::init();

    uint32_t rows = tc.rslen + 1;
    uint32_t cols = tc.haplen + 1;
    uint32_t mavx_count = (rows + s_avx_length - 1) / s_avx_length;

    /* probaility arrays */
    SimdType p_mm_arr[mavx_count], p_gapm_arr[mavx_count], p_mx_arr[mavx_count];
    SimdType p_xx_arr[mavx_count], p_my_arr[mavx_count], p_yy_arr[mavx_count];
    ProbalityArrays pa{p_mm_arr, p_gapm_arr, p_mx_arr, p_xx_arr, p_my_arr, p_yy_arr};

    /* for distm precomputation */
    SimdType distm1d_arr[mavx_count];

    /* carries the values from each stripe to the next stripe */
    MainType shift_out_m[rows + cols + s_avx_length], shift_out_x[rows + cols + s_avx_length], shift_out_y[rows + cols + s_avx_length];

    /* the vector to keep the anti-diagonals of m, x, y*/
    /* current: m_t, x_t, y_t */
    /* previous: m_t_1, x_t_1, y_t_1 */
    /* previous to previous: m_t_2, x_t_2, y_t_2 */
    SimdType m_t, m_t_1, m_t_2, x_t, x_t_1, x_t_2, y_t, y_t_1, y_t_2, m_t_y, m_t_1_y;

    /* probality vectors */
    SimdType p_gapm, p_mm, p_mx, p_xx, p_my, p_yy;
    SimdType distm, _1_distm, distm_chosen;

    MainType result_avx2;
    Context<MainType> ctx;

    MainType zero = Context<MainType>::_(0.0);
    uint32_t remaining_rows = (rows - 1) % s_avx_length;
    uint32_t stripe_cnt = ((rows - 1) / s_avx_length) + (remaining_rows != 0);

    constexpr uint32_t mask_bit_cnt = s_main_type_size;
    const uint32_t num_mask_vecs = (cols + rows + mask_bit_cnt - 1) / mask_bit_cnt;

    /* mask precomputation for distm*/
    MaskType mask_arr[num_mask_vecs][s_num_distinct_chars];
    precompute_masks_avx512_s(tc, cols, num_mask_vecs, mask_arr);

    uint8_t rs_arr[s_avx_length];
    BitMaskVec bit_mask_vec;
    BitMaskVec last_mask_shift_out;
    uint32_t shift_idx;
    uint32_t num_mask_bits_to_process;

    /* precompute initialization for probabilities and shift vector*/
    initialize_vectors_avx512_s(rows, cols, shift_out_m, shift_out_x, shift_out_y, ctx, tc, pa, distm1d_arr);

    for (uint32_t i{0}; i < stripe_cnt - 1; ++i) {
        stripe_initialization_avx512_s(i, p_gapm, p_mm, p_mx, p_xx, p_my, p_yy, distm, _1_distm, distm1d_arr, pa, m_t_1, m_t_2, x_t_1,
                                       x_t_2, y_t_1, y_t_2, m_t_1_y, shift_out_x, shift_out_m, tc);

        init_masks_for_row_avx512_s(tc, rs_arr, last_mask_shift_out, i * s_avx_length + 1, s_avx_length);

        for (uint32_t begin_d{1}; begin_d < cols + s_avx_length; begin_d += s_main_type_size) {
            num_mask_bits_to_process = std::min(s_main_type_size, cols + s_avx_length - begin_d);

            update_masks_for_cols_avx512_s((begin_d - 1) / s_main_type_size, bit_mask_vec, mask_arr, rs_arr, last_mask_shift_out);

            for (uint32_t mbi{0}; mbi < num_mask_bits_to_process; ++mbi) {
                shift_idx = begin_d + mbi + s_avx_length;

                compute_dist_vec_avx512_s(bit_mask_vec, distm_chosen, distm, _1_distm);

                compute_mxy_avx512_s(m_t, m_t_y, x_t, y_t, m_t_2, x_t_2, y_t_2, m_t_1, x_t_1, m_t_1_y, y_t_1, p_mm, p_gapm, p_mx, p_xx,
                                     p_my, p_yy, distm_chosen);

                vector_shift_avx512_s(m_t, shift_out_m[shift_idx], shift_out_m[begin_d + mbi]);
                vector_shift_avx512_s(x_t, shift_out_x[shift_idx], shift_out_x[begin_d + mbi]);
                vector_shift_avx512_s(y_t_1, shift_out_y[shift_idx], shift_out_y[begin_d + mbi]);

                m_t_2 = m_t_1;
                m_t_1 = m_t;
                x_t_2 = x_t_1;
                x_t_1 = x_t;
                y_t_2 = y_t_1;
                y_t_1 = y_t;
                m_t_1_y = m_t_y;
            }
        }
    }

    {
        uint32_t i = stripe_cnt - 1;

        stripe_initialization_avx512_s(i, p_gapm, p_mm, p_mx, p_xx, p_my, p_yy, distm, _1_distm, distm1d_arr, pa, m_t_1, m_t_2, x_t_1,
                                       x_t_2, y_t_1, y_t_2, m_t_1_y, shift_out_x, shift_out_m, tc);

        if (remaining_rows == 0) {
            remaining_rows = s_avx_length;
        }
        init_masks_for_row_avx512_s(tc, rs_arr, last_mask_shift_out, i * s_avx_length + 1, remaining_rows);

        SimdType sum_m = vec_set1_val_f(zero);
        SimdType sum_x = vec_set1_val_f(zero);

        for (uint32_t begin_d{1}; begin_d < cols + remaining_rows - 1; begin_d += s_main_type_size) {
            num_mask_bits_to_process = std::min(s_main_type_size, cols + remaining_rows - 1 - begin_d);

            update_masks_for_cols_avx512_s((begin_d - 1) / s_main_type_size, bit_mask_vec, mask_arr, rs_arr, last_mask_shift_out);

            for (uint32_t mbi{0}; mbi < num_mask_bits_to_process; ++mbi) {
                shift_idx = begin_d + mbi + s_avx_length;

                compute_dist_vec_avx512_s(bit_mask_vec, distm_chosen, distm, _1_distm);

                compute_mxy_avx512_s(m_t, m_t_y, x_t, y_t, m_t_2, x_t_2, y_t_2, m_t_1, x_t_1, m_t_1_y, y_t_1, p_mm, p_gapm, p_mx, p_xx,
                                     p_my, p_yy, distm_chosen);

                sum_m = vec_add(sum_m, m_t);
                vector_shift_last_avx512_s(m_t, shift_out_m[shift_idx]);

                sum_x = vec_add(sum_x, x_t);
                vector_shift_last_avx512_s(x_t, shift_out_x[shift_idx]);

                vector_shift_last_avx512_s(y_t_1, shift_out_y[shift_idx]);

                m_t_2 = m_t_1;
                m_t_1 = m_t;
                x_t_2 = x_t_1;
                x_t_1 = x_t;
                y_t_2 = y_t_1;
                y_t_1 = y_t;
                m_t_1_y = m_t_y;
            }
        }

        UnionType sum_mx;
        sum_mx.d = vec_add(sum_m, sum_x);
        result_avx2 = sum_mx.f[s_avx_length - remaining_rows];
    }

    return result_avx2;
}
