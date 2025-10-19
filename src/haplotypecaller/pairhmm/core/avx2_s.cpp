#include "avx_impl.h"
#include "avx_type.h"
#include "common.h"
#include "const.h"
#include "context.h"

using MainType = float;   // 精度类型
using SimdType = __m256;  // 并行向量类型
using MaskType = uint32_t;
using UnionType = mix_F256;
using BitMaskVec = __m256i;
constexpr MaskType s_mask_bit_src = 0x1;
constexpr uint32_t s_simd_data_len = 256;  // 并行向量长度
constexpr uint32_t s_precision_len = sizeof(MainType);
constexpr MaskType s_mask_all_ones = s_mask_all_ones_s;
constexpr uint32_t s_main_type_size = s_precision_len * s_bits_per_byte;
constexpr uint32_t s_avx_length = s_simd_data_len / s_main_type_size;

using ProbalityArrays = ProbalityPtrDataAvx2S;

#define vec_set1_val_i(v) _mm256_set1_epi32(v)
#define vec_set1_val_f(v) _mm256_set1_ps(v)
#define vec_add(v1, v2)   _mm256_add_ps(v1, v2)
#define vec_sub(v1, v2)   _mm256_sub_ps(v1, v2)
#define vec_mul(v1, v2)   _mm256_mul_ps(v1, v2)
#define vec_div(v1, v2)   _mm256_div_ps(v1, v2)
#define vec_set_lse(v)    _mm256_set_ps(zero, zero, zero, zero, zero, zero, zero, v)

static void precompute_masks_avx2_s(const TestCase &tc, uint32_t cols, uint32_t num_mask_vecs, MaskType (*mask_arr)[s_num_distinct_chars])
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

static void initialize_vectors_avx2_s(uint32_t rows, uint32_t cols, MainType *shift_out_m, MainType *shift_out_x, MainType *shift_out_y,
                                      Context<MainType> &ctx, const TestCase &tc, const ProbalityArrays &pa, SimdType *distm1d)
{
    const MainType zero = Context<MainType>::_(0.0);
    const MainType init_y = Context<MainType>::INITIAL_CONSTANT / static_cast<MainType>(tc.haplen);
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

    for (uint32_t r{1}; r < rows; r++) {
        const uint8_t _i = tc.i[r - 1] & 127;
        const uint8_t _d = tc.d[r - 1] & 127;
        const uint8_t _c = tc.c[r - 1] & 127;
        const uint8_t _q = tc.q[r - 1] & 127;

        *(ptr_p_mm + r - 1) = ctx.set_mm_prob(_i, _d);
        *(ptr_p_gapm + r - 1) = Context<MainType>::_(1.0) - Context<MainType>::ph2pr[_c];
        *(ptr_p_mx + r - 1) = Context<MainType>::ph2pr[_i];
        *(ptr_p_xx + r - 1) = Context<MainType>::ph2pr[_c];
        *(ptr_p_my + r - 1) = Context<MainType>::ph2pr[_d];
        *(ptr_p_yy + r - 1) = Context<MainType>::ph2pr[_c];
        *(ptr_distm1d + r - 1) = Context<MainType>::ph2pr[_q];
    }
}

static void stripe_initialization_avx2_s(uint32_t stripe_idx, SimdType &p_gapm, SimdType &p_mm, SimdType &p_mx, SimdType &p_xx,
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
        y_t_2 = vec_set_lse(init_y);
        y_t_1 = vec_set1_val_f(zero);
    }
    else {
        x_t_1 = vec_set_lse(shift_out_x[s_avx_length]);
        m_t_1 = vec_set_lse(shift_out_m[s_avx_length]);
        y_t_2 = vec_set1_val_f(zero);
        y_t_1 = vec_set1_val_f(zero);
    }

    m_t_1_y = m_t_1;
}

static inline void init_masks_for_row_avx2_s(const TestCase &tc, uint8_t *rs_arr, BitMaskVec &last_mask_shift_out, uint32_t begin_row_index,
                                             uint32_t num_rows_to_process)
{
    const uint8_t *dest = tc.rs + begin_row_index - 1;
    _mm_prefetch(dest, _MM_HINT_T0);
    _mm_prefetch(rs_arr, _MM_HINT_T0);
    _mm_prefetch(ConvertChar::k_conversion_table, _MM_HINT_T0);

    last_mask_shift_out = vec_set1_val_i(0);

    if (likely(num_rows_to_process == s_avx_length)) {
        rs_arr[0] = ConvertChar::k_conversion_table[dest[0] - 'A'];
        rs_arr[1] = ConvertChar::k_conversion_table[dest[1] - 'A'];
        rs_arr[2] = ConvertChar::k_conversion_table[dest[2] - 'A'];
        rs_arr[3] = ConvertChar::k_conversion_table[dest[3] - 'A'];
        rs_arr[4] = ConvertChar::k_conversion_table[dest[4] - 'A'];
        rs_arr[5] = ConvertChar::k_conversion_table[dest[5] - 'A'];
        rs_arr[6] = ConvertChar::k_conversion_table[dest[6] - 'A'];
        rs_arr[7] = ConvertChar::k_conversion_table[dest[7] - 'A'];
    }
    else {
        for (uint32_t ri{0}; ri < num_rows_to_process; ++ri) {
            rs_arr[ri] = ConvertChar::get(dest[ri]);
        }
    }
}

// 每个 cc 文件单独修改
static void update_masks_for_cols_avx2_s(uint32_t mask_index, BitMaskVec &bit_mask_vec, MaskType (*mask_arr)[s_num_distinct_chars],
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

    BitMaskVec src = _mm256_set_epi32(val7, val6, val5, val4, val3, val2, val1, val0);
    BitMaskVec right_shift = _mm256_set_epi32(s_mc_7, s_mc_6, s_mc_5, s_mc_4, s_mc_3, s_mc_2, s_mc_1, s_mc_0);
    BitMaskVec mask_vec = _mm256_srlv_epi32(src, right_shift);
    bit_mask_vec = _mm256_or_si256(mask_vec, last_mask_shift_out);

    BitMaskVec reserved_mask = _mm256_set_epi32(s_rm_7, s_rm_6, s_rm_5, s_rm_4, s_rm_3, s_rm_2, s_rm_1, s_rm_0);
    BitMaskVec reserved_src = _mm256_and_si256(src, reserved_mask);
    BitMaskVec left_shift = _mm256_set_epi32(s_sc_19, s_sc_1A, s_sc_1B, s_sc_1C, s_sc_1D, s_sc_1E, s_sc_1F, s_sc_20);
    last_mask_shift_out = _mm256_sllv_epi32(reserved_src, left_shift);
}

// 每个 cc 文件单独修改
static inline void compute_dist_vec_avx2_s(BitMaskVec &bit_mask_vec, SimdType &distm_chosen, const SimdType &distm,
                                           const SimdType &_1_distm)
{
    distm_chosen = _mm256_blendv_ps(distm, _1_distm, _mm256_castsi256_ps(bit_mask_vec));

    bit_mask_vec = _mm256_slli_epi32(bit_mask_vec, 1);
}

static inline void compute_mxy_avx2_s(SimdType &m_t, SimdType &m_t_y, SimdType &x_t, SimdType &y_t, const SimdType &m_t_2,
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
static inline void vector_shift_avx2_s(SimdType &x, MainType shift_in, MainType &shift_out)
{
    SimdType reversed_x = _mm256_permutevar8x32_ps(x, _mm256_set_epi32(6, 5, 4, 3, 2, 1, 0, 7));
    shift_out = _mm256_cvtss_f32(reversed_x);
    x = _mm256_blend_ps(reversed_x, _mm256_set1_ps(shift_in), 0b00000001);
}

// 每个 cc 文件单独修改
static inline void vector_shift_last_avx2_s(SimdType &x, MainType shift_in)
{
    SimdType reversed_x = _mm256_permutevar8x32_ps(x, _mm256_set_epi32(6, 5, 4, 3, 2, 1, 0, 7));
    x = _mm256_blend_ps(reversed_x, _mm256_set1_ps(shift_in), 0b00000001);
}

float compute_fp_avx2_s(const TestCase &tc)
{
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
    precompute_masks_avx2_s(tc, cols, num_mask_vecs, mask_arr);

    uint8_t rs_arr[s_avx_length];
    BitMaskVec bit_mask_vec;
    BitMaskVec last_mask_shift_out;
    uint32_t shift_idx;
    uint32_t num_mask_bits_to_process;

    /* precompute initialization for probabilities and shift vector*/
    initialize_vectors_avx2_s(rows, cols, shift_out_m, shift_out_x, shift_out_y, ctx, tc, pa, distm1d_arr);

    for (uint32_t i{0}; i < stripe_cnt - 1; ++i) {
        {
            // 预取 init_masks_for_row_avx2_s 需要的内存

            _mm_prefetch(tc.rs + i * s_avx_length, _MM_HINT_T0);
        }

        stripe_initialization_avx2_s(i, p_gapm, p_mm, p_mx, p_xx, p_my, p_yy, distm, _1_distm, distm1d_arr, pa, m_t_1, m_t_2, x_t_1, x_t_2,
                                     y_t_1, y_t_2, m_t_1_y, shift_out_x, shift_out_m, tc);

        init_masks_for_row_avx2_s(tc, rs_arr, last_mask_shift_out, i * s_avx_length + 1, s_avx_length);

        for (uint32_t begin_d{1}; begin_d < cols + s_avx_length; begin_d += s_main_type_size) {
            num_mask_bits_to_process = std::min(s_main_type_size, cols + s_avx_length - begin_d);

            update_masks_for_cols_avx2_s((begin_d - 1) / s_main_type_size, bit_mask_vec, mask_arr, rs_arr, last_mask_shift_out);

            for (uint32_t mbi{0}; mbi < num_mask_bits_to_process; ++mbi) {
                shift_idx = begin_d + mbi + s_avx_length;

                compute_dist_vec_avx2_s(bit_mask_vec, distm_chosen, distm, _1_distm);

                compute_mxy_avx2_s(m_t, m_t_y, x_t, y_t, m_t_2, x_t_2, y_t_2, m_t_1, x_t_1, m_t_1_y, y_t_1, p_mm, p_gapm, p_mx, p_xx, p_my,
                                   p_yy, distm_chosen);

                vector_shift_avx2_s(m_t, shift_out_m[shift_idx], shift_out_m[begin_d + mbi]);
                vector_shift_avx2_s(x_t, shift_out_x[shift_idx], shift_out_x[begin_d + mbi]);
                vector_shift_avx2_s(y_t_1, shift_out_y[shift_idx], shift_out_y[begin_d + mbi]);

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

        stripe_initialization_avx2_s(i, p_gapm, p_mm, p_mx, p_xx, p_my, p_yy, distm, _1_distm, distm1d_arr, pa, m_t_1, m_t_2, x_t_1, x_t_2,
                                     y_t_1, y_t_2, m_t_1_y, shift_out_x, shift_out_m, tc);

        if (remaining_rows == 0) {
            remaining_rows = s_avx_length;
        }
        init_masks_for_row_avx2_s(tc, rs_arr, last_mask_shift_out, i * s_avx_length + 1, remaining_rows);

        SimdType sum_m = vec_set1_val_f(zero);
        SimdType sum_x = vec_set1_val_f(zero);

        for (uint32_t begin_d{1}; begin_d < cols + remaining_rows - 1; begin_d += s_main_type_size) {
            num_mask_bits_to_process = std::min(s_main_type_size, cols + remaining_rows - 1 - begin_d);

            update_masks_for_cols_avx2_s((begin_d - 1) / s_main_type_size, bit_mask_vec, mask_arr, rs_arr, last_mask_shift_out);

            for (uint32_t mbi{0}; mbi < num_mask_bits_to_process; ++mbi) {
                shift_idx = begin_d + mbi + s_avx_length;

                compute_dist_vec_avx2_s(bit_mask_vec, distm_chosen, distm, _1_distm);

                compute_mxy_avx2_s(m_t, m_t_y, x_t, y_t, m_t_2, x_t_2, y_t_2, m_t_1, x_t_1, m_t_1_y, y_t_1, p_mm, p_gapm, p_mx, p_xx, p_my,
                                   p_yy, distm_chosen);

                sum_m = vec_add(sum_m, m_t);
                vector_shift_last_avx2_s(m_t, shift_out_m[shift_idx]);

                sum_x = vec_add(sum_x, x_t);
                vector_shift_last_avx2_s(x_t, shift_out_x[shift_idx]);

                vector_shift_last_avx2_s(y_t_1, shift_out_y[shift_idx]);

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
        result_avx2 = sum_mx.f[remaining_rows - 1];
    }

    return result_avx2;
}
