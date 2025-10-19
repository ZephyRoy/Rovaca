#include <assert.h>
#include <immintrin.h>
#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "avx2-functions.h"
#include "hc_marco.h"
#include "smithwaterman_common.h"

// the maximum DNA sequence length
#define MAX_SEQ_LEN       (2048)
#define MATRIX_MIN_CUTOFF (-100000000)
#define LOW_INIT_VALUE    (INT32_MIN / 2)

static void smith_waterman_back_track(p_lib_sw_avx p);
static inline uint32_t sw_avx_get_new_cigar(uint8_t cigar_op, uint32_t cigar_len);
static inline uint32_t sw_avx_make_element(uint32_t state, int length);
static void sw_avx_get_cigar(p_lib_sw_avx p);
// static int sw_avx_relloc(p_lib_sw_avx ret);

#define MAIN_CODE(bt_vec)                                                   \
    {                                                                       \
        VEC_INT_TYPE e10 = VEC_LOADU(&E[inde]);                             \
        VEC_INT_TYPE ext_score_h = VEC_ADD(e10, w_extend_vec);              \
        VEC_INT_TYPE h10 = VEC_LOADU(&H[hLeftInd]);                         \
        VEC_INT_TYPE open_score_h = VEC_ADD(h10, w_open_vec);               \
        VEC_INT_TYPE e11 = VEC_MAX(open_score_h, ext_score_h);              \
        VEC_INT_TYPE open_gt_ext_h = VEC_CMPGT(open_score_h, ext_score_h);  \
        VEC_INT_TYPE ext_vec = VEC_ANDNOT(open_gt_ext_h, ins_ext_vec);      \
        VEC_STOREU(&E[inde], e11);                                          \
        VEC_INT_TYPE f01 = VEC_LOADU(&F[indf]);                             \
        VEC_INT_TYPE f11;                                                   \
        VEC_INT_TYPE ext_score_v = VEC_ADD(f01, w_extend_vec);              \
        VEC_INT_TYPE h01 = VEC_LOADU(&H[hTopInd]);                          \
        VEC_INT_TYPE open_score_v = VEC_ADD(h01, w_open_vec);               \
        f11 = VEC_MAX(ext_score_v, open_score_v);                           \
        VEC_INT_TYPE open_gt_ext_v = VEC_CMPGT(open_score_v, ext_score_v);  \
        ext_vec = VEC_OR(ext_vec, VEC_ANDNOT(open_gt_ext_v, del_ext_vec));  \
        VEC_STOREU((VEC_INT_TYPE *)(&F[indf]), f11);                        \
        VEC_INT_TYPE h00 = VEC_LOADU((VEC_INT_TYPE *)(&H[hCurInd]));        \
        VEC_INT_TYPE s1 = VEC_LOADU((VEC_INT_TYPE *)(seq1_rev + inde));     \
        VEC_INT_TYPE s2 = VEC_LOADU((VEC_INT_TYPE *)(seq2 + seq2Ind));      \
        VEC_MASK_TYPE cmp11 = VEC_CMPEQ_MASK(s1, s2);                       \
        VEC_INT_TYPE sbt11 = VEC_BLEND(w_mismatch_vec, w_match_vec, cmp11); \
        VEC_INT_TYPE m11 = VEC_ADD(h00, sbt11);                             \
        VEC_INT_TYPE h11 = VEC_MAX(min_cutoff_vec, m11);                    \
        VEC_INT_TYPE e11_gt_h11 = VEC_CMPGT(e11, h11);                      \
        h11 = VEC_MAX(h11, e11);                                            \
        bt_vec = VEC_AND(ins_vec, e11_gt_h11);                              \
        cmp11 = VEC_CMPGT_MASK(f11, h11);                                   \
        h11 = VEC_MAX(h11, f11);                                            \
        bt_vec = VEC_BLEND(bt_vec, del_vec, cmp11);                         \
        bt_vec = VEC_OR(bt_vec, ext_vec);                                   \
        VEC_STOREU((VEC_INT_TYPE *)(&H[hCurInd]), h11);                     \
    }

static void smith_waterman_back_track(p_lib_sw_avx p)
{
    int32_t max_seq_len = p->max_seq_len;
    uint32_t *seq1_rev = p->cacl_cache.seq1_rev;
    uint32_t *seq2 = p->cacl_cache.seq2;
    int32_t nrow = p->len1;
    int32_t ncol = p->len2;
    int32_t match = p->para.match;
    int32_t mismatch = p->para.mismatch;
    int32_t open = p->para.gap_open;
    int32_t extend = p->para.gap_extend;

    VEC_INT_TYPE w_match_vec = VEC_SET1_VAL32(match);
    VEC_INT_TYPE w_mismatch_vec = VEC_SET1_VAL32(mismatch);
    VEC_INT_TYPE w_open_vec = VEC_SET1_VAL32(open);
    VEC_INT_TYPE w_extend_vec = VEC_SET1_VAL32(extend);

    int32_t i, j;
    int32_t low_init_value = LOW_INIT_VALUE;
    VEC_INT_TYPE low_init_value_vec = VEC_SET1_VAL32(low_init_value);
    VEC_INT_TYPE min_cutoff_vec = VEC_SET1_VAL32(MATRIX_MIN_CUTOFF);
    VEC_INT_TYPE ins_vec = VEC_SET1_VAL32(LIB_SW_AVX_OVERHANG_STRATEGY_INSERT);
    VEC_INT_TYPE del_vec = VEC_SET1_VAL32(LIB_SW_AVX_OVERHANG_STRATEGY_DELETE);
    VEC_INT_TYPE ins_ext_vec = VEC_SET1_VAL32(LIB_SW_AVX_OVERHANG_STRATEGY_INSERT_EXT);
    VEC_INT_TYPE del_ext_vec = VEC_SET1_VAL32(LIB_SW_AVX_OVERHANG_STRATEGY_DELETE_EXT);
    INIT_CONSTANTS;
    int32_t hwidth = max_seq_len + AVX_LENGTH;
    int32_t ewidth = max_seq_len + AVX_LENGTH;

    int32_t *E = p->mm_aloc.e;
    int32_t *F = E + 1 * ewidth;
    int32_t *H = E + 2 * ewidth;
    for (j = 0; j <= ncol; j += AVX_LENGTH) {
        VEC_STORE(F + j, low_init_value_vec);
    }
    for (i = 0; i <= nrow; i += AVX_LENGTH) {
        VEC_STORE(E + i, low_init_value_vec);
    }

    H[max_seq_len >> 1] = 0;

    for (i = 0; i < nrow; i++) {
        seq1_rev[max_seq_len - 1 - i] = p->seq1[i];
    }
    for (i = 0; i < ncol; i++) {
        seq2[i] = p->seq2[i];
    }

    int16_t *back_track = p->btrack;

    int32_t max_score = INT32_MIN;
    int32_t max_i = 0;
    int32_t max_j = 0;

    int32_t anti_diag;
    int32_t prev, cur;
    for (anti_diag = 1; anti_diag <= (nrow + ncol); anti_diag++) {
        int32_t ilo = assemble_min(anti_diag, nrow + 1);
        int32_t jhi = assemble_min(anti_diag, ncol + 1);
        int32_t ihi = anti_diag - jhi;
        int32_t jlo = anti_diag - ilo;

        prev = ((anti_diag - 1) & 1) * hwidth;
        cur = (anti_diag & 1) * hwidth;

        for (j = (jlo + 1); j < (jhi - AVX_LENGTH);) {
            int32_t back_track_ind = j - jlo - 1;
            VEC_INT_TYPE bt_vec_0, bt_vec_1;
            {
                i = anti_diag - j;
                int32_t diag = j - i;
                int32_t diagInd = max_seq_len + diag;
                int32_t inde = max_seq_len - i;
                int32_t indf = j;
                int32_t hLeftInd = prev + ((diagInd - 1) >> 1);
                int32_t hTopInd = hLeftInd + 1;
                int32_t hCurInd = cur + (diagInd >> 1);
                int32_t seq2Ind = j - 1;
                MAIN_CODE(bt_vec_0)
                j = j + AVX_LENGTH;
            }
            {
                i = anti_diag - j;
                int32_t diag = j - i;
                int32_t diagInd = max_seq_len + diag;
                int32_t inde = max_seq_len - i;
                int32_t indf = j;
                int32_t hLeftInd = prev + ((diagInd - 1) >> 1);
                int32_t hTopInd = hLeftInd + 1;
                int32_t hCurInd = cur + (diagInd >> 1);
                int32_t seq2Ind = j - 1;
                MAIN_CODE(bt_vec_1)
                j = j + AVX_LENGTH;
            }
            VEC_INT_TYPE bt_vec_2 = VEC_PERMUTE2x128_EVEN(bt_vec_0, bt_vec_1);
            VEC_INT_TYPE bt_vec_3 = VEC_PERMUTE2x128_ODD(bt_vec_0, bt_vec_1);
            VEC_INT_TYPE bt_vec = VEC_PACKS_32(bt_vec_2, bt_vec_3);
            VEC_STREAM(back_track + anti_diag * max_seq_len + back_track_ind, bt_vec);
        }
        if (j < jhi) {
            i = anti_diag - j;
            int32_t diag = j - i;
            int32_t diagInd = max_seq_len + diag;
            int32_t inde = max_seq_len - i;
            int32_t indf = j;
            int32_t hLeftInd = prev + ((diagInd - 1) >> 1);
            int32_t hTopInd = hLeftInd + 1;
            int32_t hCurInd = cur + (diagInd >> 1);
            int32_t seq2Ind = j - 1;
            VEC_INT_TYPE bt_vec_0;
            VEC_INT_TYPE bt_vec_1 = VEC_SET_ZERO();
            MAIN_CODE(bt_vec_0)
            VEC_INT_TYPE bt_vec_2 = VEC_PERMUTE2x128_EVEN(bt_vec_0, bt_vec_1);
            VEC_INT_TYPE bt_vec_3 = VEC_PERMUTE2x128_ODD(bt_vec_0, bt_vec_1);
            VEC_INT_TYPE bt_vec = VEC_PACKS_32(bt_vec_2, bt_vec_3);
            int32_t back_track_ind = j - jlo - 1;
            VEC_STREAM(back_track + anti_diag * max_seq_len + back_track_ind, bt_vec);
        }
        H[cur + ((max_seq_len + 2 * jhi - anti_diag) >> 1)] = 0;
        H[cur + ((max_seq_len + 2 * jlo - anti_diag) >> 1)] = 0;
        F[jhi] = low_init_value;
        E[max_seq_len - ilo] = low_init_value;

        if (ilo == (nrow + 1)) {
            int32_t score = H[cur + ((max_seq_len + jlo + 1 - (ilo - 1)) >> 1)];
            if ((max_score < score) || ((max_score == score) && (abs(ilo - jlo - 2) < abs(max_i - max_j)))) {
                max_score = score;
                max_i = ilo - 1;
                max_j = jlo + 1;
            }
        }
        if (jhi == (ncol + 1)) {
            int32_t score = H[cur + ((max_seq_len + jhi - 1 - (ihi + 1)) >> 1)];
            if ((max_score < score) || ((max_score == score) && ((max_j == ncol) || (abs(ihi - jhi + 2) <= abs(max_i - max_j))))) {
                max_score = score;
                max_i = ihi + 1;
                max_j = jhi - 1;
            }
        }
    }
    p->score = max_score;
    p->score = max_score;
    p->max_i = (int16_t)max_i;
    p->max_j = (int16_t)max_j;
    return;
}

/**
 * @brief 获取新的单一cigar
 *
 * @param[in]   cigar_op          cigar_op
 * @param[out]  cigar_len         cigar_len
 *
 * @retval uint32_t               cigar
 **/
static inline uint32_t sw_avx_get_new_cigar(uint8_t cigar_op, uint32_t cigar_len)
{
    uint32_t ret;

    ret = cigar_len << BAM_CIGAR_SHIFT;
    ret |= cigar_op;
    return ret;
}

static inline uint32_t sw_avx_make_element(uint32_t state, int length)
{
    uint8_t op = WORKER_SIGAR_MAX;
    switch (state) {
        case LIB_SW_AVX_OVERHANG_STRATEGY_MATCH: op = WORKER_SIGAR_STATUS_MATCH; break;
        case LIB_SW_AVX_OVERHANG_STRATEGY_INSERT: op = WORKER_SIGAR_STATUS_INSERTION; break;
        case LIB_SW_AVX_OVERHANG_STRATEGY_DELETE: op = WORKER_SIGAR_STATUS_DELETION; break;
        case LIB_SW_AVX_OVERHANG_STRATEGY_SOFT_CLIP: op = WORKER_SIGAR_SOFT_CLIP; break;
        default: return 0;
    }
    return sw_avx_get_new_cigar(op, length);
}

static void sw_avx_get_cigar(p_lib_sw_avx p)
{
    int16_t *btrack = p->btrack;
    int16_t max_i = p->max_i;
    int16_t max_j = p->max_j;
    int16_t nrow = p->len1;
    int16_t ncol = p->len2;
    int32_t max_seq_len = p->max_seq_len;
    int32_t cigar_buf_length = p->cigar_max_len;
    int16_t i, j;
    int16_t *cigar_array = p->mm_aloc.cigar_buf;

    int32_t cigar_id = 0;

    i = max_i;
    j = max_j;

    if (j < ncol) {
        cigar_array[cigar_id * 2] = LIB_SW_AVX_OVERHANG_STRATEGY_SOFT_CLIP;
        // ncol - j is automatically converted to an int because the result of int16_t - int16_t may overflow int16_t
        // That never happens in this case as both are guaranteed to be positive and j < ncol.
        cigar_array[cigar_id * 2 + 1] = (int16_t)(ncol - j);
        cigar_id++;
    }
    int state = 0;
    while ((i > 0) && (j > 0)) {
        int32_t anti_diag = i + j;
        int32_t btrInd;
        if (anti_diag <= nrow) {
            btrInd = anti_diag * max_seq_len + j - 1;
        }
        else {
            int32_t jLo = anti_diag - nrow - 1;
            btrInd = anti_diag * max_seq_len + j - jLo - 1;
        }
        int32_t btr = btrack[btrInd];
        if (state == LIB_SW_AVX_OVERHANG_STRATEGY_INSERT_EXT) {
            j--;
            cigar_array[cigar_id * 2 - 1]++;
            state = btr & LIB_SW_AVX_OVERHANG_STRATEGY_INSERT_EXT;
        }
        else if (state == LIB_SW_AVX_OVERHANG_STRATEGY_DELETE_EXT) {
            i--;
            cigar_array[cigar_id * 2 - 1]++;
            state = btr & LIB_SW_AVX_OVERHANG_STRATEGY_DELETE_EXT;
        }
        else {
            switch (btr & 3) {
                case LIB_SW_AVX_OVERHANG_STRATEGY_MATCH:
                    i--;
                    j--;
                    cigar_array[cigar_id * 2] = LIB_SW_AVX_OVERHANG_STRATEGY_MATCH;
                    cigar_array[cigar_id * 2 + 1] = 1;
                    state = 0;
                    cigar_id++;
                    break;
                case LIB_SW_AVX_OVERHANG_STRATEGY_INSERT:
                    j--;
                    cigar_array[cigar_id * 2] = LIB_SW_AVX_OVERHANG_STRATEGY_INSERT;
                    cigar_array[cigar_id * 2 + 1] = 1;
                    state = btr & LIB_SW_AVX_OVERHANG_STRATEGY_INSERT_EXT;
                    cigar_id++;
                    break;
                case LIB_SW_AVX_OVERHANG_STRATEGY_DELETE:
                    i--;
                    cigar_array[cigar_id * 2] = LIB_SW_AVX_OVERHANG_STRATEGY_DELETE;
                    cigar_array[cigar_id * 2 + 1] = 1;
                    state = btr & LIB_SW_AVX_OVERHANG_STRATEGY_DELETE_EXT;
                    cigar_id++;
                    break;
            }
        }
    }
    if (j > 0) {
        cigar_array[cigar_id * 2] = LIB_SW_AVX_OVERHANG_STRATEGY_SOFT_CLIP;
        cigar_array[cigar_id * 2 + 1] = j;
        cigar_id++;
    }
    p->alignment_offset = i;
    int16_t new_id = 0;
    int16_t prev = cigar_array[new_id * 2];
    // printf("cigar_id = %d\n", cigar_id);
    for (i = 1; i < cigar_id; i++) {
        int16_t cur = cigar_array[i * 2];
        if (cur == prev) {
            // cigar_array elements will never overflow int16_t because they're bounded by the 16 bit input length.
            cigar_array[new_id * 2 + 1] = (int16_t)(cigar_array[new_id * 2 + 1] + cigar_array[i * 2 + 1]);
        }
        else {
            new_id++;
            cigar_array[new_id * 2] = cur;
            cigar_array[new_id * 2 + 1] = cigar_array[i * 2 + 1];
            prev = cur;
        }
    }

    int cur_size = 0;
    for (i = new_id; i >= 0; i--) {
        if (__glibc_unlikely(cur_size < 0 || cigar_array[2 * i + 1] <= 0 || cur_size + 1 >= cigar_buf_length)) {
            continue;
        }
        p->cigar[cur_size] = sw_avx_make_element(cigar_array[2 * i], cigar_array[2 * i + 1]);
        cur_size += 1;
    }
    // CIGARs are specified to have their length fit in uint32_t
    p->cigar_count = cur_size;
}

void sw_avx2_run(p_lib_sw_avx p)
{
    int len = assemble_max(p->len1, p->len2);

    while (len >= p->max_seq_len) {
        p->alignment_offset = 1;
        return;
    }

    smith_waterman_back_track(p);
    sw_avx_get_cigar(p);
}

#ifdef MAIN_TEST

/**
 * Make sure that the SW didn't fail in some terrible way, and throw exception if it did
 */
static int hc_assemle_cigar_cacl_is_sw_failure(uint32_t *cigar, int32_t cigar_len, int32_t alignment_offset)
{
    uint32_t this_cigar;
    register int i;

    // check that the alignment starts at the first base, which it should given the padding
    if (alignment_offset > 0) {
        return 1;
    }

    // check that we aren't getting any S operators (which would be very bad downstream)
    for (i = 0; i < cigar_len; i++) {
        this_cigar = bam_cigar_op(cigar[i]);

        // soft clips at the end of the alignment are really insertions
        if (this_cigar == WORKER_SIGAR_SOFT_CLIP) {
            return 1;
        }
    }

    return 0;
}

void test_main_avx(void)
{
    FILE *fp = fopen("/data/rovaca-dev/caoce/code/test/smithwaterman/seq.log", "r");
    const char *cigar_string = "MIDNSHP=XB";
    unsigned int i;

    p_lib_sw_avx run = sw_avx_init();

    while (!feof(fp)) {
        fgets((char *)run->seq1, 1024, fp);
        fgets((char *)run->seq2, 1024, fp);
        run->len1 = strlen((char *)run->seq1);
        run->len2 = strlen((char *)run->seq2);

        for (i = 0; i < 10; i++) {
            sw_avx2_run(run);
        }
        continue;
        if (hc_assemle_cigar_cacl_is_sw_failure(run->cigar, run->cigar_count, run->alignment_offset)) {
            puts("failed");
            continue;
        }

        printf("%d ", run->alignment_offset);
        for (i = 0; i < run->cigar_count; i++) {
            printf("%u%c", bam_cigar_oplen(run->cigar[i]), cigar_string[bam_cigar_op(run->cigar[i])]);
        }
        puts("");
    }
    fclose(fp);
}

int main(void)
{
    test_main_avx();
    return 0;
}
#endif