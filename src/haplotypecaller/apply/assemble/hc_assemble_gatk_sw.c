/**
 * Pairwise discrete smith-waterman alignment implemented in pure java
 *
 * ************************************************************************
 * ****                    IMPORTANT NOTE:                             ****
 * ****  This class assumes that all bytes come from UPPERCASED chars! ****
 * ************************************************************************
 */

#include "hc_assemble_gatk_sw.h"

#include <limits.h>
#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "hc_func.h"
#include "hc_marco.h"

#define GATK_SW_MATRIX_MIN_CUTOFF       ((int)-1.0e8)  // never let matrix elements drop below this cutoff
#define GATK_SW_MIN(A, B)               (((A) < (B)) ? (A) : (B))
#define GATK_SW_MAX(A, B)               (((A) > (B)) ? (A) : (B))
#define GATK_SW_GET_SW_X(align, x)      (align->cacl_cache.sw + align->cacl_cache.y_len * x)
#define GATK_SW_GET_SW_POS(align, x, y) ((align->cacl_cache.sw + align->cacl_cache.y_len * x) + y)
#define GATK_SW_GET_BT_X(align, x)      (align->cacl_cache.btrack + align->cacl_cache.y_len * x)
#define GATK_SW_GET_BT_POS(align, x, y) ((align->cacl_cache.btrack + align->cacl_cache.y_len * x) + y)

#ifdef __SW_JAVA_TEST__
#define BAM_CIGAR_SHIFT (4)

// sam内cigar状态
enum WORKER_CIGAR_PARSE_STATUS {
    WORKER_SIGAR_STATUS_MATCH = 0,      // M 匹配上了
    WORKER_SIGAR_STATUS_INSERTION = 1,  // I 插入点，loop这次循环
    WORKER_SIGAR_STATUS_DELETION = 2,   // D 缺少点，跳到下一个slot写数据
    WORKER_SIGAR_STATUS_NSTATUS = 3,    // N 屏蔽点，不增加 = D
    WORKER_SIGAR_SOFT_CLIP = 4,         // S 屏蔽点，不增加
    WORKER_SIGAR_HARD_CLIP = 5,         // H 屏蔽点，不增加
    WORKER_SIGAR_PADDING = 6,           // P 屏蔽点，不增加
    WORKER_SIGAR_S_MATCH = 7,           // = = M
    WORKER_SIGAR_S_MISMATCH = 8,        // X = M
    WORKER_SIGAR_MAX = 9
};
#endif

/**
 * The state of a trace step through the matrix
 */
enum GATK_SW_STATE { GATK_SW_STATE_MATCH, GATK_SW_STATE_INSERTION, GATK_SW_STATE_DELETION, GATK_SW_STATE_CLIP };

#ifdef __SW_JAVA_TEST__
static inline uint32_t hc_assemble_utils_get_new_cigar(uint8_t cigar_op, uint32_t cigar_len);
#endif
static int gatk_sw_last_index_of(uint8_t* reference, uint8_t* query, int ref_len, int que_len);
static void gatk_sw_fill_cache(int* cache, uint32_t cache_len, int value);
static void gatk_sw_calculate_matrix(p_gatk_sw_storage align);
static void gatk_sw_calculate_cigar(p_gatk_sw_storage align);
static inline uint32_t gatk_sw_make_element(uint32_t state, int length);

/**
 * Find the last occurrence of the query sequence in the reference sequence
 *
 * Returns the index of the last occurrence or -1 if the query sequence is not found
 *
 * @param reference the reference sequence
 * @param query the query sequence
 */
static int gatk_sw_last_index_of(uint8_t* reference, uint8_t* query, int ref_len, int que_len)
{
    int queryLength = que_len;
    int r = ref_len - queryLength;

    // start search from the last possible matching position and search to the left
    for (; r >= 0; r--) {
        int q = 0;
        while (q < queryLength && reference[r + q] == query[q]) {
            q++;
        }
        if (q == queryLength) {
            return r;
        }
    }
    return -1;
}

static void gatk_sw_fill_cache(int* cache, uint32_t cache_len, int value)
{
    register uint32_t i = 0;

    for (; i < cache_len; i++) {
        cache[i] = value;
    }
}

/**
 * @brief Aligns the alternate sequence to the reference sequence
 *
 * @param align     Aligns
 *
 * @return int      1:succeed, 0:failed
 */
int hc_assemble_gatk_sw_align(p_gatk_sw_storage align)
{
    // avoid running full Smith-Waterman if there is an exact match of alternate in reference
    int matchIndex = -1;

    if (align->input.ref == NULL || align->input.ref_len == 0 || align->input.alt == NULL || align->input.alt_len == 0) {
        return 0;
    }

    if (align->input.overhang_strategy == GATK_SW_OVERHANG_STRATEGY_SOFTCLIP ||
        align->input.overhang_strategy == GATK_SW_OVERHANG_STRATEGY_IGNORE) {
        // Use a substring search to find an exact match of the alternate in the reference
        // NOTE: This approach only works for SOFTCLIP and IGNORE overhang strategies
        matchIndex = gatk_sw_last_index_of(align->input.ref, align->input.alt, align->input.ref_len, align->input.alt_len);
    }

    if (matchIndex != -1) {
        // generate the alignment result when the substring search was successful
        align->output.alignment_offset = matchIndex;
        align->output.cigar[0] = hc_assemble_utils_get_new_cigar(WORKER_SIGAR_STATUS_MATCH, align->input.alt_len);
        align->output.cigar_len = 1;
    }
    else {
        int loop_size;
        int lowInitValue = INT_MIN / 2;

        if (align->input.ref_len == 0 || align->input.alt_len == 0) {
            return 0;
        }

        align->cacl_cache.current_data = align->cacl_cache.data;
        align->cacl_cache.x_len = align->input.ref_len + 1;
        align->cacl_cache.y_len = align->input.alt_len + 1;
        loop_size = align->cacl_cache.x_len * align->cacl_cache.y_len;
        align->cacl_cache.sw = align->cacl_cache.current_data;
        align->cacl_cache.current_data += loop_size;
        align->cacl_cache.btrack = align->cacl_cache.current_data;
        align->cacl_cache.current_data += loop_size;
        align->cacl_cache.best_gap_v = align->cacl_cache.current_data;                              // align->cacl_cache.y_len + 1
        align->cacl_cache.gap_size_v = align->cacl_cache.best_gap_v + align->cacl_cache.y_len + 1;  // align->cacl_cache.y_len + 1
        align->cacl_cache.best_gap_h = align->cacl_cache.gap_size_v + align->cacl_cache.y_len + 1;  // align->cacl_cache.x_len + 1
        align->cacl_cache.gap_size_h = align->cacl_cache.best_gap_h + align->cacl_cache.x_len + 1;  // align->cacl_cache.x_len + 1

        gatk_sw_fill_cache(align->cacl_cache.sw, loop_size, 0);
        gatk_sw_fill_cache(align->cacl_cache.btrack, loop_size, 0);

        gatk_sw_fill_cache(align->cacl_cache.best_gap_v, align->cacl_cache.y_len + 1, lowInitValue);
        gatk_sw_fill_cache(align->cacl_cache.gap_size_v, align->cacl_cache.y_len + 1, 0);
        gatk_sw_fill_cache(align->cacl_cache.best_gap_h, align->cacl_cache.x_len + 1, lowInitValue);
        gatk_sw_fill_cache(align->cacl_cache.gap_size_h, align->cacl_cache.x_len + 1, 0);

        // run full Smith-Waterman
        gatk_sw_calculate_matrix(align);
        gatk_sw_calculate_cigar(align);  // length of the segment (continuous matches, insertions or deletions)
    }

    return 1;
}

/**
 * Calculates the SW matrices for the given sequences
 */
static void gatk_sw_calculate_matrix(p_gatk_sw_storage align)
{
    int32_t* curRow;
    uint32_t i, j;

    // we need to initialize the SW matrix with gap penalties if we want to keep track of indels at the edges of alignments
    if (align->input.overhang_strategy == GATK_SW_OVERHANG_STRATEGY_INDEL ||
        align->input.overhang_strategy == GATK_SW_OVERHANG_STRATEGY_LEADING_INDEL) {
        // initialize the first row
        int* topRow = align->cacl_cache.sw;
        int currentValue = align->paramates.gapOpenPenalty;

        topRow[1] = align->paramates.gapOpenPenalty;
        for (i = 2; i < align->cacl_cache.y_len; i++) {
            currentValue += align->paramates.gapExtendPenalty;
            topRow[i] = currentValue;
        }
        // initialize the first column
        *GATK_SW_GET_SW_POS(align, 1, 0) = align->paramates.gapOpenPenalty;
        currentValue = align->paramates.gapOpenPenalty;
        for (i = 2; i < align->cacl_cache.x_len; i++) {
            currentValue += align->paramates.gapExtendPenalty;
            *GATK_SW_GET_SW_POS(align, i, 0) = currentValue;
        }
    }
    // build smith-waterman matrix and keep backtrack info:
    curRow = align->cacl_cache.sw;

    // array length checks are expensive in tight loops so extract the length out
    for (i = 1; i < align->cacl_cache.x_len; i++) {
        uint8_t a_base = align->input.ref[i - 1];  // letter in a at the current pos
        int* lastRow = curRow;
        int* curBackTrackRow;

        curRow = GATK_SW_GET_SW_X(align, i);
        curBackTrackRow = GATK_SW_GET_BT_X(align, i);

        // array length checks are expensive in tight loops so extract the length out
        for (j = 1; j < align->cacl_cache.y_len; j++) {
            uint8_t b_base = align->input.alt[j - 1];  // letter in b at the current pos
            // in other words, step_diag = sw[i-1][j-1] + wd(a_base,b_base);
            int step_diag = lastRow[j - 1] + (a_base == b_base ? align->paramates.matchValue : align->paramates.mismatchPenalty);

            // optimized "traversal" of all the matrix cells above the current one (i.e. traversing
            // all 'step down' events that would end in the current cell. The optimized code
            // does exactly the same thing as the commented out loop below. IMPORTANT:
            // the optimization works ONLY for linear w(k)=wopen+(k-1)*wextend!!!!

            // if a gap (length 1) was just opened above, this is the cost of arriving to the current cell:
            int prev_gap = lastRow[j] + align->paramates.gapOpenPenalty;

            align->cacl_cache.best_gap_v[j] +=
                align->paramates.gapExtendPenalty;  // for the gaps that were already opened earlier, extending them by 1 costs
                                                    // align->paramates.gapExtendPenalty
            if (prev_gap > align->cacl_cache.best_gap_v[j]) {
                // opening a gap just before the current cell results in better score than extending by one
                // the best previously opened gap. This will hold for ALL cells below: since any gap
                // once opened always costs align->paramates.gapExtendPenalty to extend by another base, we will always get a better score
                // by arriving to any cell below from the gap we just opened (prev_gap) rather than from the previous best gap
                align->cacl_cache.best_gap_v[j] = prev_gap;
                align->cacl_cache.gap_size_v[j] = 1;  // remember that the best step-down gap from above has length 1 (we just opened it)
            }
            else {
                // previous best gap is still the best, even after extension by another base, so we just record that extension:
                align->cacl_cache.gap_size_v[j]++;
            }

            int step_down = align->cacl_cache.best_gap_v[j];
            int kd = align->cacl_cache.gap_size_v[j];

            // optimized "traversal" of all the matrix cells to the left of the current one (i.e. traversing
            // all 'step right' events that would end in the current cell. The optimized code
            // does exactly the same thing as the commented out loop below. IMPORTANT:
            // the optimization works ONLY for linear w(k)=wopen+(k-1)*wextend!!!!

            prev_gap = curRow[j - 1] +
                       align->paramates.gapOpenPenalty;  // what would it cost us to open length 1 gap just to the left from current cell
            align->cacl_cache.best_gap_h[i] +=
                align->paramates.gapExtendPenalty;  // previous best gap would cost us that much if extended by another base
            if (prev_gap > align->cacl_cache.best_gap_h[i]) {
                // newly opened gap is better (score-wise) than any previous gap with the same row index i; since
                // gap penalty is linear with k, this new gap location is going to remain better than any previous ones
                align->cacl_cache.best_gap_h[i] = prev_gap;
                align->cacl_cache.gap_size_h[i] = 1;
            }
            else {
                align->cacl_cache.gap_size_h[i]++;
            }

            int step_right = align->cacl_cache.best_gap_h[i];
            int ki = align->cacl_cache.gap_size_h[i];

            // priority here will be step diagonal, step right, step down
            if (step_diag < step_down || step_diag < step_right) {
                if (step_right >= step_down) {  // moving right is the highest

                    curRow[j] = GATK_SW_MAX(GATK_SW_MATRIX_MIN_CUTOFF, step_right);
                    curBackTrackRow[j] = -ki;  // negative = horizontal
                }
                else {
                    curRow[j] = GATK_SW_MAX(GATK_SW_MATRIX_MIN_CUTOFF, step_down);
                    curBackTrackRow[j] = kd;  // positive=vertical
                }
            }
            else {
                curRow[j] = GATK_SW_MAX(GATK_SW_MATRIX_MIN_CUTOFF, step_diag);
                curBackTrackRow[j] = 0;
            }
        }
    }
}

/**
 * Calculates the CIGAR for the alignment from the back track matrix
 *
 * @param sw                   the Smith-Waterman matrix to use
 * @param btrack               the back track matrix to use
 * @param overhangStrategy    the strategy to use for dealing with overhangs
 * @return non-null SWPairwiseAlignmentResult object
 */
static void gatk_sw_calculate_cigar(p_gatk_sw_storage align)
{
    // p holds the position we start backtracking from; we will be assembling a cigar in the backwards order
    int p1 = 0, p2 = 0, i, j;

    int maxscore = INT_MIN;  // sw scores are allowed to be negative
    int segment_length = 0;  // length of the segment (continuous matches, insertions or deletions)

    // if we want to consider overhangs as legitimate operators, then just start from the corner of the matrix
    if (align->input.overhang_strategy == GATK_SW_OVERHANG_STRATEGY_INDEL) {
        p1 = align->input.ref_len;
        p2 = align->input.alt_len;
    }
    else {
        // look for the largest score on the rightmost column. we use >= combined with the traversal direction
        // to ensure that if two scores are equal, the one closer to diagonal gets picked
        // Note: this is not technically smith-waterman, as by only looking for max values on the right we are
        // excluding high scoring local alignments
        p2 = align->input.alt_len;

        for (i = 1; i < (int)align->cacl_cache.x_len; i++) {
            if (*GATK_SW_GET_SW_POS(align, i, align->input.alt_len) >= maxscore) {
                p1 = i;
                maxscore = *GATK_SW_GET_SW_POS(align, i, align->input.alt_len);
            }
        }
        // now look for a larger score on the bottom-most row
        if (align->input.overhang_strategy != GATK_SW_OVERHANG_STRATEGY_LEADING_INDEL) {
            int* bottomRow = GATK_SW_GET_SW_X(align, align->input.ref_len);
            for (j = 1; j < (int)align->cacl_cache.y_len; j++) {
                int curScore = bottomRow[j];
                // data_offset is the offset of [n][j]
                if (curScore > maxscore || (curScore == maxscore && abs((int)align->input.ref_len - j) < abs(p1 - p2))) {
                    p1 = align->input.ref_len;
                    p2 = j;
                    maxscore = curScore;
                    segment_length = align->input.alt_len - j;  // end of sequence 2 is overhanging; we will just record it as 'M' segment
                }
            }
        }
    }
    uint32_t* lce = align->output.cigar;
    align->output.cigar_len = 0;
    if (segment_length > 0 && align->input.overhang_strategy == GATK_SW_OVERHANG_STRATEGY_SOFTCLIP) {
        lce[align->output.cigar_len] = gatk_sw_make_element(GATK_SW_STATE_CLIP, segment_length);
        align->output.cigar_len += 1;
        segment_length = 0;
    }

    // we will be placing all insertions and deletions into sequence b, so the states are named w/regard
    // to that sequence

    uint32_t state = GATK_SW_STATE_MATCH;
    do {
        int btr = *GATK_SW_GET_BT_POS(align, p1, p2);
        uint32_t new_state;
        int step_length = 1;
        if (btr > 0) {
            new_state = GATK_SW_STATE_DELETION;
            step_length = btr;
        }
        else if (btr < 0) {
            new_state = GATK_SW_STATE_INSERTION;
            step_length = (-btr);
        }
        else
            new_state = GATK_SW_STATE_MATCH;  // and step_length =1, already set above

        // move to next best location in the sw matrix:
        switch (new_state) {
            case GATK_SW_STATE_MATCH:
                p1--;
                p2--;
                break;                                               // move back along the diag in the sw matrix
            case GATK_SW_STATE_INSERTION: p2 -= step_length; break;  // move left
            case GATK_SW_STATE_DELETION: p1 -= step_length; break;   // move up
        }

        // now let's see if the state actually changed:
        if (new_state == state)
            segment_length += step_length;
        else {
            // state changed, lets emit previous segment, whatever it was (Insertion Deletion, or (Mis)Match).
            if (segment_length > 0) {
                lce[align->output.cigar_len] = gatk_sw_make_element(state, segment_length);
                align->output.cigar_len += 1;
            }
            segment_length = step_length;
            state = new_state;
        }
        // next condition is equivalent to  while ( sw[p1][p2] != 0 ) (with modified p1 and/or p2:
    } while (p1 > 0 && p2 > 0);

    // post-process the last segment we are still keeping;
    // NOTE: if reads "overhangs" the ref on the left (i.e. if p2>0) we are counting
    // those extra bases sticking out of the ref into the first cigar element if DO_SOFTCLIP is false;
    // otherwise they will be soft-clipped. For instance,
    // if read length is 5 and alignment starts at offset -2 (i.e. read starts before the ref, and only
    // last 3 bases of the read overlap with/align to the ref), the cigar will be still 5M if
    // DO_SOFTCLIP is false or 2S3M if DO_SOFTCLIP is true.
    // The consumers need to check for the alignment offset and deal with it properly.
    int alignment_offset;
    if (align->input.overhang_strategy == GATK_SW_OVERHANG_STRATEGY_SOFTCLIP) {
        lce[align->output.cigar_len] = gatk_sw_make_element(state, segment_length);
        align->output.cigar_len += 1;
        if (p2 > 0) {
            lce[align->output.cigar_len] = gatk_sw_make_element(GATK_SW_STATE_CLIP, p2);
            align->output.cigar_len += 1;
        }
        alignment_offset = p1;
    }
    else if (align->input.overhang_strategy == GATK_SW_OVERHANG_STRATEGY_IGNORE) {
        lce[align->output.cigar_len] = gatk_sw_make_element(state, segment_length + p2);
        align->output.cigar_len += 1;
        alignment_offset = p1 - p2;
    }
    else {  // align->input.overhang_strategy == GATK_SW_OVERHANG_STRATEGY_INDEL || align->input.overhang_strategy ==
            // GATK_SW_OVERHANG_STRATEGY_LEADING_INDEL

        // take care of the actual alignment
        lce[align->output.cigar_len] = gatk_sw_make_element(state, segment_length);
        align->output.cigar_len += 1;
        // take care of overhangs at the beginning of the alignment
        if (p1 > 0) {
            lce[align->output.cigar_len] = gatk_sw_make_element(GATK_SW_STATE_DELETION, p1);
            align->output.cigar_len += 1;
        }
        else if (p2 > 0) {
            lce[align->output.cigar_len] = gatk_sw_make_element(GATK_SW_STATE_INSERTION, p2);
            align->output.cigar_len += 1;
        }

        alignment_offset = 0;
    }

    // res
    for (i = 0, j = align->output.cigar_len / 2; i < j; i++) {
        uint32_t tmp = align->output.cigar[i];
        align->output.cigar[i] = align->output.cigar[align->output.cigar_len - 1 - i];
        align->output.cigar[align->output.cigar_len - 1 - i] = tmp;
    }
    align->output.alignment_offset = alignment_offset;
}

static inline uint32_t gatk_sw_make_element(uint32_t state, int length)
{
    uint8_t op = WORKER_SIGAR_MAX;
    switch (state) {
        case GATK_SW_STATE_MATCH: op = WORKER_SIGAR_STATUS_MATCH; break;
        case GATK_SW_STATE_INSERTION: op = WORKER_SIGAR_STATUS_INSERTION; break;
        case GATK_SW_STATE_DELETION: op = WORKER_SIGAR_STATUS_DELETION; break;
        case GATK_SW_STATE_CLIP: op = WORKER_SIGAR_SOFT_CLIP; break;
        default: return 0;
    }
    return hc_assemble_utils_get_new_cigar(op, length);
}

void hc_assemble_gatk_sw_init(p_gatk_sw_storage align) { (void)align; }

#ifdef __SW_JAVA_TEST__

/**
 * @brief 获取新的单一cigar
 *
 * @param[in]   cigar_op          cigar_op
 * @param[out]  cigar_len         cigar_len
 *
 * @retval uint32_t               cigar
 **/
static inline uint32_t hc_assemble_utils_get_new_cigar(uint8_t cigar_op, uint32_t cigar_len)
{
    uint32_t ret;

    ret = cigar_len << BAM_CIGAR_SHIFT;
    ret |= cigar_op;
    return ret;
}

static void gatk_sw_print_sw(p_gatk_sw_storage align)
{
    int i, j;
    for (i = 0; i < align->cacl_cache.x_len; i++) {
        for (j = 0; j < align->cacl_cache.y_len; j++) {
            printf("%d ", *GATK_SW_GET_SW_POS(align, i, j));
        }
        printf("\n");
    }
    puts("------------------------");
    for (i = 0; i < align->cacl_cache.x_len; i++) {
        for (j = 0; j < align->cacl_cache.y_len; j++) {
            printf("%d ", *GATK_SW_GET_BT_POS(align, i, j));
        }
        printf("\n");
    }
    puts("------------------------");
}

#define BAM_CIGAR_STR  "MIDNSHP=XB"
#define BAM_CIGAR_MASK 0xf
#define BAM_CIGAR_TYPE 0x3C1A7

/*! @abstract Table for converting a CIGAR operator character to BAM_CMATCH etc.
Result is operator code or -1. Be sure to cast the index if it is a plain char:
    int op = bam_cigar_table[(unsigned char) ch];
*/

#define bam_cigar_op(c)    ((c)&BAM_CIGAR_MASK)
#define bam_cigar_oplen(c) ((c) >> BAM_CIGAR_SHIFT)

void main(void)
{
    p_gatk_sw_storage sw;
    const char* cigar_string = "MIDNSHP=XB";
    int i;

    sw = malloc(sizeof(gatk_sw_storage));
    memset(sw, 0, sizeof(gatk_sw_storage));

    sw->input.ref =
        "GTCACGGGACCCGTGATACCAAAACTCGAACATGGGTCGCGACGAAAACGGAACGAGACACTGGGGTCCGTTCGACGGAGTGGAGAGACCCGGTCAAAGGGGTAACATGTCACCACGACGTGTGGGACCG"
        "GGACCGGGGCTCCACCGACCCTCCACCGAGGAGTTTGTCGGCGACAGAGTAGTCACGGGCCACGACCCAGTCCCTAGCTGACTCCGAGACTCGATTGACCCTTTGTGTCACCGGAACCTCCCGACCCCTC"
        "ACAGTACCCCCACCCCTGTCCCTCAGTGGCCAGCGTACACTGACT";
    sw->input.ref_len = strlen(
        "GTCACGGGACCCGTGATACCAAAACTCGAACATGGGTCGCGACGAAAACGGAACGAGACACTGGGGTCCGTTCGACGGAGTGGAGAGACCCGGTCAAAGGGGTAACATGTCACCACGACGTGTGGGACCG"
        "GGACCGGGGCTCCACCGACCCTCCACCGAGGAGTTTGTCGGCGACAGAGTAGTCACGGGCCACGACCCAGTCCCTAGCTGACTCCGAGACTCGATTGACCCTTTGTGTCACCGGAACCTCCCGACCCCTC"
        "ACAGTACCCCCACCCCTGTCCCTCAGTGGCCAGCGTACACTGACT");
    sw->input.alt = "GTCACGGGACCCGTGATACCAAAACCCGAACATGGGTCGCGACG";
    sw->input.alt_len = strlen("GTCACGGGACCCGTGATACCAAAACCCGAACATGGGTCGCGACG");
    sw->input.overhang_strategy = GATK_SW_OVERHANG_STRATEGY_LEADING_INDEL;
    sw->paramates = GATK_SW_STANDARD_NGS;

    printf("ret %d\n offset %d\n", hc_assemble_gatk_sw_align(sw), sw->output.alignment_offset);

    for (i = 0; i < sw->output.cigar_len; i++) {
        printf("%u%c ", bam_cigar_oplen(sw->output.cigar[i]), cigar_string[bam_cigar_op(sw->output.cigar[i])]);
    }
    puts("");
}
#endif