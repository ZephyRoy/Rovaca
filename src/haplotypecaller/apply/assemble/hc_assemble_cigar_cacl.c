/**
 * @file hc_assemble_cigar_cacl.c
 * @author Cao Ce (caoce@genomics.cn)
 * @brief Cigar修复
 * @version 0.1
 * @date 2021-12-28
 *
 * @copyright Copyright (c) 2021 Rovaca SDK
 *
 */

#include "hc_assemble_cigar_cacl.h"

#include "hc_apply.h"
#include "hc_assemble.h"
#include "hc_func.h"

/**
 * Calculate the cigar elements for this path against the reference sequence.
 *
 * This assumes that the reference and alt sequence are haplotypes derived from a de Bruijn graph or SeqGraph and have the same
 * ref source and ref sink vertices.  That is, the alt sequence start and end are assumed anchored to the reference start and end, which
 * occur at the ends of the padded assembly region.  Hence, unlike read alignment, there is no concept of a start or end coordinate here.
 * Furthermore, it is important to note that in the rare case that the alt cigar begins or ends with a deletion, we must keep the leading
 * or trailing deletion in order to maintain the original reference span of the alt haplotype.  This can occur, for example, when the ref
 * haplotype starts with N repeats of a long sequence and the alt haplotype starts with N-1 repeats.
 *
 * @param aligner
 * @param refSeq the reference sequence that all of the bases in this path should align to
 * @return a Cigar mapping this path to refSeq, or null if no reasonable alignment could be found
 */
void hc_assemle_cigar_cacl_calculate_cigar(p_hc_apply apply, p_assemble_dijkstra_one_path_node altSeq)
{
    p_assemble_dijkstra_path graph = &apply->read_graph.dijkstra_path_finder;
    int baseStart, baseEnd, trimmedLeadingDeletions, trimmedTrailingDeletions;

    if (__glibc_unlikely(altSeq->st->seq_len == 0)) {
        // horrible edge case from the unit tests, where this path has no bases
        altSeq->st->cigar[0] = hc_assemble_utils_get_new_cigar(WORKER_SIGAR_STATUS_DELETION, apply->ref_len);
        altSeq->st->cigar_len = 1;
        return;
    }

    // Note: this is a performance optimization.
    //  If two strings are equal (a O(n) check) then it's trivial to get CIGAR for them.
    //  Furthermore, if their lengths are equal and their element-by-element comparison yields two or fewer mismatches
    //  it's also a trivial M-only CIGAR, because in order to have equal length one would need at least one insertion and
    //  one deletion, in which case two substitutions is a better alignment.
    if (altSeq->st->seq_len == apply->read_thread_graph->origin_ref_len + strlen(HC_ASSEMBLE_SW_PAD_STRING) * 2) {
        int mismatchCount = 0;
        uint32_t n;

        for (n = 0; n < apply->read_thread_graph->origin_ref_len && mismatchCount <= 2; n++) {
            mismatchCount += (altSeq->st->seq[n + strlen(HC_ASSEMBLE_SW_PAD_STRING)] == apply->read_thread_graph->origin_ref[n] ? 0 : 1);
        }
        if (mismatchCount <= 2) {
            altSeq->st->cigar[0] = hc_assemble_utils_get_new_cigar(WORKER_SIGAR_STATUS_MATCH, altSeq->st->seq_len);
            graph->ref_len = strlen(HC_ASSEMBLE_SW_PAD_STRING) * 2 + apply->read_thread_graph->origin_ref_len;
            altSeq->st->cigar_len = 1;
            goto repair_cigar;
        }
    }

    graph->ref_len = strlen(HC_ASSEMBLE_SW_PAD_STRING) * 2 + apply->read_thread_graph->origin_ref_len;
    graph->ref = apply->sw_avx->seq1;

    hc_assemle_cigar_cacl_run_sw(apply, altSeq);

    if (hc_assemle_cigar_cacl_is_sw_failure(apply->sw_avx)) {
        altSeq->st->cigar_len = 0;
        return;
    }
    else {
        altSeq->st->cigar_len = apply->sw_avx->cigar_count;
        memcpy(altSeq->st->cigar, apply->sw_avx->cigar, altSeq->st->cigar_len * sizeof(uint32_t));
    }

repair_cigar:

    // cut off the padding bases
    baseStart = strlen(HC_ASSEMBLE_SW_PAD_STRING);
    baseEnd = altSeq->st->seq_len - strlen(HC_ASSEMBLE_SW_PAD_STRING) - 1;  // -1 because it's inclusive
    hc_assemle_cigar_cacl_trim_cigar_by_bases(altSeq, baseStart, baseEnd);

    // leading deletion removed by cigar trimming shift the alignment start to the right
    trimmedLeadingDeletions = altSeq->st->leadingDeletionBasesRemoved;
    // trailing deletions should be kept in order to left-align
    trimmedTrailingDeletions = altSeq->st->trailingDeletionBasesRemoved;

    if (trimmedTrailingDeletions > 0) {
        altSeq->st->cigar_back[altSeq->st->cigar_back_len] =
            hc_assemble_utils_get_new_cigar(WORKER_SIGAR_STATUS_DELETION, trimmedTrailingDeletions);
    }

    hc_assemle_cigar_cacl_left_align_indels(altSeq, trimmedLeadingDeletions, apply);

    // we must account for possible leading deletions removed when trimming the padding and when left-aligning
    // trailing deletions removed when trimming have already been restored for left-alignment, but left-alingment may have removed them
    // again.
    int totalLeadingDeletionsRemoved = trimmedLeadingDeletions + altSeq->st->leadingDeletionBasesRemoved;
    int totalTrailingDeletionsRemoved = altSeq->st->trailingDeletionBasesRemoved;

    if (totalLeadingDeletionsRemoved == 0 && totalTrailingDeletionsRemoved == 0) {
        memcpy(altSeq->st->cigar, altSeq->st->cigar_back, altSeq->st->cigar_back_len * sizeof(uint32_t));
        altSeq->st->cigar_len = altSeq->st->cigar_back_len;
        return;
    }
    else {
        altSeq->st->cigar_len = 0;
        if (totalLeadingDeletionsRemoved > 0) {
            altSeq->st->cigar[0] = hc_assemble_utils_get_new_cigar(WORKER_SIGAR_STATUS_DELETION, totalLeadingDeletionsRemoved);
            altSeq->st->cigar_len = 1;
        }
        memcpy(&altSeq->st->cigar[1], altSeq->st->cigar_back, altSeq->st->cigar_back_len * sizeof(uint32_t));
        altSeq->st->cigar_len += altSeq->st->cigar_back_len;
        if (totalTrailingDeletionsRemoved > 0) {
            altSeq->st->cigar[altSeq->st->cigar_len] =
                hc_assemble_utils_get_new_cigar(WORKER_SIGAR_STATUS_DELETION, totalTrailingDeletionsRemoved);
            altSeq->st->cigar_len += 1;
        }
    }
}

/**
 * Takes the alignment of the read sequence <code>readSeq</code> to the reference sequence <code>refSeq</code>
 * starting at 0-based position <code>readStart</code> on the <code>ref</code> and specified by its <code>cigar</code>.
 * <p/>
 * If the alignment has one or more indels, this method attempts to move them left across a stretch of repetitive bases.
 * For instance, if the original cigar specifies that (any) one AT is deleted from a repeat sequence TATATATA, the output
 * cigar will always mark the leftmost AT as deleted. If there is no indel in the original cigar or if the indel position
 * is determined unambiguously (i.e. inserted/deleted sequence is not repeated), the original cigar is returned.
 *
 * Soft-clipped bases in the cigar are presumed to correspond to bases in the byte[] of read sequence.  That is, this method
 * assumes that inputs are precise about the distinction between hard clips (removed from the read sequence) and soft clips
 * (kept in the read sequence but not aligned).  For example, with the inputs {cigar: 2S2M2I, read sequence: TTAAAA, ref sequence: GGAA,
 * read start: 2} the method lines up the AAAA (2M2I) of the read with the AA of the ref and left-aligns the indel to yield a cigar of
 * 2S2I2M.
 *
 * @param cigar     structure of the original alignment
 * @param ref    reference sequence the read is aligned to
 * @param read   read sequence
 * @param readStart  0-based position on ref of the first aligned base in the read sequence
 * @return a non-null cigar, in which the indels are guaranteed to be placed at the leftmost possible position across a repeat (if any)
 */
static void hc_assemle_cigar_cacl_left_align_indels(p_assemble_haplotype haplot, int readStart, p_hc_apply apply)
{
    int latest_indel, refLength, newMatchOnLeftDueToTrimming, remainingBasesOnLeft;
    uint32_t* resultRightToLeft = haplot->st->cigar;
    uint32_t* resultRightToLeftLen = &haplot->st->cigar_len;
    uint32_t n = 0;
    *resultRightToLeftLen = 0;

    // we need reference bases from the start of the read to the rightmost indel
    latest_indel = hc_assemle_cigar_cacl_cigar_get_latest_indel(haplot->st->cigar_back, haplot->st->cigar_back_len);
    if (latest_indel == -1) {
        haplot->st->leadingDeletionBasesRemoved = 0;
        haplot->st->trailingDeletionBasesRemoved = 0;
        return;
    }
    int r =
        hc_assemle_cigar_cacl_cigar_get_part_ref_length(haplot->st->cigar_back, latest_indel + 1, haplot->st->cigar_back_len, &refLength);
    if (readStart + r <= (int)apply->ref_len) {
        return;
    }

    // at this point, we are one base past the end of the read.  Now we traverse the cigar from right to left
    haplot->st->refIndelRange.from = readStart + refLength;
    haplot->st->refIndelRange.to = readStart + refLength;
    haplot->st->readIndelRange.from = haplot->st->seq_len - strlen(HC_ASSEMBLE_SW_PAD_STRING) * 2;
    haplot->st->readIndelRange.to = haplot->st->seq_len - strlen(HC_ASSEMBLE_SW_PAD_STRING) * 2;

    for (n = haplot->st->cigar_back_len - 1; n != 0; n--) {
        uint32_t element = haplot->st->cigar_back[n];
        uint32_t element_op = bam_cigar_op(element);
        uint32_t element_opl = bam_cigar_oplen(element);

        // if it's an indel, just accumulate the read and ref bases consumed.  We won't shift the indel until we hit an alignment
        // block or the read start.
        if (hc_assemle_cigar_cacl_is_indel(element)) {
            hc_assemle_cigar_cacl_range_shift_start_left(&haplot->st->refIndelRange, hc_assemle_cigar_cacl_length_on_reference(element));
            hc_assemle_cigar_cacl_range_shift_start_left(&haplot->st->readIndelRange, hc_assemle_cigar_cacl_length_on_read(element));
        }
        else if (hc_assemle_cigar_cacl_range_size(&haplot->st->refIndelRange) == 0 &&
                 hc_assemle_cigar_cacl_range_size(&haplot->st->readIndelRange) == 0) {
            // no indel, just add the cigar element to the result
            resultRightToLeft[*resultRightToLeftLen] = element;
            *resultRightToLeftLen += 1;
            hc_assemle_cigar_cacl_range_shift_left(&haplot->st->refIndelRange, hc_assemle_cigar_cacl_length_on_reference(element));
            hc_assemle_cigar_cacl_range_shift_left(&haplot->st->readIndelRange, hc_assemle_cigar_cacl_length_on_read(element));
        }
        else {
            // there's an indel that we need to left-align
            // we can left-align into match cigar elements but not into clips
            int maxShift = hc_assemle_cigar_cacl_is_alignment(element) ? bam_cigar_oplen(element) : 0;

            hc_assemle_cigar_cacl_normalize_alleles(haplot, &apply->read_graph.dijkstra_path_finder, maxShift);

            // Pair<Integer, Integer> shifts = hc_assemle_cigar_cacl_normalize_alleles(Arrays.asList(ref, read),
            // Arrays.asList(refIndelRange, readIndelRange), maxShift, true);
            //  startShift, endShift

            // account for new match alignments on the right due to left-alignment
            resultRightToLeft[*resultRightToLeftLen] = hc_assemble_utils_get_new_cigar(WORKER_SIGAR_STATUS_MATCH, haplot->st->end_shift);
            *resultRightToLeftLen += 1;

            // emit if we didn't go all the way to the start of an alignment block OR we have reached clips OR we have reached the start of
            // the cigar we may have actually shifted right to make the alleles parsimonious
            newMatchOnLeftDueToTrimming = haplot->st->start_shift < 0 ? -haplot->st->start_shift : 0;
            remainingBasesOnLeft = haplot->st->start_shift < 0 ? element_opl : (element_opl - haplot->st->start_shift);

            // some of this alignment block remains after left-alignment -- emit the indel
            if (n == 0 || haplot->st->start_shift < maxShift || !hc_assemle_cigar_cacl_is_alignment(element)) {
                resultRightToLeft[*resultRightToLeftLen] = hc_assemble_utils_get_new_cigar(
                    WORKER_SIGAR_STATUS_DELETION, hc_assemle_cigar_cacl_range_size(&haplot->st->refIndelRange));
                *resultRightToLeftLen += 1;
                resultRightToLeft[*resultRightToLeftLen] = hc_assemble_utils_get_new_cigar(
                    WORKER_SIGAR_STATUS_INSERTION, hc_assemle_cigar_cacl_range_size(&haplot->st->readIndelRange));
                *resultRightToLeftLen += 1;

                // ref indel range is now empty and points to start of left-aligned indel
                hc_assemle_cigar_cacl_range_shift_end_left(&haplot->st->refIndelRange,
                                                           hc_assemle_cigar_cacl_range_size(&haplot->st->refIndelRange));
                // read indel range is now empty and points to start of left-aligned indel
                hc_assemle_cigar_cacl_range_shift_end_left(&haplot->st->readIndelRange,
                                                           hc_assemle_cigar_cacl_range_size(&haplot->st->readIndelRange));

                hc_assemle_cigar_cacl_range_shift_left(
                    &haplot->st->refIndelRange, newMatchOnLeftDueToTrimming + (consumes_ref_bases(element_op) ? remainingBasesOnLeft : 0));
                hc_assemle_cigar_cacl_range_shift_left(
                    &haplot->st->readIndelRange,
                    newMatchOnLeftDueToTrimming + (consumes_read_bases(element_op) ? remainingBasesOnLeft : 0));
                // now read and ref indel ranges are empty and point to end of element preceding this block
            }
            resultRightToLeft[*resultRightToLeftLen] =
                hc_assemble_utils_get_new_cigar(WORKER_SIGAR_STATUS_MATCH, newMatchOnLeftDueToTrimming);
            *resultRightToLeftLen += 1;
            resultRightToLeft[*resultRightToLeftLen] = hc_assemble_utils_get_new_cigar(element_opl, remainingBasesOnLeft);
            *resultRightToLeftLen += 1;
        }
    }

    // account for any indels at the start of the cigar that weren't processed because they have no adjacent non-indel element to the left
    resultRightToLeft[*resultRightToLeftLen] =
        hc_assemble_utils_get_new_cigar(WORKER_SIGAR_STATUS_DELETION, hc_assemle_cigar_cacl_range_size(&haplot->st->refIndelRange));
    *resultRightToLeftLen += 1;
    resultRightToLeft[*resultRightToLeftLen] =
        hc_assemble_utils_get_new_cigar(WORKER_SIGAR_STATUS_INSERTION, hc_assemle_cigar_cacl_range_size(&haplot->st->readIndelRange));
    *resultRightToLeftLen += 1;

    // rebuild
    hc_assemble_utils_reverse_stack(resultRightToLeft, *resultRightToLeftLen);
    hc_assemle_cigar_cacl_init_cigar_add(haplot);
    haplot->st->cigar_back_len = 0;
    for (n = 0; n < haplot->st->cigar_len; n++) {
        hc_assemle_cigar_cacl_add_cigar(haplot->st->cigar[n], haplot->st);
    }
    hc_assemle_cigar_cacl_make_and_record_deletions_removed_result(haplot->st);
}

static void hc_assemle_cigar_cacl_init_cigar_add(p_assemble_haplotype haplot)
{
    haplot->st->trailingDeletionBasesRemoved = 0;
    haplot->st->leadingDeletionBasesRemoved = 0;
    haplot->st->section = HC_ASSEMBLE_CIGAR_CACL_LEFT_HARD_CLIP;
}

static int hc_assemle_cigar_cacl_get_min_bound(p_hc_assemble_dijkstra_index_range first, p_hc_assemble_dijkstra_index_range second)
{
    int first_size = hc_assemle_cigar_cacl_range_size(first);
    int second_szie = hc_assemle_cigar_cacl_range_size(second);

    return assemble_min(first_size, second_szie);
}

/**
 *  Example usage:  reference = GAAT, read = GAAAT (insertion of one A) and we initially consider the insertion of the A to occur before
 *  the T.  Thus the reference range of this allele is [3,3) (no bases) and the read range is [3,4).  This will be left-aligned so that
 *  the insertion occurs after the G, so that the ranges become [1,1) and [1,2) and the returned shifts are 2 bases for both the start and
 * end of the range.
 *
 *  If the given allele ranges are not parsimonious, for example [3,4) and [3,5) in the above example to include the common T in both
 * alleles, the resulting ranges will be shifted by different amounts.  In this case, the shifts are 2 bases in the front and 3 at the end.
 *
 *  Note that we use the convention that the ref allele in the case of an alt insertion, or the alt allele in case of a deletion, is
 * represented by [n,n) where n is the last aligned coordinate before the indel.  This makes sense when you think in terms of alignment
 * CIGARs: suppose for example we have a 5M5I5M read with start 100.  Then the match bases are [100,105) on the ref and [0,5) on the read
 * and the inserted bases are [5,10) on the read and [5,5) on the reference.
 *
 * @param sequences bases of sequences containing different alleles -- could be reference, a haplotype, a read, or subsequences thereof
 * @param bounds    initial ranges (inclusive start, exclusive end) of alleles in same order as {@code sequences}
 *                  Note that this method adjusts these ranges as a side effect
 * @param maxShift  maximum allowable shift left.  This may be less than the amount demanded by the array bounds.  For example, when
 *                  left-aligning a read with multiple indels, we don't want to realign one indel past another (if they "collide" we merge
 *                  them into a single indel and continue -- see {@link AlignmentUtils::hc_assemle_cigar_cacl_left_align_indels}
 * @param trim      If true, remove redundant shared bases at the start and end of alleles
 * @return          The number of bases the alleles were shifted left such that they still represented the same event.
 */
static void hc_assemle_cigar_cacl_normalize_alleles(p_assemble_haplotype haplot, p_assemble_dijkstra_path haplot_st, int maxShift)
{
    int startShift = 0;
    int endShift = 0;

    // consume any redundant shared bases at the end of the alleles
    int minSize = hc_assemle_cigar_cacl_get_min_bound(&haplot->st->readIndelRange, &haplot->st->refIndelRange);

    while (minSize > 0 && hc_assemle_cigar_cacl_last_base_on_right_is_same(haplot, haplot_st)) {
        hc_assemle_cigar_cacl_range_shift_end_left(&haplot->st->refIndelRange, 1);
        hc_assemle_cigar_cacl_range_shift_end_left(&haplot->st->readIndelRange, 1);
        minSize--;
        endShift++;
    }

    while (minSize > 0 && hc_assemle_cigar_cacl_first_base_on_left_is_same(haplot, haplot_st)) {
        hc_assemle_cigar_cacl_range_shift_start(&haplot->st->refIndelRange, 1);
        hc_assemle_cigar_cacl_range_shift_start(&haplot->st->readIndelRange, 1);
        minSize--;
        startShift--;
    }

    // we shift left as long as the last bases on the right are equal among all sequences and the next bases on the left are all equal.
    // if a sequence is empty (eg the reference relative to an insertion alt allele) the last base on the right is the next base on the left
    while (startShift < maxShift && hc_assemle_cigar_cacl_next_base_on_left_is_same(haplot, haplot_st) &&
           hc_assemle_cigar_cacl_last_base_on_right_is_same(haplot, haplot_st)) {
        hc_assemle_cigar_cacl_range_shift_left(&haplot->st->refIndelRange, 1);
        hc_assemle_cigar_cacl_range_shift_left(&haplot->st->readIndelRange, 1);
        startShift++;
        endShift++;
    }

    haplot->st->start_shift = startShift;
    haplot->st->end_shift = endShift;
}

// do all sequences share a common base at the end of the given index range
static inline int hc_assemle_cigar_cacl_last_base_on_right_is_same(p_assemble_haplotype haplot, p_assemble_dijkstra_path haplot_st)
{
    uint8_t* ref = haplot_st->ref + strlen(HC_ASSEMBLE_SW_PAD_STRING);
    uint8_t* seq = haplot->st->seq + strlen(HC_ASSEMBLE_SW_PAD_STRING);

    if (seq[haplot->st->readIndelRange.to - 1] != ref[haplot->st->refIndelRange.to - 1]) {
        return 0;
    }

    return 1;
}

// do all sequences share a common first base
static inline int hc_assemle_cigar_cacl_first_base_on_left_is_same(p_assemble_haplotype haplot, p_assemble_dijkstra_path haplot_st)
{
    uint8_t* ref = haplot_st->ref + strlen(HC_ASSEMBLE_SW_PAD_STRING);
    uint8_t* seq = haplot->st->seq + strlen(HC_ASSEMBLE_SW_PAD_STRING);

    if (seq[haplot->st->readIndelRange.from] != ref[haplot->st->refIndelRange.from]) {
        return 0;
    }

    return 1;
}

// do all sequences share a common base just before the given index ranges
static inline int hc_assemle_cigar_cacl_next_base_on_left_is_same(p_assemble_haplotype haplot, p_assemble_dijkstra_path haplot_st)
{
    uint8_t* ref = haplot_st->ref + strlen(HC_ASSEMBLE_SW_PAD_STRING);
    uint8_t* seq = haplot->st->seq + strlen(HC_ASSEMBLE_SW_PAD_STRING);

    if (seq[haplot->st->readIndelRange.from - 1] != ref[haplot->st->refIndelRange.from - 1]) {
        return 0;
    }

    return 1;
}

/** Returns true if the operator is a Insertion or Deletion operator */
static int hc_assemle_cigar_cacl_is_indel(uint32_t t) { return t == WORKER_SIGAR_STATUS_INSERTION || t == WORKER_SIGAR_STATUS_DELETION; }

static int hc_assemle_cigar_cacl_cigar_get_latest_indel(uint32_t* cigar, uint32_t cigar_len)
{
    int ret_pos = -1;
    uint32_t i = 0;
    uint32_t opr, opl, ret_length = 0;

    for (; i < cigar_len; i++) {
        opr = bam_cigar_op(cigar[i]);
        if (opr == WORKER_SIGAR_STATUS_INSERTION || opr == WORKER_SIGAR_STATUS_DELETION) {
            opl = bam_cigar_oplen(cigar[i]);
            if (opl > ret_length) {
                ret_length = opl;
                ret_pos = i;
            }
        }
    }
    return ret_pos;
}

static int hc_assemle_cigar_cacl_cigar_get_part_ref_length(uint32_t* cigar, uint32_t cigar_len, uint32_t all_len, int* all_ref)
{
    int ret_len = 0;
    uint32_t i = 0;
    uint32_t opr, opl;

    *all_ref = 0;

    for (; i < all_len; i++) {
        opr = bam_cigar_op(cigar[i]);
        opl = bam_cigar_oplen(cigar[i]);

        if (consumes_ref_bases(opr)) {
            if (i <= cigar_len + 1) {
                ret_len += opl;
            }
            *all_ref += opl;
        }
    }
    return ret_len;
}

static void hc_assemle_cigar_cacl_run_sw(p_hc_apply apply, p_assemble_dijkstra_one_path_node path)
{
    p_lib_sw_avx sw = apply->sw_avx;

    sprintf((char*)sw->seq1, "%s%.*s%s", HC_ASSEMBLE_SW_PAD_STRING, apply->read_thread_graph->origin_ref_len,
            apply->read_thread_graph->origin_ref, HC_ASSEMBLE_SW_PAD_STRING);
    sw->len1 = apply->read_graph.dijkstra_path_finder.ref_len;
    memcpy(sw->seq2, path->st->seq, path->st->seq_len);
    sw->len2 = path->st->seq_len;

    sw->avx_function(sw);
}

/**
 * Make sure that the SW didn't fail in some terrible way, and throw exception if it did
 */
static int hc_assemle_cigar_cacl_is_sw_failure(p_lib_sw_avx alignment)
{
    uint32_t* cigar = alignment->cigar;
    uint32_t this_cigar, cigar_count = alignment->cigar_count;
    register uint32_t i;

    // check that the alignment starts at the first base, which it should given the padding
    if (alignment->alignment_offset > 0) {
        return 1;
    }

    // check that we aren't getting any S operators (which would be very bad downstream)
    for (i = 0; i < cigar_count; i++) {
        this_cigar = bam_cigar_op(cigar[i]);

        // soft clips at the end of the alignment are really insertions
        if (this_cigar == WORKER_SIGAR_SOFT_CLIP) {
            return 1;
        }
    }

    return 0;
}

/**
 * Workhorse for hc_assemle_cigar_cacl_trim_cigar_by_bases and trimCigarByReference
 *
 * @param cigar a non-null Cigar to trim down
 * @param start Where should we start keeping bases in the cigar (inclusive)?  The first position is 0
 * @param end Where should we stop keeping bases in the cigar (inclusive)?  The maximum value is cigar.getLength() - 1
 * @param byReference should start and end be interpreted as position in the reference or the read to trim to/from?
 * @return a non-null cigar
 */
static void hc_assemle_cigar_cacl_trim_cigar_by_bases(p_assemble_dijkstra_one_path_node st, int start, int end)
{
    // these variables track the inclusive start and exclusive end of the current cigar element in reference (if byReference) or read
    // (otherwise) coordinates
    int elementStart;    // inclusive
    int elementEnd = 0;  // exclusive -- start of next element
    uint32_t overlapLength;
    uint32_t i;
    uint32_t* cigar = st->st->cigar;

    hc_assemle_cigar_cacl_init_cigar_add(st);
    for (i = 0; i < st->st->cigar_len; i++) {
        elementStart = elementEnd;
        elementEnd = elementStart + hc_assemle_cigar_cacl_length_on_read(cigar[i]);

        // we are careful to include zero-length elements at both ends, that is, elements with elementStart == elementEnd == start and
        // elementStart == elementEnd == end + 1
        if (elementEnd < start || (elementEnd == start && elementStart < start)) {
            continue;
        }
        else if (elementStart > end && elementEnd > end + 1) {
            break;
        }

        overlapLength = elementEnd == elementStart ? bam_cigar_oplen(cigar[i])
                                                   : (uint32_t)(assemble_min(end + (int)1, elementEnd) - assemble_max(start, elementStart));
        hc_assemle_cigar_cacl_add_cigar(hc_assemble_utils_get_new_cigar(bam_cigar_op(cigar[i]), overlapLength), st->st);
    }

    hc_assemle_cigar_cacl_make_and_record_deletions_removed_result(st->st);
}

/**
 * Return a Result object containing the output of make() as well as the number of leading and trailing deletion bases
 * removed relative to the cigar elements that were add()ed.  This is very useful when in addition to transforming a cigar we must also
 * keep track of an alignment start or end.
 */
static int hc_assemle_cigar_cacl_make_and_record_deletions_removed_result(p_hc_assemble_dijkstra_out_path_storage st)
{
    if (st->section == HC_ASSEMBLE_CIGAR_CACL_LEFT_SOFT_CLIP && bam_cigar_op(st->cigar_back[0]) == WORKER_SIGAR_SOFT_CLIP) {
        return 0;
    }

    st->trailingDeletionBasesRemovedInMake = 0;
    if (bam_cigar_op(st->lastOperator) == WORKER_SIGAR_STATUS_DELETION) {
        st->trailingDeletionBasesRemovedInMake = bam_cigar_oplen(st->cigar_back[st->cigar_back_len - 1]);
        st->cigar_back_len -= 1;
    }
    else if (hc_assemle_cigar_cacl_last_two_elements_were_deletion_and_insertion(st)) {
        st->trailingDeletionBasesRemovedInMake = bam_cigar_oplen(st->cigar_back[st->cigar_back_len - 2]);
        st->cigar_back[st->cigar_back_len - 2] = st->cigar_back[st->cigar_back_len - 1];
        st->cigar_back_len -= 1;
    }
    if (st->cigar_back_len == 0) {
        return 0;
    }
    else {
        return 1;
    }
}

static void hc_assemle_cigar_cacl_add_cigar(uint32_t element, p_hc_assemble_dijkstra_out_path_storage st)
{
    uint32_t element_op = bam_cigar_op(element);
    uint32_t element_op_len = bam_cigar_oplen(element);
    uint32_t lastOperator_op = bam_cigar_op(st->lastOperator);

    if (element_op_len == 0) {
        return;
    }

    // skip a deletion after clipping ie at the beginning of the read
    // note the edge case of a deletion following a leading insertion, which we also skip
    if (element_op == WORKER_SIGAR_STATUS_DELETION) {
        if (lastOperator_op == WORKER_SIGAR_MAX || hc_assemle_cigar_cacl_is_clipping(lastOperator_op)) {
            st->leadingDeletionBasesRemoved += element_op_len;
            return;
        }
        else if (lastOperator_op == WORKER_SIGAR_STATUS_INSERTION &&
                 (st->cigar_back_len == 1 || hc_assemle_cigar_cacl_is_clipping(bam_cigar_op(st->cigar_back[st->cigar_back_len - 2])))) {
            st->leadingDeletionBasesRemoved += element_op_len;
            return;
        }
    }

    hc_assemle_cigar_cacl_advance_section_and_validate_cigar_order(element_op, st);

    // merge consecutive elements with the same element_op
    if (element_op == lastOperator_op) {
        if (st->cigar_back_len) {
            st->cigar_back[st->cigar_back_len - 1] =
                hc_assemble_utils_get_new_cigar(element_op, bam_cigar_oplen(st->cigar_back[st->cigar_back_len - 1]) + element_op_len);
        }
        else {
            st->cigar_back[0] = hc_assemble_utils_get_new_cigar(element_op, element_op_len);
            st->cigar_back_len += 1;
        }
    }
    else {
        if (lastOperator_op == WORKER_SIGAR_MAX) {
            st->cigar_back[st->cigar_back_len] = element;
            st->cigar_back_len += 1;
            st->lastOperator = element;
        }
        else if (hc_assemle_cigar_cacl_is_clipping(element_op)) {
            // if we have just started clipping on the right and realize the last element_op was a deletion, remove it
            // if we have just started clipping on the right and the last two operators were a deletion and insertion, remove the deletion
            if (!consumes_read_bases(lastOperator_op) && !hc_assemle_cigar_cacl_is_clipping(lastOperator_op)) {
                st->trailingDeletionBasesRemoved += bam_cigar_oplen(st->cigar_back[st->cigar_back_len - 1]);
                st->cigar_back[st->cigar_back_len - 1] = element;
                st->lastOperator = element;
            }
            else if (hc_assemle_cigar_cacl_last_two_elements_were_deletion_and_insertion(st)) {
                st->trailingDeletionBasesRemoved += bam_cigar_oplen(st->cigar_back[st->cigar_back_len - 2]);

                st->cigar_back[st->cigar_back_len - 2] = st->cigar_back[st->cigar_back_len - 1];
                st->cigar_back[st->cigar_back_len - 1] =
                    hc_assemble_utils_get_new_cigar(element_op, bam_cigar_oplen(st->cigar_back[st->cigar_back_len - 1]));
            }
            else {
                st->cigar_back[st->cigar_back_len] = element;
                st->cigar_back_len += 1;
                st->lastOperator = element;
            }
        }
        else if (element_op == WORKER_SIGAR_STATUS_DELETION && lastOperator_op == WORKER_SIGAR_STATUS_INSERTION) {
            // The order of deletion and insertion elements is arbitrary, so to standardize we shift deletions to the left
            // that is, we place the deletion before the insertion and shift the insertion right
            // if the element before the insertion is another deletion, we merge in the new deletion
            // note that the last element_op remains an insertion
            int size = st->cigar_back_len;
            if (size > 1 && bam_cigar_op(st->cigar_back[size - 2]) == WORKER_SIGAR_STATUS_DELETION) {
                st->cigar_back[size - 2] = hc_assemble_utils_get_new_cigar(WORKER_SIGAR_STATUS_DELETION,
                                                                           bam_cigar_oplen(st->cigar_back[size - 2]) + element_op_len);
            }
            else {
                st->cigar_back[st->cigar_back_len] = st->cigar_back[st->cigar_back_len - 1];
                st->cigar_back[st->cigar_back_len - 1] = element;
                st->cigar_back_len += 1;
            }
        }
        else {
            st->cigar_back[st->cigar_back_len] = element;
            st->cigar_back_len += 1;
            st->lastOperator = element_op;
        }
    }
}

static inline int hc_assemle_cigar_cacl_last_two_elements_were_deletion_and_insertion(p_hc_assemble_dijkstra_out_path_storage st)
{
    return st->lastOperator == WORKER_SIGAR_STATUS_INSERTION && st->cigar_back_len > 1 &&
           bam_cigar_op(st->cigar_back[st->cigar_back_len - 2]) == WORKER_SIGAR_STATUS_DELETION;
}

// validate that cigar structure is hard clip, soft clip, unclipped, soft clip, hard clip
static void hc_assemle_cigar_cacl_advance_section_and_validate_cigar_order(uint32_t o, p_hc_assemble_dijkstra_out_path_storage st)
{
    if (o == WORKER_SIGAR_HARD_CLIP) {
        if (st->section == HC_ASSEMBLE_CIGAR_CACL_LEFT_SOFT_CLIP || st->section == HC_ASSEMBLE_CIGAR_CACL_MIDDLE ||
            st->section == HC_ASSEMBLE_CIGAR_CACL_RIGHT_SOFT_CLIP) {
            st->section = HC_ASSEMBLE_CIGAR_CACL_RIGHT_HARD_CLIP;
        }
    }
    else if (o == WORKER_SIGAR_SOFT_CLIP) {
        if (st->section == HC_ASSEMBLE_CIGAR_CACL_LEFT_HARD_CLIP) {
            st->section = HC_ASSEMBLE_CIGAR_CACL_LEFT_SOFT_CLIP;
        }
        else if (st->section == HC_ASSEMBLE_CIGAR_CACL_MIDDLE) {
            st->section = HC_ASSEMBLE_CIGAR_CACL_RIGHT_SOFT_CLIP;
        }
    }
    else {
        if (st->section == HC_ASSEMBLE_CIGAR_CACL_LEFT_HARD_CLIP || st->section == HC_ASSEMBLE_CIGAR_CACL_LEFT_SOFT_CLIP) {
            st->section = HC_ASSEMBLE_CIGAR_CACL_MIDDLE;
        }
    }
}

/** Returns true if the operator is a clipped (hard or soft) operator */
static inline int hc_assemle_cigar_cacl_is_clipping(uint32_t cigar_op)
{
    return cigar_op == WORKER_SIGAR_SOFT_CLIP || cigar_op == WORKER_SIGAR_HARD_CLIP;
}
