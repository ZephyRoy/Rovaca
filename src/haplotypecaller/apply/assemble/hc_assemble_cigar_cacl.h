#ifndef __PRIVATE_HC_ASSEMBLE_CIGAR_CACL_H__
#define __PRIVATE_HC_ASSEMBLE_CIGAR_CACL_H__

#include "hc_assemble_dijkstra_shortest_path.h"
#include "hc_func.h"
#include "hc_marco.h"
#include "smithwaterman_common.h"
enum HC_ASSEMBLE_CIGAR_CACL_ALLELES_NORMALIZE_T {
    HC_ASSEMBLE_SEQ_FINDER_ALLELES_NORMALIZE_REF = 0,
    HC_ASSEMBLE_SEQ_FINDER_ALLELES_NORMALIZE_SEQ = 1,
    HC_ASSEMBLE_SEQ_FINDER_ALLELES_NORMALIZE_SUM = 2
};

enum HC_ASSEMBLE_CIGAR_CACL_SECTION {
    HC_ASSEMBLE_CIGAR_CACL_LEFT_HARD_CLIP,
    HC_ASSEMBLE_CIGAR_CACL_LEFT_SOFT_CLIP,
    HC_ASSEMBLE_CIGAR_CACL_MIDDLE,
    HC_ASSEMBLE_CIGAR_CACL_RIGHT_SOFT_CLIP,
    HC_ASSEMBLE_CIGAR_CACL_RIGHT_HARD_CLIP
};

/**
 * Represents 0-based integer index range.
 *
 * <p>
 * It represents an integer index range as the pair values:
 * <dl>
 *     <dt>{@link #from}</dt>
 *     <dd>- index of the first element in range (i.e. inclusive).</dd>
 *     <dt>{@link #to}</dt>
 *     <dd>- index of the element following the last element in range (i.e. exclusive).</dd>
 * </dl>
 * </p>
 *
 * <p>
 *     This class is intended to specify a valid index range in arrays or ordered collections.
 * </p>
 *
 * <p>
 *     All instances are constraint so that neither <code>from</code> nor <code>to</code> can
 *     be negative nor <code>from</code> can be larger than <code>to</code>.
 * </p>
 *
 * <p>
 *     You can use {@link #isValidLength(int) isValidFor(length)} to verify that a range instance represents a valid
 *     range for an 0-based indexed object with {@code length} elements.
 * </p>
 */

static inline void hc_assemle_cigar_cacl_range_shift(p_hc_assemble_dijkstra_index_range t, int length)
{
    t->to += length;
    t->from += length;
}

static inline void hc_assemle_cigar_cacl_range_shift_left(p_hc_assemble_dijkstra_index_range t, int length)
{
    hc_assemle_cigar_cacl_range_shift(t, -length);
}

static inline void hc_assemle_cigar_cacl_range_shift_start(p_hc_assemble_dijkstra_index_range t, int length) { t->from += length; }

static inline void hc_assemle_cigar_cacl_range_shift_start_left(p_hc_assemble_dijkstra_index_range t, int length)
{
    hc_assemle_cigar_cacl_range_shift_start(t, -length);
}

static inline void hc_assemle_cigar_cacl_range_shift_end(p_hc_assemble_dijkstra_index_range t, int length) { t->to += length; }

static inline void hc_assemle_cigar_cacl_range_shift_end_left(p_hc_assemble_dijkstra_index_range t, int length)
{
    hc_assemle_cigar_cacl_range_shift_end(t, -length);
}

/**
 * Returns number indexes expanded by this range.
 *
 * @return 0 or greater.
 */
static inline int hc_assemle_cigar_cacl_range_size(p_hc_assemble_dijkstra_index_range t) { return t->to - t->from; }

static inline int hc_assemle_cigar_cacl_length_on_read(uint32_t element)
{
    return consumes_read_bases(bam_cigar_op(element)) ? bam_cigar_oplen(element) : 0;
}

static inline int hc_assemle_cigar_cacl_length_on_reference(uint32_t element)
{
    return consumes_ref_bases(bam_cigar_op(element)) ? bam_cigar_oplen(element) : 0;
}

/** Returns true if the operator is a M, a X or a EQ */
static inline int hc_assemle_cigar_cacl_is_alignment(uint32_t element)
{
    uint32_t t = bam_cigar_op(element);
    return t == WORKER_SIGAR_STATUS_MATCH || t == WORKER_SIGAR_S_MISMATCH || t == WORKER_SIGAR_S_MATCH;
}

static int hc_assemle_cigar_cacl_cigar_get_latest_indel(uint32_t* cigar, uint32_t cigar_len);
static int hc_assemle_cigar_cacl_cigar_get_part_ref_length(uint32_t* cigar, uint32_t cigar_len, uint32_t all_len, int* all_ref);
static void hc_assemle_cigar_cacl_run_sw(p_hc_apply apply, p_assemble_dijkstra_one_path_node path);
static int hc_assemle_cigar_cacl_is_sw_failure(p_lib_sw_avx alignment);
static void hc_assemle_cigar_cacl_left_align_indels(p_assemble_haplotype haplot, int readStart, p_hc_apply apply);
static int hc_assemle_cigar_cacl_is_indel(uint32_t t);
static void hc_assemle_cigar_cacl_trim_cigar_by_bases(p_assemble_dijkstra_one_path_node st, int start, int end);
static int hc_assemle_cigar_cacl_get_min_bound(p_hc_assemble_dijkstra_index_range first, p_hc_assemble_dijkstra_index_range second);
static int hc_assemle_cigar_cacl_make_and_record_deletions_removed_result(p_hc_assemble_dijkstra_out_path_storage st);
static void hc_assemle_cigar_cacl_add_cigar(uint32_t element, p_hc_assemble_dijkstra_out_path_storage st);
static void hc_assemle_cigar_cacl_advance_section_and_validate_cigar_order(uint32_t o, p_hc_assemble_dijkstra_out_path_storage st);
static inline int hc_assemle_cigar_cacl_last_two_elements_were_deletion_and_insertion(p_hc_assemble_dijkstra_out_path_storage st);
static inline int hc_assemle_cigar_cacl_is_clipping(uint32_t cigar_op);
static void hc_assemle_cigar_cacl_normalize_alleles(p_assemble_haplotype haplot, p_assemble_dijkstra_path haplot_st, int maxShift);
static inline int hc_assemle_cigar_cacl_last_base_on_right_is_same(p_assemble_haplotype haplot, p_assemble_dijkstra_path haplot_st);
static inline int hc_assemle_cigar_cacl_first_base_on_left_is_same(p_assemble_haplotype haplot, p_assemble_dijkstra_path haplot_st);
static inline int hc_assemle_cigar_cacl_next_base_on_left_is_same(p_assemble_haplotype haplot, p_assemble_dijkstra_path haplot_st);
static void hc_assemle_cigar_cacl_init_cigar_add(p_assemble_haplotype haplot);

#endif