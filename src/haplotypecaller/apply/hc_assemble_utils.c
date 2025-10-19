/**
 * @file hc_apply_utils.c
 * @author Cao Ce (caoce@genomics.cn)
 * @brief Apply工具
 * @version 0.1
 * @date 2021-05-18
 *
 * @copyright Copyright (c) 2021 ROVACA SDK
 *
 */

#include "debug.h"
#include "hc_func.h"
#include "hc_marco.h"
#include "mem_pool_auto_enlarge.h"

#define HC_APPLY_UTILS_PCR_SNV_ERROR_QUAL (20)
#define HC_APPLY_UTILS_GET_POS(read)      (read->pos_start + 1)
#define HC_APPLY_UTILS_GET_MPOS(read)     (read->read.core.mpos + 1)
#define HC_APPLY_UTILS_GET_LSEQ(read)     (read->read_len)
#define HC_APPLY_UTILS_GET_ISIZE(read)    (read->insert_size)

static inline int hc_apply_utils_read_is_unmapped(p_hc_apply_one_read read);
inline int hc_apply_utils_mate_read_is_unmapped(p_hc_apply_one_read read);
static hts_pos_t hc_apply_utils_read_get_end(p_hc_apply_one_read read);
static int hc_apply_utils_get_adaptor_boundary(p_hc_apply_one_read read);
static int get_read_index_for_reference_coordinate(int alignmentStart, int refCoord, p_hc_apply_one_read read, int* out_cigar);

/**
 * @brief read 是否unmapped
 *
 * @param[in]   read    read
 *
 * @return int  1:true, 0:false
 */
static inline int hc_apply_utils_read_is_unmapped(p_hc_apply_one_read read)
{
    if (__glibc_unlikely((UNMAPPED(&read->read) || read->read.core.tid < 0 || read->pos_start))) {
        return 1;
    }
    return 0;
}

/**
 * @brief mate read 是否unmapped
 *
 * @param[in]   read    read
 *
 * @return int  1:true, 0:false
 */
inline int hc_apply_utils_mate_read_is_unmapped(p_hc_apply_one_read read)
{
    if (__glibc_unlikely(getBit(read->read.core.flag, WGS_SAM_FLAG_MATE_UNMAPPED) || read->read.core.mtid < 0 ||
                         HC_APPLY_UTILS_GET_MPOS(read) == 0)) {
        return 1;
    }
    return 0;
}

/**
 * Can the adaptor sequence of read be reliably removed from the read based on the alignment of
 * read and its mate?
 *
 * @param read the read to check
 * @return true if it can, false otherwise
 */
int hc_apply_utils_has_well_defined_fragment_size(p_hc_apply_one_read read)
{
    if (__glibc_unlikely(HC_APPLY_UTILS_GET_ISIZE(read) == 0)) {
        return 0;
    }
    if (__glibc_unlikely(!PAIRED(&read->read))) {
        return 0;
    }
    if (__glibc_unlikely(UNMAPPED(&read->read) || MATE_UNMAPPED(&read->read))) {
        return 0;
    }
    if (REVERSE_STRAND(&read->read) == MATE_REVERSE_STRAND(&read->read)) {
        return 0;
    }
    if (REVERSE_STRAND(&read->read)) {
        return bam_endpos(&read->read) > HC_APPLY_UTILS_GET_MPOS(read);
    }
    else {
        return read->pos_start <= HC_APPLY_UTILS_GET_MPOS(read) + HC_APPLY_UTILS_GET_ISIZE(read);
    }
}

/**
 * @brief HC Assemble 重建缓存read
 *
 * @param[in] read                  read
 **/
void hc_apply_utils_hard_clip_adaptor_sequence(p_hc_apply_one_read read)
{
    int adaptorBoundary = hc_apply_utils_get_adaptor_boundary(read);

    if (adaptorBoundary == INT_MIN || !(adaptorBoundary >= read->pos_start && adaptorBoundary <= hc_apply_utils_read_get_end(read))) {
        return;
    }
    if (REVERSE_STRAND(&read->read)) {
        hc_apply_utils_clip_by_reference_coordinates(-1, adaptorBoundary, read);
    }
    else {
        hc_apply_utils_clip_by_reference_coordinates(adaptorBoundary, -1, read);
    }
}

/**
 * @return 1-based closed-ended position, undefined if getContig() == null
 */
static hts_pos_t hc_apply_utils_read_get_end(p_hc_apply_one_read read)
{
    uint32_t* cigar = read->cigar;
    hts_pos_t rlen;
    for (uint32_t k = rlen = 0; k < read->cigar_len; ++k) {
        if (bam_cigar_type(bam_cigar_op(cigar[k])) & 2) rlen += bam_cigar_oplen(cigar[k]);
    }
    if (rlen == 0) rlen = 1;
    return read->pos_start + rlen - 1;  //-1是因为pos_start初始化时已经tranfer to 1base了。
}

/**
 * Finds the adaptor boundary around the read and returns the first base inside the adaptor that is closest to
 * the read boundary. If the read is in the positive strand, this is the first base after the end of the
 * fragment (Picard calls it 'insert'), if the read is in the negative strand, this is the first base before the
 * beginning of the fragment.
 *
 * There are two cases we need to treat here:
 *
 * 1) Our read is in the reverse strand :
 *
 *     <----------------------| *
 *   |--------------------->
 *
 *   in these cases, the adaptor boundary is at the mate start (minus one)
 *
 * 2) Our read is in the forward strand :
 *
 *   |---------------------->   *
 *     <----------------------|
 *
 *   in these cases the adaptor boundary is at the start of the read plus the inferred insert size (plus one)
 *
 * @param read the read being tested for the adaptor boundary
 * @return the reference coordinate for the adaptor boundary (effectively the first base IN the adaptor, closest to the read.
 * CANNOT_COMPUTE_ADAPTOR_BOUNDARY if the read is unmapped or the mate is mapped to another contig.
 */
static int hc_apply_utils_get_adaptor_boundary(p_hc_apply_one_read read)
{
    if (!hc_apply_utils_has_well_defined_fragment_size(read)) {
        return INT_MIN;
    }
    else if (REVERSE_STRAND(&read->read)) {
        return HC_APPLY_UTILS_GET_MPOS(read) - 1;
    }
    else {
        // insert_size换了数据类型，防止溢出，改labs
        int insertSize = labs(read->insert_size);  // the inferred insert size can be negative if the mate is mapped before the read (so we
                                                   // take the absolute value)
        return read->pos_start + insertSize;
    }
}

/**
 * Generic functionality to  clip a read, used internally by hardClipByReferenceCoordinatesLeftTail
 * and hardClipByReferenceCoordinatesRightTail. Should not be used directly.
 *
 * Note, it REQUIRES you to give the directionality of your hard clip (i.e. whether you're clipping the
 * left of right tail) by specifying either refStart < 0 or refStop < 0.
 *
 * @param refStart  first base to clip (inclusive)
 * @param refStop last base to clip (inclusive)
 * @param clippingOp clipping operation to be performed
 * @return a new read, without the clipped bases (May return empty, unclipped reads)
 */
void hc_apply_utils_clip_by_reference_coordinates(int refStart, int refStop, p_hc_apply_one_read read)
{
    int start;
    int stop;
    int res_cigar = WORKER_SIGAR_MAX;

    // Determine the read coordinate to start and stop hard clipping
    if (refStart < 0) {
        if (refStop < 0) {
            return;
        }
        start = 0;

        int stopPosAndOperator = get_read_index_for_reference_coordinate(hc_apply_utils_get_soft_start(read), refStop, read, &res_cigar);

        // if the refStop falls in a deletion, the above method returns the position after the deletion.  Since the stop we return here
        // is inclusive, we decrement the stop to avoid overclipping by one base.  As a result we do not clip the deletion, which is fine.
        stop = stopPosAndOperator - (consumes_read_bases(res_cigar) ? 0 : 1);
    }
    else {
        if (refStop >= 0) {
            return;
        }
        // unlike the above case where we clip the start fo the read, here we clip the end and returning the base to the right of a deletion
        // avoids overclipping
        start = get_read_index_for_reference_coordinate(hc_apply_utils_get_soft_start(read), refStart, read, &res_cigar);
        stop = HC_APPLY_UTILS_GET_LSEQ(read) - 1;
    }

    if (start == -1 || stop == -1) {
        return;
    }

    if (start < 0 || stop > (int)HC_APPLY_UTILS_GET_LSEQ(read) - 1) {
        return;
    }

    if (start > stop) {
        return;
    }

    if (start > 0 && stop < (int)HC_APPLY_UTILS_GET_LSEQ(read) - 1) {
        return;
    }
    if (__glibc_unlikely(stop > (int)HC_APPLY_UTILS_GET_LSEQ(read))) {
        stop = (int)HC_APPLY_UTILS_GET_LSEQ(read) - 1;
    }
    hc_assemble_utils_clip_read(read, start, stop, HARDCLIP);
}

/**
 * Find the 0-based index within a read base array corresponding to a given 1-based position in the reference, along with the cigar operator
 * of the element containing that base.  If the reference coordinate occurs within a deletion, the first index after the deletion is
 * returned. Note that this treats soft-clipped bases as if they align with the reference, which is useful for hard-clipping reads with soft
 * clips.
 *
 * @param alignmentStart        The soft start of the read on the reference
 * @param cigar                 The read's cigar
 * @param refCoord              The target reference coordinate
 * @return                      If the reference coordinate occurs before the read start or after the read end {@code
 * CLIPPING_GOAL_NOT_REACHED}; if the reference coordinate falls within an alignment block of the read's cigar, the corresponding read
 * coordinate; if the reference coordinate falls within a deletion, the first read coordinate after the deletion.  Note: if the last cigar
 * element is a deletion (which isn't meaningful), it returns {@code CLIPPING_GOAL_NOT_REACHED}.
 */
static int get_read_index_for_reference_coordinate(int alignmentStart, int refCoord, p_hc_apply_one_read read, int* out_cigar)
{
    uint32_t* cigar = read->cigar;
    int firstReadPosOfElement = 0;              // inclusive
    int firstRefPosOfElement = alignmentStart;  // inclusive
    int lastReadPosOfElement = 0;               // exclusive
    int lastRefPosOfElement = alignmentStart;   // exclusive
    uint32_t this_cigar, this_cigar_len;
    register uint32_t i;

    if (refCoord < alignmentStart) {
        return -1;
    }

    // advance forward through all the cigar elements until we bracket the reference coordinate
    for (i = 0; i < read->cigar_len; i++) {
        this_cigar = bam_cigar_op(cigar[i]);
        this_cigar_len = bam_cigar_oplen(cigar[i]);

        firstReadPosOfElement = lastReadPosOfElement;
        firstRefPosOfElement = lastRefPosOfElement;
        lastReadPosOfElement += consumes_read_bases(this_cigar) ? this_cigar_len : 0;
        lastRefPosOfElement += consumes_ref_bases(this_cigar) || this_cigar == BAM_CSOFT_CLIP ? this_cigar_len : 0;

        if (firstRefPosOfElement <= refCoord && refCoord < lastRefPosOfElement) {  // refCoord falls within this cigar element
            int readPosAtRefCoord = firstReadPosOfElement + (consumes_read_bases(this_cigar) ? (refCoord - firstRefPosOfElement) : 0);
            *out_cigar = this_cigar;
            return readPosAtRefCoord;
        }
    }
    return -1;
}

/**
 * Calculates the reference coordinate for the beginning of the read taking into account soft clips but not hard clips.
 *
 * Note: getUnclippedStart() adds soft and hard clips, this function only adds soft clips.
 *
 * @return the unclipped start of the read taking soft clips (but not hard clips) into account
 */
int hc_apply_utils_get_soft_start(p_hc_apply_one_read read)
{
    int softStart = read->pos_start;

    uint32_t* cigar = read->cigar;
    uint32_t this_cigar, this_cigar_len;
    register uint32_t i;

    // 已知位点
    for (i = 0; i < read->cigar_len; i++) {
        this_cigar = bam_cigar_op(cigar[i]);
        this_cigar_len = bam_cigar_oplen(cigar[i]);
        switch (this_cigar) {
            case WORKER_SIGAR_SOFT_CLIP: softStart -= this_cigar_len; break;
            case WORKER_SIGAR_HARD_CLIP: break;
            default: return softStart;
        }
    }
    return softStart;
}

/**
 * @brief 如果可能，通过调整碱基质量来修复同一片段的两个重叠读段
 *
 * @param read1     read
 * @param read2     read
 */
void hc_assemble_utils_adjust_overlapping_paired_qual(p_hc_apply_one_read first, p_hc_apply_one_read second)
{
    if (!first || !second) {
        statics_print_error("Unexpected null read");
        return;
    }
    if (strcmp(bam_get_qname(&first->read), bam_get_qname(&second->read))) {
        statics_print_error("Attempting to merge two reads with different names");
        exit(-1);
    }

    bool inOrder = hc_apply_utils_get_soft_start(first) < hc_apply_utils_get_soft_start(second);
    p_hc_apply_one_read firstRead = inOrder ? first : second;
    p_hc_apply_one_read secondRead = inOrder ? second : first;

    int out_cigar, i;
    hts_pos_t first_soft_start = hc_apply_utils_get_soft_start(firstRead);
    hts_pos_t second_soft_start = hc_apply_utils_get_soft_start(secondRead);
    hts_pos_t first_end = firstRead->pos_end;
    hts_pos_t second_end = secondRead->pos_end;
    if (first_end < secondRead->pos_start || firstRead->read.core.tid != secondRead->read.core.tid) {
        return;
    }
    int offset = get_read_index_for_reference_coordinate(first_soft_start, secondRead->pos_start, firstRead, &out_cigar);
    if (offset == -1 || (out_cigar & (BAM_CHARD_CLIP | BAM_CSOFT_CLIP))) {
        return;
    }
    int first_end_base = get_read_index_for_reference_coordinate(first_soft_start, first_end, firstRead, &out_cigar);
    int second_end_base = get_read_index_for_reference_coordinate(second_soft_start, second_end, secondRead, &out_cigar);

    int first_read_stop = offset;
    int second_offset = get_read_index_for_reference_coordinate(second_soft_start, secondRead->pos_start, secondRead, &out_cigar);
    int overlappling_count = assemble_min(first_end_base - first_read_stop, second_end_base - second_offset) + 1;

    uint8_t* firstReadBases = firstRead->seq;
    uint8_t* firstReadQuals = firstRead->qual;
    uint8_t* secondReadBases = secondRead->seq;
    uint8_t* secondReadQuals = secondRead->qual;

    for (i = 0; i < overlappling_count; i++) {
        int firstReadIndex = first_read_stop + i;
        int secondReadIndex = second_offset + i;

        char base1 = firstReadBases[firstReadIndex];
        char base2 = secondReadBases[secondReadIndex];

        if (base1 == base2) {
            firstReadQuals[firstReadIndex] = assemble_min(firstReadQuals[firstReadIndex], HC_APPLY_UTILS_PCR_SNV_ERROR_QUAL);
            secondReadQuals[secondReadIndex] = assemble_min(secondReadQuals[secondReadIndex], HC_APPLY_UTILS_PCR_SNV_ERROR_QUAL);
        }
        else {
            firstReadQuals[firstReadIndex] = 0;
            secondReadQuals[secondReadIndex] = 0;
        }
    }
    // 此处省略计算PCR indel qualities代码.
}

/**
 * @brief cigar获取ref长度
 *
 * @param[in] cigar         cigar
 * @param[in] cigar_len     cigar长度
 *
 * @return int ref长度
 */
int hc_assemble_utils_get_ref_length(uint32_t* cigar, uint32_t cigar_len)
{
    int ret_len = 0;
    uint32_t i = 0;
    uint32_t opr, opl;

    for (; i < cigar_len; i++) {
        opr = bam_cigar_op(cigar[i]);
        opl = bam_cigar_oplen(cigar[i]);

        if (consumes_ref_bases(opr)) {
            ret_len += opl;
        }
    }
    return ret_len;
}

/**
 * @brief 初始化apply region
 *
 * @param apply             apply
 * @param region            region
 * @param chr_ref           染色体参考Base
 * @param chr_ref_length    染色体参考Base长度
 */
void hc_apply_utils_init_region(p_hc_apply apply, p_hc_region_active_storage region, const uint8_t* chr_ref, uint32_t chr_ref_length)
{
    apply->region = region;
    apply->ref = chr_ref;
    apply->ref_chr_len = chr_ref_length;
}