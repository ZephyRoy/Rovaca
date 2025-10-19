#ifndef READS_FILTER_LIB_H
#define READS_FILTER_LIB_H

#include "htslib/sam.h"

class ReadFilterLib
{
public:
    static bool valid_alignment_start(bam1_t* read);
    static bool valid_alignment_end(bam1_t* read);
    static bool alignment_agree_with_hdr(bam1_t* read, bam_hdr_t* header);
    static bool has_read_group(bam1_t* read);
    static bool matching_base_and_qual(bam1_t* read);
    static bool read_len_equal_cigar_len(bam1_t* read);
    static bool seq_is_stored(bam1_t* read);
    static bool is_well_formed(bam1_t* read, bam_hdr_t* header);
    static bool is_good_cigar(bam1_t* read);
    static bool is_first_last_cigar_deletion(bam1_t* read);
};
/**
 * @brief 判断read是否占用reference.即是否有M、D、N、=和X.
 * @param cigar
 * @return int
 */
static int cigar_consumes_reference(uint32_t cigar) { return bam_cigar_type(cigar) & 0x2; }

/**
 * @brief 过滤非法首尾，是否有RG，碱基与质量匹配否，染色体信息与header冲突否。
 * @param read
 * @return true
 * @return false
 */
bool ReadFilterLib::is_well_formed(bam1_t* read, bam_hdr_t* header)
{
    return valid_alignment_start(read) && valid_alignment_end(read) && has_read_group(read) && read_len_equal_cigar_len(read) &&
           matching_base_and_qual(read) && seq_is_stored(read) && alignment_agree_with_hdr(read, header);
}

/**
 * @brief 判断是否有连续的Indel，有则过滤。
 *        判断是否有N，有则过滤。
 *        判断read是否比对上reference.即是否有M、D、N、=和X.没有则过滤。
 * @param read
 * @return true
 * @return false
 */
bool ReadFilterLib::is_good_cigar(bam1_t* read)
{
    uint32_t* cigar = bam_get_cigar(read);
    int n_cigar = read->core.n_cigar;
    bool ref_consumes = false;
    bool prev_indel = false;
    bool is_indel = false;
    bool valid_deletion = is_first_last_cigar_deletion(read);
    if (!valid_deletion) {
        return false;
    }
    for (int i = 0; i < n_cigar; i++) {
        int op = bam_cigar_op(cigar[i]);
        // consective insertion.
        is_indel = ((op & (BAM_CINS | BAM_CDEL)) > 0);
        if (is_indel && prev_indel) {
            return false;
        }
        // whether the list has any N operators.
        if (op == BAM_CREF_SKIP) {
            return false;
        }
        // consumes reference && oplen > 0, pass when at least one op meet.
        if (cigar_consumes_reference(op)) {
            ref_consumes = true;
        }
        prev_indel = is_indel;
    }
    return ref_consumes;
}

/**
 * @brief 判断比对起始位置是否合法.
 * @param read
 * @return true
 * @return false
 */
bool ReadFilterLib::valid_alignment_start(bam1_t* read) { return (read->core.pos >= 0); }

/**
 * @brief 判断比对终止位置是否合法.
 * @param read
 * @return true
 * @return false
 */
bool ReadFilterLib::valid_alignment_end(bam1_t* read) { return (bam_endpos(read) >= read->core.pos); }

/**
 * @brief 判断read的染色体id正常并且起始点没超出的染色体长度范围.
 * @param read
 * @param header
 * @return true
 * @return false
 */
bool ReadFilterLib::alignment_agree_with_hdr(bam1_t* read, bam_hdr_t* header)
{
    return (read->core.tid >= 0) && (read->core.tid < header->n_targets) && (read->core.pos <= *(header->target_len + read->core.tid));
}

/**
 * @brief 判断是否有read group信息.
 * @param read
 * @return true
 * @return false
 */
bool ReadFilterLib::has_read_group(bam1_t* read)
{
    uint8_t* rg = bam_aux_get(read, "RG");
    return (rg != NULL);
}

/**
 * @brief 判断碱基和质量得分是否匹配.
 * @param read
 * @return true
 * @return false
 */
bool ReadFilterLib::matching_base_and_qual(bam1_t* read) { return bam_get_qual(read) != nullptr; }

/**
 * @brief 判断read长度是否等于cigar长度.
 * @param read
 * @return true
 * @return false
 */
bool ReadFilterLib::read_len_equal_cigar_len(bam1_t* read)
{
    return (read->core.l_qseq == bam_cigar2qlen(read->core.n_cigar, bam_get_cigar(read)));
}

/**
 * @brief 判断是否存储了序列信息.
 * @param read
 * @return true
 * @return false
 */
bool ReadFilterLib::seq_is_stored(bam1_t* read) { return (read->core.l_qseq > 0); }

/**
 * @brief 判断cigar首尾是否是deletion，不考虑clip.
 * @param read
 * @return true, 正常
 * @return false，需要过滤
 */
bool ReadFilterLib::is_first_last_cigar_deletion(bam1_t* read)
{
    uint32_t* cigar = bam_get_cigar(read);
    int n_cigar = read->core.n_cigar;
    uint32_t* p1 = cigar;
    uint32_t* p2 = cigar + n_cigar - 1;
    p1 += ((bam_cigar_op(*p1) & (BAM_CSOFT_CLIP | BAM_CHARD_CLIP)) > 0);
    p2 -= ((bam_cigar_op(*p2) & (BAM_CSOFT_CLIP | BAM_CHARD_CLIP)) > 0);
    return !((bam_cigar_op(*p1) | bam_cigar_op(*p2)) & BAM_CDEL);
}
#endif  // READS_FILTER_LIB_H