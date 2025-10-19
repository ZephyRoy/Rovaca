#include "hc_assemble_cigar_utils.h"

#include "hc_apply.h"

/**
 * @brief get end for user defined read. 1base.
 *
 * @param read
 * @return hts_pos_t
 */
hts_pos_t hc_assemble_utils_read_get_end(p_hc_apply_one_read read)
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
 * @brief Merge consecutive identical cigars.
 * @note Only for debug, may incur certain expenses.
 * @param read
 */
void merge_consecutive_identical_cigar(p_hc_apply_one_read read)
{
    uint32_t* cigar = read->cigar;
    uint32_t next_oplen, op_len;
    uint32_t next_op, this_op;
    uint32_t i, j;
    for (i = 0, j = 1; j < read->cigar_len; ++j) {
        this_op = bam_cigar_op(cigar[i]);
        next_op = bam_cigar_op(cigar[j]);
        if (next_op == this_op) {
            op_len = bam_cigar_oplen(cigar[i]);
            next_oplen = bam_cigar_oplen(cigar[j]);
            uint32_t new_op = (this_op & BAM_CIGAR_MASK) | ((op_len + next_oplen) << BAM_CIGAR_SHIFT);
            cigar[i] = new_op;
        }
        else {
            i++;
            if (i != j) {
                cigar[i] = cigar[j];
            }
        }
    }
    read->cigar_len = i + 1;
    read->pos_end = hc_assemble_utils_read_get_end(read);
    return;
}
