#ifndef ROVACA_HC_CIGAR_UTILS_H_
#define ROVACA_HC_CIGAR_UTILS_H_
#include <string>

#include "forward.h"
#include "genotype_macors.h"
#include "htslib/sam.h"

namespace rovaca
{

namespace CigarUtils
{

/*!
 * @brief 判断
 */
#define consumes_read_bases(op) (bam_cigar_type(op) & 0x1)
#define consumes_ref_bases(op)  (bam_cigar_type(op) & 0x2)
#define consumes_all_bases(op)  (bam_cigar_type(op) & 0x3)

// MIDNSHP=XB
#define cigar_op_is_match(op)     (BAM_CMATCH == op)
#define cigar_op_is_ins(op)       (BAM_CINS == op)
#define cigar_op_is_del(op)       (BAM_CDEL == op)
#define cigar_op_is_ref_skip(op)  (BAM_CREF_SKIP == op)
#define cigar_op_is_soft_clip(op) (BAM_CSOFT_CLIP == op)
#define cigar_op_is_hard_clip(op) (BAM_CHARD_CLIP == op)
#define cigar_op_is_pad(op)       (BAM_CPAD == op)
#define cigar_op_is_equal(op)     (BAM_CEQUAL == op)
#define cigar_op_is_diff(op)      (BAM_CDIFF == op)
#define cigar_op_is_back(op)      (BAM_CBACK == op)

#define cigar_op_is_uninit(op)    (BAM_CUNINITIALIZE == op)
#define cigar_op_is_indel(op)     (cigar_op_is_ins(op) || cigar_op_is_del(op))
#define cigar_op_is_clipping(op)  (cigar_op_is_hard_clip(op) || cigar_op_is_soft_clip(op))
#define cigar_op_is_alignment(op) (cigar_op_is_match(op) || cigar_op_is_equal(op) || cigar_op_is_diff(op))

/*! @brief Given a cigar string, soft clip up to leftClipEnd and soft clip starting at rightClipBegin */
pCigar clip_cigar(const uint32_t* cigar, uint32_t cigar_num, int64_t start, int64_t stop, uint32_t clip_op, pMemoryPool pool);

/*! @brief How many bases to the right does a read's alignment start shift given its cigar and the number of left soft clips */
int64_t alignment_start_shift(const uint32_t* cigar, uint32_t cigar_num, int64_t num_clipped);

/*! @brief debug使用，字符串cigar转bam1_t类型cigar */
pCigar str2uint(const char* cigar, uint32_t cigar_num, pMemoryPool pool);
std::string uint2str(const uint32_t* cigar, uint32_t cigar_num);

uint32_t count_ref_bases_and_clips(const uint32_t* cigar, uint32_t cigar_num, uint32_t start_idx, uint32_t end_idx);
uint32_t count_ref_bases_and_soft_clips(const uint32_t* cigar, uint32_t cigar_num, uint32_t start_idx, uint32_t end_idx);
uint32_t count_ref_bases_and_maybe_also_clips(const uint32_t* cigar, uint32_t cigar_num, uint32_t start_idx, uint32_t end_idx,
                                              bool include_soft_clips, bool include_hard_clips);

}  // namespace CigarUtils

}  // namespace rovaca

#endif  // ROVACA_HC_CIGAR_UTILS_H_
