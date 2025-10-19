#include "clipping_op.h"

#include "genotype_macors.h"
#include "genotype_struct.h"
#include "rovaca_logger.h"
#include "read_record.h"
#include "utils/cigar_utils.h"
#include "utils/read_record_utils.h"

namespace rovaca
{

pReadRecord ClippingOp::apply(pReadRecord const& original, ClippingRepresentation op)
{
    switch (op) {
        case HARDCLIP_BASES: {
            return apply_hard_clip_bases(original, _start, _stop);
        }
        case WRITE_NS:
        case WRITE_Q0S:
        case WRITE_NS_Q0S:
        case SOFTCLIP_BASES:
        case REVERT_SOFTCLIPPED_BASES:
        default: {
            RovacaLogger::error("unexpected Clipping operator type: {}", (uint32_t)op);
            exit(EXIT_FAILURE);
        }
    }
}

pReadRecord ClippingOp::apply_hard_clip_bases(pReadRecord const& read, int64_t start, int64_t stop) const
{
    int64_t new_length = read->seq_length() - (stop - start + 1);
    if (new_length == 0) {
        return nullptr;
    }

    uint32_t* old_cigar = read->cigar();
    uint32_t old_cigar_num = read->cigar_length();

    pCigar cigar = read->is_unmapped() ? new ALLOC_FLEXIBLE_IN_POOL(_pool, Cigar, 0, uint32_t) Cigar{0}
                                       : CigarUtils::clip_cigar(old_cigar, old_cigar_num, start, stop + 1, BAM_CHARD_CLIP, _pool);

    auto* new_bases = new ALLOC_MEM_IN_POOL(_pool, sizeof(uint8_t) * new_length) uint8_t[new_length]{};
    auto* new_quals = new ALLOC_MEM_IN_POOL(_pool, sizeof(uint8_t) * new_length) uint8_t[new_length]{};
    int64_t copy_start = (start == 0) ? stop + 1 : 0;

    // @note: read的原始数据可能源于bam1_t，也可能源于p_hc_apply_one_read，由于编码方式差异，无法直接获取char*指针，只能一个个赋值
    for (int64_t idx = 0; idx < new_length; ++idx) {
        new_bases[idx] = read->seq_i(idx + copy_start);
    }
    memcpy(new_quals, read->qual() + copy_start, new_length);

    pReadRecord hard_clipped_read = ReadRecordUtils::copy_read(read, new_bases, new_quals, new_length, cigar, _bam_pool, _pool);
    if (start == 0 && !read->is_unmapped()) {
        hard_clipped_read->set_start(read->get_start() + CigarUtils::alignment_start_shift(old_cigar, old_cigar_num, stop + 1));
    }

    // if (ReadUtils.hasBaseIndelQualities(read)) {
    //     final byte[] newBaseInsertionQuals = new byte[newLength];
    //     final byte[] newBaseDeletionQuals = new byte[newLength];
    //     System.arraycopy(ReadUtils.getBaseInsertionQualities(read), copyStart, newBaseInsertionQuals, 0, newLength);
    //     System.arraycopy(ReadUtils.getBaseDeletionQualities(read), copyStart, newBaseDeletionQuals, 0, newLength);
    //     ReadUtils.setInsertionBaseQualities(hardClippedRead, newBaseInsertionQuals);
    //     ReadUtils.setDeletionBaseQualities(hardClippedRead, newBaseDeletionQuals);
    // }

    return hard_clipped_read;
}

}  // namespace rovaca