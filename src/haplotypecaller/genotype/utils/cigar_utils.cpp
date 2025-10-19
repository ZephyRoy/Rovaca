#include "cigar_utils.h"

#include "cigar_builder.h"
#include "rovaca_logger.h"
#include "genotype_struct.h"

namespace rovaca
{

pCigar CigarUtils::clip_cigar(const uint32_t *cigar, uint32_t cigar_num, int64_t start, int64_t stop, uint32_t clip_op, pMemoryPool pool)
{
    CHECK_CONDITION_EXIT(!cigar_op_is_soft_clip(clip_op) && !cigar_op_is_hard_clip(clip_op), "not a clipping operator");

    pCigarBuilder cigar_builder = CigarBuilder::create(pool);

    bool clip_left = (start == 0);
    int64_t element_start = 0, element_end, unclipped_length, clipped_length;
    uint32_t idx, cigar_element, current_op, current_op_len;

    for (idx = 0; idx < cigar_num; ++idx) {
        cigar_element = cigar[idx];
        current_op = bam_cigar_op(cigar_element);
        current_op_len = bam_cigar_oplen(cigar_element);

        // copy hard clips
        if (cigar_op_is_hard_clip(current_op)) {
            cigar_builder->add(current_op_len << BAM_CIGAR_SHIFT | current_op);
            continue;
        }
        element_end = element_start + (consumes_read_bases(current_op) ? current_op_len : 0);

        // element precedes start or follows end of clip, copy it to new cigar
        if (element_end <= start || element_start >= stop) {
            // edge case: deletions at edge of clipping are meaningless and we skip them
            if (consumes_read_bases(current_op) || (element_start != start && element_start != stop)) {
                cigar_builder->add(current_op_len << BAM_CIGAR_SHIFT | current_op);
            }
        }
        else {  // otherwise, some or all of the element is soft-clipped
            unclipped_length = clip_left ? element_end - stop : start - element_start;
            clipped_length = current_op_len - unclipped_length;
            if (unclipped_length <= 0) {  // totally clipped
                if (consumes_read_bases(current_op)) {
                    cigar_builder->add(current_op_len << BAM_CIGAR_SHIFT | clip_op);
                }
            }
            else if (clip_left) {
                cigar_builder->add(clipped_length << BAM_CIGAR_SHIFT | clip_op);
                cigar_builder->add(unclipped_length << BAM_CIGAR_SHIFT | current_op);
            }
            else {
                cigar_builder->add(unclipped_length << BAM_CIGAR_SHIFT | current_op);
                cigar_builder->add(clipped_length << BAM_CIGAR_SHIFT | clip_op);
            }
        }
        element_start = element_end;
    }

    return cigar_builder->make();
}

int64_t CigarUtils::alignment_start_shift(const uint32_t *cigar, uint32_t cigar_num, int64_t num_clipped)
{
    int64_t ref_bases_clipped = 0, element_start = 0, element_end;
    uint32_t i, ce, op, op_len;
    for (i = 0; i < cigar_num; ++i) {
        ce = cigar[i];
        op = bam_cigar_op(ce);
        op_len = bam_cigar_oplen(ce);
        if (cigar_op_is_hard_clip(op)) {
            continue;
        }
        element_end = element_start + (consumes_read_bases(op) ? op_len : 0);

        if (element_end <= num_clipped) {
            // totally within clipped span -- this includes deletions immediately following clipping
            ref_bases_clipped += consumes_ref_bases(op) ? op_len : 0;
        }
        else if (element_start < num_clipped) {
            // clip in middle of element, which means the element necessarily consumes read bases
            int64_t clipped_length = num_clipped - element_start;
            ref_bases_clipped += consumes_ref_bases(op) ? clipped_length : 0;
            break;
        }
        element_start = element_end;
    }

    return ref_bases_clipped;
}

pCigar CigarUtils::str2uint(const char *cigar, uint32_t cigar_num, pMemoryPool pool)
{
    auto *new_ciagr = new ALLOC_FLEXIBLE_IN_POOL(pool, Cigar, cigar_num, uint32_t) Cigar{cigar_num};
    uint32_t *buffer = new_ciagr->data;
    auto num = (size_t)cigar_num;

    if (sam_parse_cigar(cigar, nullptr, &buffer, &num) != cigar_num) {
        RovacaLogger::error("sam_parse_cigar error");
        return nullptr;
    }
    return new_ciagr;
}

std::string CigarUtils::uint2str(const uint32_t *cigar, uint32_t cigar_num)
{
    std::string ret;
    for (uint32_t i = 0, op, op_pen; i < cigar_num; ++i) {
        op = bam_cigar_op(cigar[i]);
        op_pen = bam_cigar_oplen(cigar[i]);
        ret.append(std::to_string(op_pen));
        ret.append(1, BAM_CIGAR_STR[op]);
    }
    return ret;
}

uint32_t CigarUtils::count_ref_bases_and_clips(const uint32_t *cigar, uint32_t cigar_num, uint32_t start_idx, uint32_t end_idx)
{
    return count_ref_bases_and_maybe_also_clips(cigar, cigar_num, start_idx, end_idx, true, true);
}

uint32_t CigarUtils::count_ref_bases_and_soft_clips(const uint32_t *cigar, uint32_t cigar_num, uint32_t start_idx, uint32_t end_idx)
{
    return count_ref_bases_and_maybe_also_clips(cigar, cigar_num, start_idx, end_idx, true, false);
}

uint32_t CigarUtils::count_ref_bases_and_maybe_also_clips(const uint32_t *cigar, uint32_t cigar_num, uint32_t start_idx, uint32_t end_idx,
                                                          bool include_soft_clips, bool include_hard_clips)
{
    CHECK_CONDITION_EXIT(end_idx > cigar_num || end_idx < start_idx, "invalid index");
    uint32_t op, result = 0;
    for (uint32_t i = start_idx; i < end_idx; ++i) {
        op = bam_cigar_op(cigar[i]);
        if (consumes_ref_bases(op) || (include_soft_clips && cigar_op_is_soft_clip(op)) || (include_hard_clips && cigar_op_is_hard_clip(op))) {
            result += bam_cigar_oplen(cigar[i]);
        }
    }
    return result;
}

}  // namespace rovaca