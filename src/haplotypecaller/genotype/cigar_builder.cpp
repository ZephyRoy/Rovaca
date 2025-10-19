#include "cigar_builder.h"

#include "genotype_struct.h"
#include "htslib/sam.h"
#include "rovaca_logger.h"
#include "utils/cigar_utils.h"

namespace rovaca
{

pCigarBuilder CigarBuilder::create(pMemoryPool pool) { return new ALLOC_TYPE_IN_POOL(pool, CigarBuilder) CigarBuilder{pool}; }

pCigarBuilder CigarBuilder::create(bool remove_deletions_at_ends, pMemoryPool pool)
{
    return new ALLOC_TYPE_IN_POOL(pool, CigarBuilder) CigarBuilder{remove_deletions_at_ends, pool};
}

CigarBuilder* CigarBuilder::add(uint32_t cigar_element)
{
    uint32_t current_op = bam_cigar_op(cigar_element);
    uint32_t current_op_len = bam_cigar_oplen(cigar_element);
    if (current_op_len == 0) {
        return this;
    }

    // skip a deletion after clipping ie at the beginning of the read
    // note the edge case of a deletion following a leading insertion, which we also skip
    if (_remove_deletions_at_ends && cigar_op_is_del(current_op)) {
        if (skip1() || skip2()) {
            _leading_deletion_bases_removed += current_op_len;
            return this;
        }
    }

    advance_section_and_validate_cigar_order(current_op);

    // merge consecutive elements with the same operator
    if (current_op == _last_operator) {
        uint32_t loop_op_len = bam_cigar_oplen(_cigar_elements.data[_cigar_elements.num - 1]);
        _cigar_elements.data[_cigar_elements.num - 1] = ((loop_op_len + current_op_len) << BAM_CIGAR_SHIFT) | current_op;
    }
    else {
        if (cigar_op_is_uninit(_last_operator)) {
            _cigar_elements.data[_cigar_elements.num++] = cigar_element;
            _last_operator = current_op;
        }
        else if (cigar_op_is_clipping(current_op)) {
            // if we have just started clipping on the right and realize the last operator was a deletion, remove it
            // if we have just started clipping on the right and the last two operators were a deletion and insertion, remove the deletion
            if (_remove_deletions_at_ends && !consumes_read_bases(_last_operator) && !cigar_op_is_clipping(_last_operator)) {
                _trailing_deletion_bases_removed += bam_cigar_oplen(_cigar_elements.data[_cigar_elements.num - 1]);
                _cigar_elements.data[_cigar_elements.num - 1] = cigar_element;
                _last_operator = current_op;
            }
            else if (_remove_deletions_at_ends && last_two_elements_were_deletion_and_insertion()) {
                _trailing_deletion_bases_removed += bam_cigar_oplen(_cigar_elements.data[_cigar_elements.num - 2]);
                _cigar_elements.data[_cigar_elements.num - 2] = _cigar_elements.data[_cigar_elements.num - 1];
                _cigar_elements.data[_cigar_elements.num - 1] = cigar_element;
            }
            else {
                _cigar_elements.data[_cigar_elements.num++] = cigar_element;
                _last_operator = current_op;
            }
        }
        else if (cigar_op_is_del(current_op) && cigar_op_is_ins(_last_operator)) {
            // The order of deletion and insertion elements is arbitrary, so to standardize we shift deletions to the left that is, we place
            // the deletion before the insertion and shift the insertion right if the element before the insertion is another deletion, we
            // merge in the new deletion note that the last operator remains an insertion
            uint32_t cigar_num = _cigar_elements.num;
            if (cigar_num > 1 && cigar_op_is_del(bam_cigar_op(_cigar_elements.data[cigar_num - 2]))) {
                uint32_t sub2_len = bam_cigar_oplen(_cigar_elements.data[cigar_num - 2]);
                _cigar_elements.data[cigar_num - 2] = ((sub2_len + current_op_len) << BAM_CIGAR_SHIFT) | BAM_CDEL;
            }
            else {
                _cigar_elements.num++;
                _cigar_elements.data[cigar_num] = _cigar_elements.data[cigar_num - 1];
                _cigar_elements.data[cigar_num - 1] = cigar_element;
            }
        }
        else {
            _cigar_elements.data[_cigar_elements.num] = cigar_element;
            _cigar_elements.num++;
            _last_operator = current_op;
        }
    }
    return this;
}

CigarBuilder* CigarBuilder::add_all(uint32_t* cigar, uint32_t cigar_num)
{
    CHECK_CONDITION_EXIT(cigar == nullptr, "cigar pointer cannot be nullptr");
    uint32_t i;
    for (i = 0; i < cigar_num; ++i) {
        add(cigar[i]);
    }
    return this;
}

pCigar CigarBuilder::make(bool allow_empty)
{
    if (is_left_soft_clip() && cigar_op_is_soft_clip(bam_cigar_op(_cigar_elements.data[0]))) {
        RovacaLogger::error("cigar is completely soft-clipped");
        return nullptr;
    }
    _trailing_deletion_bases_removed_in_make = 0;
    if (_remove_deletions_at_ends && cigar_op_is_del((_last_operator))) {
        _trailing_deletion_bases_removed_in_make = bam_cigar_oplen(_cigar_elements.data[_cigar_elements.num - 1]);
        _cigar_elements.num--;
    }
    else if (_remove_deletions_at_ends && last_two_elements_were_deletion_and_insertion()) {
        _trailing_deletion_bases_removed_in_make = bam_cigar_oplen(_cigar_elements.data[_cigar_elements.num - 2]);
        _cigar_elements.data[_cigar_elements.num - 2] = _cigar_elements.data[_cigar_elements.num - 1];
        _cigar_elements.num--;
    }

    CHECK_CONDITION_EXIT(!(allow_empty || _cigar_elements.num != 0),
                         "no cigar elements left after removing leading and trailing deletions");

    pCigar cigar = new ALLOC_FLEXIBLE_IN_POOL(_pool, Cigar, _cigar_elements.num, uint32_t) Cigar{_cigar_elements.num};
    memcpy(cigar->data, _cigar_elements.data, _cigar_elements.num * sizeof(uint32_t));
    return cigar;
}

CigarBuilder::Result CigarBuilder::make_and_record_deletions_removed_result()
{
    CigarBuilder::Result r{};
    r.cigar = make();
    r.leading_deletion_bases_removed = leading_deletion_bases_removed();
    r.trailing_deletion_bases_removed = trailing_deletion_bases_removed();
    return r;
}

CigarBuilder::CigarBuilder(bool remove_deletions_at_ends, pMemoryPool gpool)
    : _pool(gpool)
    , _cigar_elements({0, {}})
    , _remove_deletions_at_ends(remove_deletions_at_ends)
{}

bool CigarBuilder::last_two_elements_were_deletion_and_insertion()
{
    return cigar_op_is_ins(_last_operator) && _cigar_elements.num > 1 &&
           cigar_op_is_del(bam_cigar_op(_cigar_elements.data[_cigar_elements.num - 2]));
}

void CigarBuilder::advance_section_and_validate_cigar_order(uint32_t op)
{
    if (cigar_op_is_hard_clip(op)) {
        if (is_left_soft_clip() || is_middle() || is_right_soft_clip()) {
            _section = Section::RIGHT_HARD_CLIP;
        }
    }
    else if (cigar_op_is_soft_clip(op)) {
        CHECK_CONDITION_EXIT(is_right_hard_clip(), "cigar has already reached its right hard clip");
        if (is_left_hard_clip()) {
            _section = Section::LEFT_SOFT_CLIP;
        }
        else if (is_middle()) {
            _section = Section::RIGHT_SOFT_CLIP;
        }
    }
    else {
        CHECK_CONDITION_EXIT(is_right_soft_clip() || is_right_hard_clip(), "cigar has already reached right clip");
        if (is_left_hard_clip() || is_left_soft_clip()) {
            _section = Section::MIDDLE;
        }
    }
}

bool CigarBuilder::skip1() const { return cigar_op_is_uninit(_last_operator) || cigar_op_is_clipping(_last_operator); }

bool CigarBuilder::skip2() const
{
    return cigar_op_is_ins(_last_operator) &&
           (_cigar_elements.num == 1 || cigar_op_is_clipping(bam_cigar_op(_cigar_elements.data[_cigar_elements.num - 2])));
}

}  // namespace rovaca