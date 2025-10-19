#include "pileup_element.h"

#include <memory_resource>

#include "read_record.h"
#include "utils/cigar_utils.h"

namespace rovaca
{

// static constexpr uint8_t DELETION_QUAL = (char)16;
// static constexpr char DELETION_BASE = (char)'D';

pPileupElement PileupElement::create(pReadRecord read, uint32_t offset, uint32_t cigar, uint32_t cigar_offset, uint32_t offset_in_cigar,
                                     pMemoryPool pool)
{
    return new ALLOC_TYPE_IN_POOL(pool, PileupElement) PileupElement{read, offset, cigar, cigar_offset, offset_in_cigar};
}

PileupElement::PileupElement(pReadRecord read, uint32_t offset, uint32_t cigar, uint32_t cigar_offset, uint32_t offset_in_cigar)
    : _read(read)
    , _offset(offset)
    , _current_cigar_offset(cigar_offset)
    , _offset_in_current_cigar(offset_in_cigar)
    , _current_cigar_element(cigar)
{
    // _qual = is_deletion() ? DELETION_QUAL : _read->qual_i(_offset);
    // _base = is_deletion() ? DELETION_BASE : _read->seq_i(_offset);
    _is_at_end_of_current_cigar = _offset_in_current_cigar == bam_cigar_oplen(_current_cigar_element) - 1;
}

bool PileupElement::is_deletion() const { return cigar_op_is_del(bam_cigar_op(_current_cigar_element)); }

pReadRecord PileupElement::get_read() { return _read; }

uint32_t PileupElement::get_next_on_genome_cigar_element() { return get_nearest_on_genome_cigar_element(Direction::next); }

uint32_t PileupElement::get_nearest_on_genome_cigar_element(Direction d)
{
    // chatGPT优化后的代码，虽然有点小问题，但我很认可他的思路
    uint32_t n_cigar_elements = _read->cigar_length();
    uint32_t i = _current_cigar_offset + d;
    for (; i < n_cigar_elements; i += d) {
        if (consumes_ref_bases(bam_cigar_op(_read->cigar_i(i)))) {
            return _read->cigar_i(i);
        }
    }
    // getting here means that you didn't find anything
    return BAM_CUNINITIALIZE;
}

bool PileupElement::is_before_deletion_start()
{
    return !is_deletion() && at_end_of_current_cigar() && has_operator(get_next_on_genome_cigar_element(), BAM_CDEL);
}

bool PileupElement::at_end_of_current_cigar() const { return _is_at_end_of_current_cigar; }

bool PileupElement::has_operator(uint32_t maybe_cigar_element, uint32_t to_match_operator)
{
    return 0 != bam_cigar_oplen(maybe_cigar_element) && bam_cigar_op(maybe_cigar_element) == to_match_operator;
}

bool PileupElement::is_before_insertion() { return is_immediately_before(BAM_CINS); }

bool PileupElement::is_immediately_before(uint32_t cigar_operator)
{
    return at_end_of_current_cigar() && cigar_operator == get_adjacent_operator(Direction::next);
}

uint32_t PileupElement::get_adjacent_operator(Direction d)
{
    int32_t i = _current_cigar_offset + d;
    if (i < 0 || i >= (int32_t)_read->cigar_length()) {
        return BAM_CUNINITIALIZE;
    }
    return bam_cigar_op(_read->cigar_i(i));
}

}  // namespace rovaca