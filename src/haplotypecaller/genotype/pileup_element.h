#ifndef ROVACA_HC_PLIEUP_ELEMENT_H_
#define ROVACA_HC_PLIEUP_ELEMENT_H_

#include <cstdint>
#include <set>

#include "forward.h"
#include "genotype_macors.h"

namespace rovaca
{

class PileupElement
{
    /* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
private:
    pReadRecord _read;
    enum Direction { next = 1, prev = -1 };

    /* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
public:
    uint32_t _offset : 11;
    uint32_t _current_cigar_offset : 10;
    uint32_t _offset_in_current_cigar : 10;
    uint32_t _is_at_end_of_current_cigar : 1;
    uint32_t _current_cigar_element;
    // uint8_t _qual;
    // uint8_t _base;

    /* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
public:
    DISALLOW_COPY_AND_ASSIGN(PileupElement);

    static pPileupElement create(pReadRecord read, uint32_t offset, uint32_t cigar, uint32_t cigar_offset, uint32_t offset_in_cigar,
                                 pMemoryPool pool);

    /**
     * @brief Is this element a deletion w.r.t. the reference genome?
     * @return true if this is a deletion, false otherwise
     */
    bool is_deletion() const;

    /**
     * Get the read for this pileup element. Returns the live object stored in this pileup element.
     * @return a non-null Read
     */
    pReadRecord get_read();

#if 0
    uint8_t get_qual() const;
    uint8_t get_base() const;
#endif

    /**
     * Is the current element immediately before a deletion, but itself not a deletion? Suppose we are aligning a read
     * with cigar 3M2D1M.  This function is true if we are in the last cigar position of the 3M, but not if we are in
     * the 2D itself.
     * @return true if the next alignment position is a deletion w.r.t. the reference genome
     */
    bool is_before_deletion_start();

    /**
     * Is the current position at the end of the current cigar? For example, if we are in element 3M, this function
     * returns true if we are at offsetInCurrentCigar of 2, but not 0 or 1.
     * @return true if we're at the end of the current cigar
     */
    bool at_end_of_current_cigar() const;

    /**
     * Get the cigar element of the next/previous genomic aligned position
     * @see #getPreviousOnGenomeCigarElement() for more details
     * @return a CigarElement, or null (indicating that no next element exists)
     */
    uint32_t get_next_on_genome_cigar_element();

    /**
     * Does an insertion occur immediately after the current position on the genome?
     * @return true an insertion occurs immediately after the current position on the genome, false
     * otherwise
     */
    bool is_before_insertion();

    /* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
private:
    /**
     * @brief Construct a new Pileup Element object

     * @param read a non-null read to pileup
     * @param base_offset the offset uint32_to the read's base / qual vector aligned to this position on the genome. If
     * the current cigar element is a deletion, offset should be the offset of the
     * last M/=/X position.
     * @param ce a non-null CigarElement that indicates the cigar element aligning the read to the genome
     * @param ce_index the offset of currentElement in read.getCigar().getElement(currentCigarOffset) == currentElement)
     * @param offset_in_ce how far uint32_to the currentElement are we in our alignment to the genome?
     */
    PileupElement(pReadRecord read, uint32_t base_offset, uint32_t ce, uint32_t ce_index, uint32_t offset_in_ce);

    static bool has_operator(uint32_t p_maybe_cigar_element, uint32_t to_match_operator);

    /**
     * @return true if a given operator event occurs immediately after the current position on the
     * genome, false otherwise.
     */
    bool is_immediately_before(uint32_t cigar_operator);

    /**
     * Helper function to get cigar operator right next to this position.
     * @param direction PREVIOUS if we want before, NEXT if we want after
     * @return the next Cigar operator in the given direction, or null if there is none
     * TODO: 失败时返回 999
     */
    uint32_t get_adjacent_operator(Direction d);

    /**
     * Helper function to get the cigar element of the next or previous genomic position
     * @param direction the direction to look in
     * @return nearest on-genome CigarElement or null if no such element exists
     */
    uint32_t get_nearest_on_genome_cigar_element(Direction d);

};  // PileupElement end

}  // namespace rovaca

#endif