#ifndef ROVACA_HC_CIGAR_BUILDER_H_
#define ROVACA_HC_CIGAR_BUILDER_H_
#include "forward.h"
#include "genotype_macors.h"

namespace rovaca
{

class CigarBuilder
{
    /* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
private:
    enum Section { LEFT_HARD_CLIP = 0, LEFT_SOFT_CLIP, MIDDLE, RIGHT_SOFT_CLIP, RIGHT_HARD_CLIP };

    /* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
public:
    struct Result
    {
        pCigar cigar;
        uint32_t leading_deletion_bases_removed;
        uint32_t trailing_deletion_bases_removed;
    };

    /* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
private:
    pMemoryPool _pool;
    struct
    {
        uint32_t num;
        uint32_t data[128];
    } _cigar_elements;
    uint32_t _last_operator{BAM_CUNINITIALIZE};
    uint32_t _leading_deletion_bases_removed{0};
    uint32_t _trailing_deletion_bases_removed{0};
    uint32_t _trailing_deletion_bases_removed_in_make{0};
    Section _section{Section::LEFT_HARD_CLIP};
    bool _remove_deletions_at_ends;

    /* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
public:
    DISALLOW_COPY_AND_ASSIGN(CigarBuilder);

    static pCigarBuilder create(pMemoryPool pool);
    static pCigarBuilder create(bool remove_deletions_at_ends, pMemoryPool pool);

    CigarBuilder* add(uint32_t cigar_element);
    CigarBuilder* add_all(uint32_t* cigar, uint32_t cigar_num);
    pCigar make(bool allow_empty);
    pCigar make() { return make(false); }
    uint32_t leading_deletion_bases_removed() const { return _leading_deletion_bases_removed; }
    uint32_t trailing_deletion_bases_removed() const { return _trailing_deletion_bases_removed + _trailing_deletion_bases_removed_in_make; }
    Result make_and_record_deletions_removed_result();

    /* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
private:
    explicit CigarBuilder(pMemoryPool pool)
        : CigarBuilder(true, pool)
    {}
    CigarBuilder(bool remove_deletions_at_ends, pMemoryPool pool);

    bool last_two_elements_were_deletion_and_insertion();
    void advance_section_and_validate_cigar_order(uint32_t op);

    /*! @brief skip a deletion after clipping ie at the beginning of the read */
    bool skip1() const;
    /*! @brief note the edge case of a deletion following a leading insertion, which we also skip */
    bool skip2() const;

    inline bool is_middle() const { return _section == Section::MIDDLE; }
    inline bool is_left_hard_clip() const { return _section == Section::LEFT_HARD_CLIP; }
    inline bool is_left_soft_clip() const { return _section == Section::LEFT_SOFT_CLIP; }
    inline bool is_right_hard_clip() const { return _section == Section::RIGHT_HARD_CLIP; }
    inline bool is_right_soft_clip() const { return _section == Section::RIGHT_SOFT_CLIP; }
};

}  // namespace rovaca

#endif  // ROVACA_HC_CIGAR_BUILDER_H_
