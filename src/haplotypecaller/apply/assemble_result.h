#ifndef ASSEMBLE_RESULT_H
#define ASSEMBLE_RESULT_H

#include <memory_resource>
#include <unordered_map>

#include "assemble_argument.h"
#include "debug.h"
#include "genotype/forward.h"
#include "hc_assemble_main.h"

using namespace rovaca;
using HaplotypeHashMap = std::unordered_map<pHaplotype, pHaplotype>;
using HaplotypeHashMapIterator = std::unordered_map<pHaplotype, pHaplotype>::iterator;

class AssembleRegion
{
public:
    bool isActive;
    pSimpleInterval activeSpan;
    pSimpleInterval paddedSpan;
    pMemoryPool resource;

    AssembleRegion(p_hc_region_active_storage region, pMemoryPool mem);
    AssembleRegion(bool active, rovaca::pSimpleInterval original, rovaca::pSimpleInterval padded);
    AssembleRegion(const AssembleRegion& one, pMemoryPool mem);

    ~AssembleRegion() { resource->~memory_resource(); }
};

using pAssembleRegion = AssembleRegion*;
class AssembleResult
{
public:
    static AssembleResult* create(pMemoryPool resource);
    AssembleResult(const AssembleResult& one);
    ~AssembleResult();
    inline void set_reads(const std::pmr::vector<pReadRecord>& assemble_reads) { reads_ = assemble_reads; }
    inline void set_haplotypes(const std::pmr::vector<pHaplotype>& assemble_hapltypes) { hapltypes_ = assemble_hapltypes; }
    inline void set_region(pAssembleRegion assemble_region) { genotyping_region_ = assemble_region; }
    inline void set_padded_ref(uint8_t* ref_with_padded) { padded_ref_ = ref_with_padded; }
    inline void set_padded_refloc(pSimpleInterval refloc_with_padded) { padded_ref_loc_ = refloc_with_padded; }
    inline void set_ref_haplotype(pHaplotype refHaplotype) { refHaplotype_ = refHaplotype; }

    inline const std::pmr::vector<pReadRecord>& get_reads() { return reads_; }
    inline const std::pmr::vector<pHaplotype>& get_haplotypes() { return hapltypes_; }
    inline pAssembleRegion get_region() { return genotyping_region_; }
    inline uint8_t* get_padded_ref() { return padded_ref_; }
    inline pSimpleInterval get_padded_refloc() { return padded_ref_loc_; }
    inline pHaplotype refHaplotype() { return refHaplotype_; }
    bool add_result(pHaplotype result);

    /**
     * @brief 根据trim以后的region来trimAssembleResult,涉及genotyping region、单倍体以及reference.
     * @param trimmedAssemblyRegion
     * @param arguments
     * @param resource
     * @return AssembleResult*
     */
    AssembleResult* trim_to(pAssembleRegion trimmedAssemblyRegion, const AssembleArgument& arguments, pMemoryPool resource);

private:
    pMemoryPool resource_;
    uint8_t* padded_ref_;
    pHaplotype refHaplotype_;
    pSimpleInterval padded_ref_loc_;
    pAssembleRegion genotyping_region_;
    std::pmr::vector<pReadRecord> reads_;
    std::pmr::vector<pHaplotype> hapltypes_;
    std::vector<pVariant> variant_events_;
    std::vector<uint32_t> kemers_;

private:
    AssembleResult(pMemoryPool resource);
    AssembleResult(pAssembleRegion assemble_region, uint8_t* ref_with_padded, pSimpleInterval refloc_with_padded,
                   const std::pmr::vector<pReadRecord>& assemble_reads, const std::pmr::vector<pHaplotype>& assemble_hapltypes);
    /**
     * @brief
     * @param span
     * @param arguments
     * @return HaplotypeHashMap
     */
    HaplotypeHashMap calculateOriginalByTrimmedHaplotypes(pSimpleInterval span, const AssembleArgument& arguments, pMemoryPool resource);

    /**
     * @brief trim haplotypes to span merging haplotypes that are equal in bases.
     * @param span
     * @param haplotypeList
     * @return HaplotypeHashMap
     */
    HaplotypeHashMap trimDownHaplotypes(pSimpleInterval span, const std::pmr::vector<pHaplotype>& haplotypeList, pMemoryPool resource);

    /**
     * @brief 原为class Haplotype 的一个成员函数,为保持编译独立性，外部实现于此.
     * @param original
     * @param loc
     * @param ignoreRefState
     * @param pool
     * @return pHaplotype
     */
    pHaplotype trim_haplotype(pHaplotype original, pSimpleInterval loc, bool ignoreRefState, pMemoryPool pool);

    /**
     * @brief Get the bases covering ref interval object
     * 原为AlignmentUtils的一个方法,拷贝来用
     * @param ref_start
     * @param ref_end
     * @param bases
     * @param bases_start_on_ref
     * @param bases_to_ref_cigar
     * @param pool
     * @return pBases
     */
    pBases get_bases_covering_ref_interval(int64_t ref_start, int64_t ref_end, pBases bases, int64_t bases_start_on_ref,
                                           pCigar bases_to_ref_cigar, pMemoryPool pool);

    /**
     * @brief 原为AlignmentUtils的一个方法,拷贝来用
     *
     * @param cigar
     * @param start
     * @param end
     * @param pool
     * @return pCigar
     */
    pCigar trim_cigar_by_reference(pCigar cigar, int64_t start, int64_t end, pMemoryPool pool);

    /**
     * @brief 用以更新单倍体的reference信息.
     * @param hap
     */
    void updateReferenceHaplotype(pHaplotype hap);

    void mapOriginalToTrimmed(const HaplotypeHashMap& originalByTrimmedHaplotypes, HaplotypeHashMap& sortedOriginalByTrimmedHaplotypes,
                              const std::pmr::vector<pHaplotype>& trimmedHaplotypes);
};

#endif  // ASSEMBLE_RESULT_H