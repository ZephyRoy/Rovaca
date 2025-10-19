#ifndef ROVACA_HC_ROVACA_VARIANT_CONTEXT_UTILS_H_
#define ROVACA_HC_ROVACA_VARIANT_CONTEXT_UTILS_H_
#include <algorithm>
#include <unordered_map>

#include "allele.h"
#include "forward.h"
#include "genotype_enum.h"

namespace rovaca
{

typedef struct RepeatUnitsResult
{
    Int32Vector lengths;
    RefFragment bases;
    explicit RepeatUnitsResult(pMemoryPool pool)
        : lengths{pool}
        , bases{0, nullptr}
    {}
} RepeatUnitsResult, *pRepeatUnitsResult;

/*!
 * @brief 为了满足以下几点要求
 * 1. 无重复
 * 2. 保留插入顺序
 * 3. 需要保留old->new的映射关系
 */
struct AlleleMap
{
    std::pmr::unordered_map<pAllele, pAllele, AlleleHash, AlleleEqual> old2new_map;
    std::pmr::vector<pAllele> new_arr;

    explicit AlleleMap(pMemoryPool pool)
        : old2new_map(pool)
        , new_arr(pool)
    {
        new_arr.reserve(10);
    }

    AlleleVector remap(const AlleleVector& as) const
    {
        AlleleVector result{new_arr.get_allocator()};
        result.reserve(new_arr.size());
        std::for_each(as.begin(), as.end(), [&](pAllele a) { result.push_back(old2new_map.count(a) ? old2new_map.at(a) : a); });
        return result;
    }

    void clear()
    {
        old2new_map.clear();
        new_arr.clear();
    }
};

class Event
{
private:
    int32_t tid_;
    int64_t start_;
    int64_t stop_;
    pAllele ref_;
    pAllele alt_;

public:
    Event(int32_t tid, int64_t start, pAllele ref, pAllele alt, pMemoryPool pool)
        : tid_(tid)
        , start_(start)
        , stop_(0)
        , ref_(nullptr)
        , alt_(nullptr)
    {
        make_minimal_representation(ref, alt, pool);
    }

    void make_minimal_representation(pAllele ref, pAllele alt, pMemoryPool pool);

    bool operator==(const Event& o) const;
};

namespace ROVACAVariantContextUtils
{
/*!
 * @brief create an allele mapping for the given context where its reference allele must (potentially) be extended to the given allele
 * The refAllele is the longest reference allele seen at this start site.
 * So imagine it is:
 * ref_allele: ACGTGA
 * in_ref:     ACGT
 * in_alt:     A
 * We need to remap all of the alleles in vc to include the extra GA so that in_ref => ref_allele and in_alt => AGA
 */
void create_allele_mapping(pAllele ref_allele, pAllele in_ref, const AlleleVector& in_alts, AlleleMap& mm);
void create_allele_mapping(pAllele ref_allele, pVariant vc, AlleleMap& mm);

void resolve_incompatible_alleles(pAllele ref_allele, pVariant vc, AlleleMap& mm);

pAllele determine_reference_allele(pAllele ref1, pAllele ref2);
pAllele determine_reference_allele(const VariantVector& vcs, pInterfaceLocatable loc);
pAllele determine_reference_allele(const VariantVector& vcs);

/*! @brief sorts a vector of variant contexts by priority based on a given priority list and genotype merge type */
VariantVector sort_variant_contexts_by_priority(const VariantVector& unsorted_vcs, GenotypeMergeType g);

/*! @brief merges variant_contexts into a single hybrid. takes genotypes for common samples in priority order, if provided. */
pVariant simple_merge(const VariantVector& unsorted_vcs, FilteredRecordMergeType f, GenotypeMergeType g, bool filtered_are_uncalled);

/*!
 * @brief 此函数实质就是给vc中的GenotypesContext添加alleles
 * @param vc
 * @param default_ploidy
 * @param pool
 * @return
 */
pGenotypesContext subset_to_ref_only(pVariant vc, int32_t default_ploidy, pMemoryPool pool);

/**
 * add the genotype call (gt) field to genotype_builder using the requested {@link genotype_assignment_method}
 * @param gb the builder where we should put our newly called alleles, cannot be null
 * @param assignment_method the method to use to do the assignment, cannot be null
 * @param genotype_likelihoods a vector of likelihoods to use if the method requires pls, should be log10 likelihoods, cannot be null
 * @param alleles_to_use the alleles with respect to which the likelihoods are defined
 */
void make_genotype_call(int32_t ploidy, pGenotype g, GenotypeAssignmentMethod method, const DoubleVector& likelihoods,
                        const AlleleVector& alleles_to_use, const AlleleVector& original_gt, pGenotypePriorCalculator gpc);

bool is_informative(const DoubleVector& gls);

AlleleVector best_match_to_original_gt(const AlleleVector& alleles_to_use, const AlleleVector& original_gt);

/*!
 * @brief 对输入的 Variant 对象反向修剪
 * 对所有等位基因进行比较，找到最长的公共前缀和最长的公共后缀，并将其修剪掉。这样可以将等位基因的长度缩短，减少数据量和提高效率
 * @param untrimmed_result
 * @return
 */
pVariant reverse_trim_alleles(pVariant untrimmed_result);

/*!
 * @brief 对输入的 Variant 对象进行正向或者反向修剪操作，将其不必要的等位基因修剪掉，以减少数据量和提高效率
 * @param input_vc
 * @param trim_forward 否从正向修剪
 * @param trim_reverse 是否反向修剪
 * @return
 */
pVariant trim_alleles(pVariant input_vc, bool trim_forward, bool trim_reverse);

/*!
 * @brief 将 Variant 对象的所有等位基因前面 fwd_trim_end 个字符和后面 rev_trim 个字符修剪掉
 * @param input_vc
 * @param fwd_trim_end
 * @param rev_trim
 * @return
 */
pVariant trim_alleles(pVariant input_vc, int32_t fwd_trim_end, int32_t rev_trim);

pGenotypesContext update_genotypes_with_mapped_alleles(pGenotypesContext genotypes, const AlleleMap& mm);

/*!
 * @brief 获取一个变异位点的串联重复序列单元（tandem repeat units）的长度和重复单元的序列
 * @return 0成功，1失败
 */
bool get_num_tandem_repeat_units(pVariant vc, pRefFragment ref_bases_starting_at_vc_without_pad, pRepeatUnitsResult result,
                                 pMemoryPool pool);
bool get_num_tandem_repeat_units(pRefFragment ref_bases, pRefFragment alt_bases, pRefFragment remaining_ref_context,
                                 pRepeatUnitsResult result, pMemoryPool pool);

/*!
 * @brief 判断一个字符串是否可以表示为一个子字符串的串联。如果可以，返回重复单元的长度；如果不行，返回输入字符串的长度
 */
uint32_t find_repeated_substring(pRefFragment bases);

/*!
 * @brief 在一个字符串中查找一个重复单元出现的次数
 */
int32_t find_number_of_repetitions(pRefFragment repeat_unit, pRefFragment test_string, bool leading_repeats);
int32_t find_number_of_repetitions(pRefFragment repeat_unit_full, int32_t offset_in_repeat_unit_full, int32_t repeat_unit_length,
                                   pRefFragment test_string_full, int32_t offset_in_test_string_full, int32_t test_string_length,
                                   bool leading_repeats);

bool ref_fragment_equal(pRefFragment first, pRefFragment second);

bool equal_range(pRefFragment left, int32_t left_offset, pRefFragment right, int32_t right_offset, int32_t length);

/*!
 * @brief return the rightmost Variant in maybe_overlapping that overlaps cur_pos
 * @param cur_pos
 * @param maybe_overlapping
 * @return
 */
pVariant get_overlapping_variant_context(InterfaceLocatable* cur_pos, const VariantVector& maybe_overlapping);

AlleleVector homozygous_allele_list(pAllele a, int32_t ploidy, pMemoryPool pool);

int32_t calculate_gq_from_pls(const Int32Vector& pls);
int32_t calculate_gq_from_pls(const std::vector<int32_t>& pls);

std::pmr::vector<Event> split_variant_context_to_biallelics_event(pVariant vc, bool trim_left, pMemoryPool pool);
VariantVector split_variant_context_to_biallelics(pVariant vc, bool trim_left, pMemoryPool pool);

}  // namespace ROVACAVariantContextUtils

}  // namespace rovaca

#endif  // ROVACA_HC_ROVACA_VARIANT_CONTEXT_UTILS_H_