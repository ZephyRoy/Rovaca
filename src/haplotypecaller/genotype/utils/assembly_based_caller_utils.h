#ifndef ROVACA_HC_ASSEMBLY_BASED_CALLER_UTILS_H_
#define ROVACA_HC_ASSEMBLY_BASED_CALLER_UTILS_H_
#include <unordered_map>

#include "allele.h"
#include "forward.h"

namespace rovaca
{

enum PhaseGroup { PHASE_10 = 0, PHASE_01 = 1 };

typedef std::pair<int32_t, PhaseGroup> Int32PhaseGroupPair;
typedef std::pmr::map<pVariant, std::pair<int32_t, PhaseGroup>> VariantToInt32PhaseGroupPairMap;

namespace AssemblyBasedCallerUtils
{

#define REF_MODEL_DELETION_QUAL (30)
#define BASE_QUAL_THRESHOLD     (6)

/*!
 * @brief Returns a mapping from Allele in the merged, which represents all of the alleles being genotyped at loc, to a list of Haplotypes
 * that support that allele. If the merged includes a spanning deletion allele, all haplotypes that support spanning deletions will be
 * assigned to that allele in the map
 */
AlleleToHaplotypeVectorMap create_allele_mapper(const HaplotypeVector& haplotypes, pVariant merged, int64_t loc, bool emit_spanning_dels);

/*!
 * @brief 该函数返回在组装的单倍型中发现的在当前位置的变异集合
 * @param loc
 * @param haplotypes
 * @param include_spanning_events 如果参数设置为 true，则结果将包括跨越 loc 的事件
 * @return
 */
VariantVector get_variant_contexts_from_active_haplotypes(int64_t loc, const HaplotypeVector& haplotypes, bool include_spanning_events);

/*!
 * @brief 多个变异位点合并成一个 Variant 对象
 * @note 返回的 Variant 中等位基因顺序与输入的 vcs 一致
 * @param vcs
 * @return
 */
pVariant make_merged_variant_context(const VariantVector& vcs);

AlleleSet get_alleles_consistent_with_given_alleles(const AlleleVector& given_alleles, pVariant merged_vc);

void realign_reads_to_their_best_haplotype(pRHLikelihoods read_likelihoods, pHaplotype reference, int64_t padded_reference_loc_start,
                                           p_lib_sw_avx sw, pMemoryPool pool, pBamDataPool bam_pool);

/*!
 * @brief 迭代每个Read，生成pileup
 * @param un_sorted_reads
 * @param active_loc
 * @return
 */
ReadPileupVector get_pileups_over_reference(const ReadHashSet& un_sorted_reads, pSimpleInterval active_loc, pBases ref,
                                            int64_t bases_offset, pMemoryPool pool);

bool dont_include_read_in_pileup(int32_t adaptor_boundary, int64_t insert_size, bool is_reverse_strand, int64_t pos);

/*!
 * @brief 对单个等位基因进行分相(phasing)。分相是指确定一对等位基因中的每个等位基因来自哪个亲本
 * @param calls
 * @param called_haplotypes
 * @return
 */
VariantVector phase_calls(const VariantVector& calls, const HaplotypeList& called_haplotypes);

/*!
 * @brief 构建从变异位点的等位基因到包含该等位基因的单倍型的映射
 * @param calls
 * @param called_haplotypes
 * @return
 */
VariantToHaplotypeSetMap construct_haplotype_mapping(const VariantVector& calls, const HaplotypeList& called_haplotypes);

/*!
 * @brief 构建从 Variant 对象到相位集合ID的映射
 * 遍历 calls 中的每个 Variant 对象，如果该 Variant 对象不是双等位基因或不止一个位点特异性的替代等位基因，
 * 则将该 Variant 对象存储在 phase_set_mapping中，并将其对应的值设置为空pair。
 * 否则，它会找到该 Variant 对象的位点特异性的替代等位基因，并使用该等位基因从 haplotype_map 中获取包含该等位基因的 haplotype 对象集合。
 * 然后，它会将这些 haplotype 对象集合分成不同的相位集合，并为每个相位集合分配一个唯一的相位集合id。
 * 最后，它将该 Variant 对象和相位集合id存储在 phase_set_mapping 中，并返回该map。
 * @param calls
 * @param haplotype_map
 * @return
 */
VariantToInt32PhaseGroupPairMap construct_phase_set_mapping(const VariantVector& calls, const VariantToHaplotypeSetMap& haplotype_map);

/*!
 * @brief 将相位集合组装在一起，并相应地更新原始的 Variant 对象
 * 方法的主要操作包括：
 *    遍历所有相位组，对于每个相位组，找到属于该组的所有变异位点，并按照原始顺序记录它们的索引。
 *    如果一个相位组中的变异位点数小于2，则抛出异常。
 *    对于每个相位组，为其内部的所有变异位点创建一个唯一的ID，并将该相位组的第一个变异位点的位置作为相位组ID。
 *    更新每个变异位点的信息，将其相位化，并将更新后的变异位点添加到输出列表中。
 * @param calls
 * @param phase_set_mapping
 * @param index_to
 * @return
 */
VariantVector construct_phase_groups(const VariantVector& calls, const VariantToInt32PhaseGroupPairMap& phase_set_mapping,
                                     int32_t index_to);

/*!
 * @brief Create a unique identifier given the variant context
 * @param vc
 * @return
 */
pBases create_unique_id(pVariant vc);

/*!
 * @brief Add physical phase information to the provided variant context
 * @param vc
 * @param id
 * @param phase_gt
 * @param phase_set_id
 * @return
 */
pVariant phase_vc(pVariant vc, pBases id, PhaseGroup phase_gt, int32_t phase_set_id);

/*!
 * @brief 检测是否为变异等位基因
 * @param a
 * @return
 */
bool is_site_specific_alt_allele(pAllele a);

/*!
 * @brief 检测是否仅有一个合格的变异
 * @param vc
 * @return
 */
bool is_biallelic_with_one_site_specific_alternate_allele(pVariant vc);

/**
 * If at least one exists, returns a concrete (not NONREF) site-specific (starting at the current POS) alternate allele
 * from within the current variant context.
 */
pAllele get_site_specific_alternate_allele(pVariant vc);

bool contains_all(const HaplotypeSet& big, const HaplotypeSet& small);
void retain_all(HaplotypeSet& src, const HaplotypeSet& target);

/*!
 * @brief
 * @param ref ref_loc 对应的 ref
 * @param ref_loc original_padded 基础上左右扩 500
 * @param padded
 * @param pool
 * @return
 */
pHaplotype create_reference_haplotype(pRefFragment ref, pSimpleInterval ref_loc, pSimpleInterval padded, pMemoryPool pool);

}  // namespace AssemblyBasedCallerUtils

}  // namespace rovaca

#endif  // ROVACA_HC_ASSEMBLY_BASED_CALLER_UTILS_H_
