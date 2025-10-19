#ifndef ROVACA_HC_ALLELE_SUBSETTING_UTILS_H_
#define ROVACA_HC_ALLELE_SUBSETTING_UTILS_H_
#include "forward.h"
#include "genotype_enum.h"
#include "interface/interface_allele_list.hpp"

namespace rovaca
{

namespace AlleleSubsettingUtils
{

/*!
 * @brief 给定一组原始等位基因和一组新等位基因，找出旧PL索引，以及新PL索引之间的对应关系
 * @param ploidy
 * @param original_alleles
 * @param new_alleles
 * @return
 */
Int32Vector subsetted_plindices(int32_t ploidy, const AlleleVector& old_alleles, const AlleleVector& new_alleles,
                                InterfaceAlleleListPermutation<pAllele>* allele_permutation, pMemoryPool pool);

/*!
 * @brief 使用子集的PLs和ADs创建一个新的GenotypesContext，根据alleles_to_keep提供的排序重新排序子集的等位基因。
 * @param old 原始的GenotypesContext
 * @param default_ploidy
 * @param original_alleles 原始的等位基因
 * @param alleles_to_keep 新Genotypes要使用的子集等位基因
 * @param gpc
 * @param assignment_method PLs的分配策略
 * @param allele_based_length_annots 基于等位基因数量（A，R和G长度类型）的注释列表
 * @param pool
 * @return
 */
pGenotypesContext subset_alleles(pGenotypesContext old, int32_t default_ploidy, const AlleleVector& original_alleles,
                                 const AlleleVector& alleles_to_keep, pGenotypePriorCalculator gpc,
                                 GenotypeAssignmentMethod assignment_method, const StringSet& allele_based_length_annots, pMemoryPool pool);
pGenotypesContext subset_alleles(pGenotypesContext old, int32_t default_ploidy, const AlleleVector& original_alleles,
                                 const AlleleVector& alleles_to_keep, pGenotypePriorCalculator gpc,
                                 GenotypeAssignmentMethod assignment_method, pMemoryPool pool);

/*!
 * @brief 根据样本基因型中的数量和基因型的置信度计算每个等位基因的可能性分数
 */
AlleleVector calculate_most_likely_alleles(pVariant vc, int32_t default_ploidy, int32_t num_alt_alleles_to_keep,
                                           bool ensure_return_contains_alt, pMemoryPool pool);

/*!
 * @brief 计算变异上备用等位基因的可能性总和。
 * 对于每个备用等位基因，该方法计算所有样本中最可能的基因型包含该等位基因的样本的GL差异和最可能的基因型与同源参考基因型之间的GL差异之和。
 * 由于GL是对数似然度，因此这个数量可以表示为所有最可能的基因型包含该备用等位基因的样本的对数似然度之和与所有最可能的基因型包含同源参考基因型的样本的对数似然度
 * 之和的差值。该方法还考虑了样本的杂合性和多倍体性。
 */
DoubleVector calculate_likelihood_sums(pVariant vc, int32_t default_ploidy, bool count_alleles_without_hom_ref, pMemoryPool pool);
AlleleVector filter_to_max_number_of_alt_alleles_based_on_scores(int32_t num_alt_alleles_to_keep, const AlleleVector& alleles,
                                                                 const DoubleVector& likelihoods, pMemoryPool pool);

}  // namespace AlleleSubsettingUtils

}  // namespace rovaca

#endif  // ROVACA_HC_ALLELE_SUBSETTING_UTILS_H_
