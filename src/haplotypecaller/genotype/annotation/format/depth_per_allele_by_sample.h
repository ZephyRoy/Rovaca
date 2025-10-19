#ifndef ROVACA_HC_DEPTH_PER_ALLELE_BY_SAMPLE_H_
#define ROVACA_HC_DEPTH_PER_ALLELE_BY_SAMPLE_H_
#include "allele.h"
#include "allele_likelihoods.hpp"
#include "rovaca_logger.h"
#include "genotype.h"
#include "interface/interface_genotype_annotation.hpp"
#include "variant.h"

namespace rovaca
{

/*!
 * @brief 计算Format列：AD
 * 在gVCF文件中，AD列是一个必需的FORMAT列，表示每个位点的等位基因深度。
 * AD的全称是Allelic Depth，表示每个等位基因在该位点的测序深度，即每个等位基因被测序的次数。
 *
 * AD列的值是一个整数数组，数组的长度等于该位点的等位基因数。
 * 例如，一个二等位基因位点的AD值为10,20，表示第一个等位基因被测序了10次，第二个等位基因被测序了20次。
 *
 * AD列的值可以用于检测等位基因的分型质量和杂合性等问题。
 * 如果一个等位基因的AD值过低，那么可能需要增加测序深度，以提高分型质量。如果一个等位基因的AD值过高，那么可能存在杂合性等问题，需要进一步检查。
 *
 * 需要注意的是，AD列是一个必需的FORMAT列，所有的gVCF文件都包含该列。
 * 在使用gVCF文件进行变异检测和分型时，AD列的值是非常重要的，可以影响变异检测的灵敏度和特异性。
 */
class DepthPerAlleleBySample : public InterfaceGenotypeAnnotation
{
    /* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
public:
    void annotate(pRefFragment ref, pVariant vc, pGenotype g, pRALikelihoods likelihoods, pMemoryPool pool) override;

    ~DepthPerAlleleBySample() override = default;

    template <typename E>
    static Int32Vector annotate_with_likelihoods(pVariant vc, pGenotype g, const AlleleVector& alleles,
                                                 AlleleLikelihoods<E, pAllele>* likelihoods, pMemoryPool pool)
    {
        std::pmr::map<pAllele, int32_t> allele_counts{pool};
        std::pmr::map<pAllele, std::pmr::vector<pAllele>> allele_subset{pool};
        std::for_each(vc->alleles().begin(), vc->alleles().end(), [&](pAllele a) {
            allele_counts.insert({a, 0});
            allele_subset.insert({a, {a}});
        });

        auto* subsetted_likelihoods = likelihoods->marginalize(alleles, allele_subset);
        auto best_alleles = subsetted_likelihoods->best_alleles_breaking_ties(g->sample_id());

        std::for_each(best_alleles.begin(), best_alleles.end(), [&](typename AlleleLikelihoods<E, pAllele>::BestAllele* b) {
            if (b->is_informative()) {
                allele_counts.at(b->best_allele) += 1;
            }
        });

        Int32Vector counts(allele_counts.size(), pool);
        counts[0] = allele_counts.at(vc->ref_allele());
        for (size_t i = 0, len = vc->allele_num() - 1; i < len; ++i) {
            counts[i + 1] = allele_counts.at(vc->alternate_allele_at(i));
        }

        return counts;
    }
};

}  // namespace rovaca

#endif  // ROVACA_HC_DEPTH_PER_ALLELE_BY_SAMPLE_H_
