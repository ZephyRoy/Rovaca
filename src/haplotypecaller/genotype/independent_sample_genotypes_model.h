#ifndef ROVACA_HC_INDEPENDENT_SAMPLE_GENOTYPES_MODEL_H_
#define ROVACA_HC_INDEPENDENT_SAMPLE_GENOTYPES_MODEL_H_
#include <array>

// 相对顺序不要变，否则会导致另一个无法引入
// clang-format off
#include "allele_likelihoods.hpp"
#include "allele_likelihood_matrix_mapper.hpp"
// clang-format on

#include "forward.h"
#include "genotype_likelihood_calculator.h"
#include "genotype_macors.h"
#include "genotyping_likelihoods.hpp"

namespace rovaca
{

class IndependentSampleGenotypesModel
{
    /* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
public:
    DISALLOW_COPY_AND_ASSIGN(IndependentSampleGenotypesModel);
    IndependentSampleGenotypesModel() = default;

    template <typename A>
    GenotypingLikelihoods<A>* calculate_likelihoods(InterfaceAlleleList<A>* allele_list, pInterfacePloidyModel ploidy_model,
                                                    AlleleLikelihoods<pReadRecord, A>* likelihoods, pRefFragment ref, int64_t offset,
                                                    pMemoryPool pool)
    {
        un_used(ref);
        un_used(offset);

        InterfaceAlleleListPermutation<A>* permutation = likelihoods->permutation(allele_list, pool);
        AlleleLikelihoodMatrixMapper<A> allele_likelihood_matrix_mapper(permutation);

        size_t sample_count = ploidy_model->number_of_samples();
        if (sample_count == 0) {
            return nullptr;
        }

        size_t allele_count = allele_list->number_of_alleles();
        int32_t sample_ploidy, first_ploidy = ploidy_model->sample_ploidy(0);

        std::pmr::vector<pGenotypeLikelihoods> genotype_likelihoods(pool);
        genotype_likelihoods.reserve(sample_count);

        pGenotypeLikelihoodCalculator calculator = GenotypeLikelihoodCalculator::create(first_ploidy, (int32_t)allele_count, pool);

        for (size_t i = 0; i < sample_count; ++i) {
            sample_ploidy = ploidy_model->sample_ploidy(i);

            // get a new likelihoodsCalculator if this sample's ploidy differs from the previous sample's
            if (sample_ploidy != calculator->ploidy()) {
                calculator = GenotypeLikelihoodCalculator::create(sample_ploidy, (int32_t)allele_count, pool);
            }

            auto* sample_likelihoods = allele_likelihood_matrix_mapper.map_alleles(likelihoods->sample_matrix(i), pool);
            genotype_likelihoods.push_back(calculator->genotype_likelihoods(sample_likelihoods));
        }

        return GenotypingLikelihoods<A>::create(allele_list, ploidy_model, std::move(genotype_likelihoods));
    }
};

}  // namespace rovaca

#endif  // ROVACA_HC_INDEPENDENT_SAMPLE_GENOTYPES_MODEL_H_
