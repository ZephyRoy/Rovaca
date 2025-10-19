#ifndef ROVACA_HC_GENOTYPE_UTILS_H_
#define ROVACA_HC_GENOTYPE_UTILS_H_
#include "forward.h"

namespace rovaca
{

typedef struct GenotypeCounts
{
    double ref, het, hom;
} GenotypeCounts, *pGenotypeCounts;

namespace GenotypeUtils
{

bool genotype_is_usable_for_afcalculation(pGenotype g);
bool is_diploid_with_likelihoods(pGenotype g);
bool is_called_and_diploid_with_likelihoods_or_with_gq(pGenotype g);

/*!
 * @brief Make approximate likelihoods for a diploid genotype without PLs
 * For a hom-ref, as long as we have GQ we can make a very accurate QUAL calculation
 * since the hom-var likelihood should make a minuscule contribution
 * @param g
 * @param num_alleles
 * @param pool
 * @return
 */
DoubleVector make_approximate_diploid_log10likelihoods_from_gq(pGenotype g, int32_t num_alleles, pMemoryPool pool);

/**
 * returns a triple of ref/het/hom genotype "counts".
 * @param vc the Variant that the GenotypesContext originated from, non-null
 * @param genotypes a GenotypesContext containing genotypes to count, these must be a subset of vc->genotypes()
 * @param round_contribution_from_each_genotype if this is true, the normalized likelihood from each genotype will be rounded before adding
 *                                              to the total count
 */
pGenotypeCounts compute_diploid_genotype_counts(pVariant vc, pGenotypesContext gc, bool round_contribution_from_each_genotype,
                                                pMemoryPool pool);

}  // namespace GenotypeUtils

}  // namespace rovaca

#endif  // ROVACA_HC_GENOTYPE_UTILS_H_
