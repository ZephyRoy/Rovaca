#include "genotype_utils.h"

#include <cmath>

#include "rovaca_logger.h"
#include "genotype.h"
#include "genotype_likelihood_calculator.h"
#include "genotypes_context.hpp"
#include "math_utils.h"
#include "quality_utils.h"
#include "variant.h"

namespace rovaca
{

static constexpr int32_t s_typical_base_quality = 30;
static int32_t s_ploidy_2_hom_var_scale_factor = (int32_t)std::round(s_typical_base_quality / -10.0 / std::log10(0.5));

bool GenotypeUtils::genotype_is_usable_for_afcalculation(pGenotype g)
{
    return g->has_likelihoods() || (g->is_hom_ref() && g->has_gq() && 2 == g->get_ploidy());
}

bool GenotypeUtils::is_diploid_with_likelihoods(pGenotype g)
{
    CHECK_CONDITION_EXIT(nullptr == g, "nullptr == g");
    return g->has_likelihoods() && 2 == g->get_ploidy();
}

bool GenotypeUtils::is_called_and_diploid_with_likelihoods_or_with_gq(pGenotype g)
{
    CHECK_CONDITION_EXIT(nullptr == g, "nullptr == g");
    return g->is_called() && 2 == g->get_ploidy() && (g->has_likelihoods() || g->has_gq());
}

DoubleVector GenotypeUtils::make_approximate_diploid_log10likelihoods_from_gq(pGenotype g, int32_t num_alleles, pMemoryPool pool)
{
    CHECK_CONDITION_EXIT(2 != g->get_ploidy(), "This method can only be used to approximate likelihoods for diploid genotypes");
    CHECK_CONDITION_EXIT(!g->has_gq(), "Genotype must have GQ in order to approximate PLs");
    Int32Vector per_sample_indexes_of_relevant_alleles(num_alleles, 1, pool);
    per_sample_indexes_of_relevant_alleles[0] = 0;

    int32_t gq = g->get_gq();
    int32_t ploidy = g->get_ploidy();
    Int32Vector approx_likelihoods{{0, gq, s_ploidy_2_hom_var_scale_factor * gq}, pool};

    pGenotypeLikelihoodCalculator calculator = GenotypeLikelihoodCalculator::create(ploidy, num_alleles, pool);
    Int32Vector genotype_index_map_by_ploidy = calculator->genotype_index_map(per_sample_indexes_of_relevant_alleles);

    Int32Vector pls(genotype_index_map_by_ploidy.size(), pool);

    for (size_t i = 0, len = pls.size(); i < len; ++i) {
        pls.at(i) = approx_likelihoods.at(genotype_index_map_by_ploidy.at(i));
    }

    return std::move(GenotypeLikelihoods::create(pls, pool)->_log10likelihoods);
}

pGenotypeCounts GenotypeUtils::compute_diploid_genotype_counts(pVariant vc, pGenotypesContext gc,
                                                               bool round_contribution_from_each_genotype, pMemoryPool pool)
{
    CHECK_CONDITION_EXIT(nullptr == vc, "nullptr == vc");
    CHECK_CONDITION_EXIT(nullptr == gc, "nullptr == gc");

    int32_t idx_aa = 0;
    int32_t idx_ab = 1;
    int32_t idx_bb = 2;
    double genotype_with_two_refs_count = 0;  // i.e. 0/0
    double genotypes_with_one_ref_count = 0;  // e.g. 0/1, 0/2, etc.
    double genotypes_with_no_refs_count = 0;  // e.g. 1/1, 1/2, 2/2, etc.

    pGenotype g;
    for (size_t i = 0, len = gc->size(); i < len; ++i) {
        g = gc->at(i);
        if (!is_diploid_with_likelihoods(g) && !is_called_and_diploid_with_likelihoods_or_with_gq(g)) {
            continue;
        }

        if (!g->has_likelihoods() && g->is_hom_ref()) {
            if (round_contribution_from_each_genotype) {
                genotype_with_two_refs_count += 1;
                continue;
            }
            else if (g->get_gq() == 0) {
                genotype_with_two_refs_count += 1.0 / 3;
                genotypes_with_one_ref_count += 1.0 / 3;
                genotypes_with_no_refs_count += 1.0 / 3;
                continue;
            }
            else {
                double gq = (double)g->get_gq();
                genotype_with_two_refs_count += QualityUtils::qual_to_prob(gq);
                genotypes_with_one_ref_count += 1 - QualityUtils::qual_to_prob(gq);
                // assume last likelihood is negligible
                continue;
            }
        }

        CHECK_CONDITION_EXIT(!g->has_likelihoods(), "Genotype has no likelihoods");
        DoubleVector normalized_likelihoods = math_utils::normalize_from_log10to_linear_space(g->get_likelihoods()->_log10likelihoods);

        DoubleVector biallelic_likelihoods{pool};
        if (vc->alt_alleles().size() > 1) {
            int32_t max_index = (int32_t)math_utils::max_element_index(normalized_likelihoods);
            const GenotypeLikelihoodsAllelePair& allele_pair = GenotypeLikelihoods::get_allele_pair(max_index);
            if (allele_pair.allele_index1 != 0 && allele_pair.allele_index2 != 0) {
                // all likelihoods go to genotypes_with_no_refs_count because no ref allele is called
                genotypes_with_no_refs_count++;
                continue;
            }

            double max_likelihood = normalized_likelihoods.at(idx_ab);
            int32_t het_index = idx_ab;
            int32_t var_index = idx_bb;

            const AlleleVector& alt_alleles = vc->alt_alleles();
            for (pAllele curr_alt : alt_alleles) {
                Int32Vector idx_vector = vc->get_gl_indices_of_alternate_allele(curr_alt, pool);
                int32_t temp_index = idx_vector.at(1);
                if (normalized_likelihoods.at(temp_index) > max_likelihood) {
                    max_likelihood = normalized_likelihoods.at(temp_index);
                    het_index = temp_index;
                    var_index = idx_vector.at(2);
                }
            }
            DoubleVector value{
                {normalized_likelihoods.at(idx_aa), normalized_likelihoods.at(het_index), normalized_likelihoods.at(var_index)}, pool};
            biallelic_likelihoods = math_utils::normalize_sum_to_one(value);
        }
        else {
            biallelic_likelihoods = std::move(normalized_likelihoods);
        }

        double ref_likelihood = biallelic_likelihoods.at(idx_aa);
        double het_likelihood = biallelic_likelihoods.at(idx_ab);
        double var_likelihood = biallelic_likelihoods.at(idx_bb);
        // note: rounding is special cased for [0,0,x] and [x,0,0] pls because likelihoods can come out as [0.5, 0.5, 0] and both counts
        // round up
        if (round_contribution_from_each_genotype) {
            genotype_with_two_refs_count += math_utils::fast_round(ref_likelihood);
            if (ref_likelihood != het_likelihood) {
                // if gq = 0 (specifically [0,0,x] pls) count as hom_ref and don't add to the other counts
                genotypes_with_one_ref_count += math_utils::fast_round(het_likelihood);
            }
            if (var_likelihood != het_likelihood) {
                // if gq = 0 (specifically [x,0,0] pls) count as het and don't count as variant
                // need specific genotypes_with_no_refs_count (rather than complement of the others) for pl[0,0,0] case
                genotypes_with_no_refs_count += math_utils::fast_round(var_likelihood);
            }
        }
        else {
            genotype_with_two_refs_count += ref_likelihood;
            genotypes_with_one_ref_count += het_likelihood;
            genotypes_with_no_refs_count += var_likelihood;
        }
    }
    return new ALLOC_TYPE_IN_POOL(pool, GenotypeCounts)
        GenotypeCounts{genotype_with_two_refs_count, genotypes_with_one_ref_count, genotypes_with_no_refs_count};
}

}  // namespace rovaca
