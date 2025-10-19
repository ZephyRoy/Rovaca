#include "excess_het.h"

#include <numeric>

#include "rovaca_logger.h"
#include "genotypes_context.hpp"
#include "info_data.hpp"
#include "utils/genotype_utils.h"
#include "variant.h"

namespace rovaca
{

static constexpr bool s_round_genotype_counts = true;
static constexpr double s_min_needed_value = 10.e-16;
static constexpr double s_phred_scaled_min_p_value = 160.00;

void ExcessHet::annotate([[maybe_unused]] pRefFragment ref, pVariant vc, [[maybe_unused]] pRALikelihoods likelihoods, pInfoData target,
                         pMemoryPool pool)
{
    pGenotypesContext genotypes = vc->genotype();
    if (nullptr == genotypes || genotypes->empty() || !vc->is_variant()) {
        return;
    }

    std::pair<int32_t, double> sample_count_eh = calculate_eh(vc, genotypes, pool);
    if (sample_count_eh.first < 1) {
        return;
    }
    target->set_excess_het(sample_count_eh.second);
}

std::pair<int32_t, double> ExcessHet::calculate_eh(pGenotypeCounts t, int32_t sample_count, pMemoryPool pool)
{
    int32_t ref_count = int32_t(t->ref);
    int32_t het_count = int32_t(t->het);
    int32_t hom_count = int32_t(t->hom);

    // if the actual phred_pval would be infinity we will probably still filter out just a very large number
    // since the method does not guarantee precision for any p-value smaller than 1e-16, we can return the phred scaled version
    double pval = exact_test(het_count, ref_count, hom_count, pool);
    if (pval < 10e-60) {
        return {sample_count, s_phred_scaled_min_p_value};
    }

    double phred_pval = -10.0 * std::log10(pval) + 0.0;
    return {sample_count, phred_pval};
}

std::pair<int32_t, double> ExcessHet::calculate_eh(pVariant vc, pGenotypesContext genotypes, pMemoryPool pool)
{
    pGenotypeCounts t = GenotypeUtils::compute_diploid_genotype_counts(vc, genotypes, s_round_genotype_counts, pool);

    int32_t sample_count = 0;
    for (size_t i = 0, len = genotypes->size(); i < len; ++i) {
        if (GenotypeUtils::is_called_and_diploid_with_likelihoods_or_with_gq(genotypes->at(i))) {
            sample_count++;
        }
    }
    return calculate_eh(t, sample_count, pool);
}

double ExcessHet::exact_test(int32_t het_count, int32_t ref_count, int32_t hom_count, pMemoryPool pool)
{
    // split into observed common allele and rare allele
    int32_t obs_hom_r = ref_count < hom_count ? ref_count : hom_count;
    int32_t obs_hom_c = ref_count < hom_count ? hom_count : ref_count;

    int32_t rare_copies = 2 * obs_hom_r + het_count;
    int32_t n = het_count + obs_hom_c + obs_hom_r;

    DoubleVector probs(rare_copies + 1, pool);

    // find (something close to the) mode for the midpoint
    int32_t mid = (int32_t)std::floor(rare_copies * (2.0 * n - rare_copies) / (2.0 * n - 1.0));
    if ((mid % 2) != (rare_copies % 2)) {
        mid++;
    }

    probs[mid] = 1.0;
    double mysum = 1.0;

    // calculate probabilities from midpoint down
    int curr_hets = mid;
    int curr_hom_r = (rare_copies - mid) / 2;
    int curr_hom_c = n - curr_hets - curr_hom_r;

    while (curr_hets >= 2) {
        double potential_prob = probs[curr_hets] * curr_hets * (curr_hets - 1.0) / (4.0 * (curr_hom_r + 1.0) * (curr_hom_c + 1.0));
        if (potential_prob < s_min_needed_value) {
            break;
        }

        probs[curr_hets - 2] = potential_prob;
        mysum = mysum + probs[curr_hets - 2];

        // 2 fewer hets means one additional hom_r and hom_c each
        curr_hets = curr_hets - 2;
        curr_hom_r = curr_hom_r + 1;
        curr_hom_c = curr_hom_c + 1;
    }

    // calculate probabilities from midpoint up
    curr_hets = mid;
    curr_hom_r = (rare_copies - mid) / 2;
    curr_hom_c = n - curr_hets - curr_hom_r;

    while (curr_hets <= rare_copies - 2) {
        double potential_prob = probs[curr_hets] * 4.0 * curr_hom_r * curr_hom_c / ((curr_hets + 2.0) * (curr_hets + 1.0));
        if (potential_prob < s_min_needed_value) {
            break;
        }

        probs[curr_hets + 2] = potential_prob;
        mysum = mysum + probs[curr_hets + 2];

        // 2 more hets means 1 fewer hom_r and hom_c each
        curr_hets = curr_hets + 2;
        curr_hom_r = curr_hom_r - 1;
        curr_hom_c = curr_hom_c - 1;
    }

    double right_pval = probs[het_count] / mysum;
    // check if we observed the highest possible number of hets
    if (het_count == rare_copies) {
        return std::max(0.0, std::min(1.0, right_pval));
    }

    double sum = std::accumulate(probs.begin() + het_count + 1, probs.end(), 0.0);
    return std::max(0.0, std::min(1.0, right_pval + sum / mysum));
}

}  // namespace rovaca