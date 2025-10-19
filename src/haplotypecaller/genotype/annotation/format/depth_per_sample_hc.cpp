#include "depth_per_sample_hc.h"

#include "allele_likelihoods.hpp"
#include "genotype.h"
#include "rovaca_logger.h"
#include "variant.h"

namespace rovaca
{

void DepthPerSampleHC::annotate([[maybe_unused]] pRefFragment ref, pVariant vc, pGenotype g, pRALikelihoods likelihoods, pMemoryPool pool)
{
    if (!g->is_called()) {
        return;
    }

    int32_t sample_index = g->sample_id();
    if (likelihoods->sample_evidence_count((size_t)sample_index) == 0) {
        g->set_dp(0);
        return;
    }

    const AlleleVector& alleles = vc->alleles();
    for (const auto& a : alleles) {
        if (!likelihoods->contains_allele(a)) {
            RovacaLogger::warn("Variant alleles not a subset of AlleleLikelihoods alleles");
            return;
        }
    }

    std::pmr::map<pAllele, std::pmr::vector<pAllele>> allele_subset{pool};
    std::for_each(alleles.begin(), alleles.end(), [&](pAllele a) { allele_subset.insert({a, {a}}); });
    auto* subsetted_likelihoods = likelihoods->marginalize(alleles, allele_subset);

    auto best_alleles = subsetted_likelihoods->best_alleles_breaking_ties(g->sample_id());

    int32_t depth = 0;
    std::for_each(best_alleles.begin(), best_alleles.end(), [&](AlleleLikelihoods<pReadRecord, pAllele>::BestAllele* b) {
        if (b->is_informative()) {
            depth += 1;
        }
    });

    g->set_dp(depth);
}

}  // namespace rovaca