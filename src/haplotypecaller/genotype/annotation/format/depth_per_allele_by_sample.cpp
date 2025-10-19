#include "depth_per_allele_by_sample.h"

#include "genotype_macors.h"

namespace rovaca
{

void DepthPerAlleleBySample::annotate([[maybe_unused]] pRefFragment ref, pVariant vc, pGenotype g, pRALikelihoods likelihoods,
                                      pMemoryPool pool)
{
    if (nullptr == g || !g->is_called() || nullptr == likelihoods) {
        return;
    }

    const AlleleVector& alleles = vc->alleles();
    for (const auto& a : alleles) {
        CHECK_CONDITION_EXIT(!likelihoods->contains_allele(a), "Variant alleles not a subset of AlleleLikelihoods alleles");
    }

    Int32Vector ad = annotate_with_likelihoods(vc, g, alleles, likelihoods, pool);
    g->set_ad(std::move(ad));
}

}  // namespace rovaca