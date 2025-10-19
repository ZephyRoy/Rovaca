#include "coverage.h"

#include "allele_likelihoods.hpp"
#include "genotype_macors.h"
#include "info_data.hpp"

namespace rovaca
{

void Coverage::annotate([[maybe_unused]] pRefFragment ref, [[maybe_unused]] pVariant vc, pRALikelihoods likelihoods, pInfoData target,
                        [[maybe_unused]] pMemoryPool pool)
{
    int32_t evidence_count = (int32_t)likelihoods->evidence_count();
    if (evidence_count == 0) {
        return;
    }

    target->set_dp(evidence_count);
}

}  // namespace rovaca