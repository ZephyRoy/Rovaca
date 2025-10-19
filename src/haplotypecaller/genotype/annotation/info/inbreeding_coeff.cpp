#include "inbreeding_coeff.h"

#include "genotype_macors.h"
#include "genotypes_context.hpp"
#include "rovaca_logger.h"
#include "variant.h"

namespace rovaca
{

static constexpr size_t s_min_samples = 10;

void InbreedingCoeff::annotate([[maybe_unused]] pRefFragment ref, pVariant vc, [[maybe_unused]] pRALikelihoods likelihoods,
                               [[maybe_unused]] pInfoData target, [[maybe_unused]] pMemoryPool pool)
{
    pGenotypesContext genotypes = vc->genotype();
    if (nullptr == genotypes || genotypes->size() < s_min_samples || !vc->is_variant()) {
        return;
    }

    // note：hc 中 genotypes->size() 仅为 1
    RovacaLogger::error("invalid code");
    exit(EXIT_FAILURE);
}

}  // namespace rovaca