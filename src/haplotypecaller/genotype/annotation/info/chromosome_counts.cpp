#include "chromosome_counts.h"

#include "info_data.hpp"
#include "variant.h"

namespace rovaca
{

void ChromosomeCounts::annotate([[maybe_unused]] pRefFragment ref, pVariant vc, [[maybe_unused]] pRALikelihoods likelihoods,
                                pInfoData target, pMemoryPool pool)
{
    if (!vc->has_genotypes()) {
        return;
    }

    int32_t an = vc->get_called_chr_count();
    const AlleleVector& alt_alleles = vc->alt_alleles();
    if (0 == an || alt_alleles.empty()) {
        return;
    }

    Int32Vector ac(pool);
    DoubleVector af(pool);
    int32_t founders_alt_chromosomes;
    for (pAllele alt : alt_alleles) {
        founders_alt_chromosomes = vc->get_called_chr_count(alt);
        ac.push_back(founders_alt_chromosomes);
        af.emplace_back((double)founders_alt_chromosomes / (double)an);
    }
    target->set_chromosome_counts(an, std::move(ac), std::move(af));
}

}  // namespace rovaca