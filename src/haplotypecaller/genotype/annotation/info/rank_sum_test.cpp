#include "rank_sum_test.h"

#include "allele_likelihoods.hpp"
#include "rovaca_logger.h"
#include "genotype.h"
#include "genotypes_context.hpp"
#include "quality_utils.h"
#include "read_record.h"
#include "utils/mann_whitney_u.h"
#include "variant.h"

namespace rovaca
{

void RankSumTest::annotate([[maybe_unused]] pRefFragment ref, pVariant vc, pRALikelihoods likelihoods, pInfoData target, pMemoryPool pool)
{
    pGenotypesContext genotypes = vc->genotype();
    if (nullptr == genotypes || genotypes->empty()) {
        return;
    }

    DoubleVector ref_quals(pool);
    DoubleVector alt_quals(pool);
    if (nullptr != likelihoods) {
        fill_quals_from_likelihood(vc, likelihoods, ref_quals, alt_quals);
    }

    if (ref_quals.empty() || alt_quals.empty()) {
        return;
    }

    pMannWhitneyU m = MannWhitneyU::create(pool);
    MannWhitneyU::pResult result = m->test(alt_quals, ref_quals, MannWhitneyU::TestType::FIRST_DOMINATES);
    double z_score = result->z;
    if (!std::isnan(z_score)) {
        set_values_to_info(target, &z_score);
    }
}

void RankSumTest::fill_quals_from_likelihood(pVariant vc, pRALikelihoods likelihoods, DoubleVector& ref_quals, DoubleVector& alt_quals)
{
    auto best_alleles = likelihoods->best_alleles_breaking_ties();
    pReadRecord read;
    pAllele allele;
    for (auto pbest : best_alleles) {
        read = pbest->evidence;
        allele = pbest->best_allele;
        if (pbest->is_informative() && is_usable_read(read, vc)) {
            OptionalDouble d = get_element_for_read(read, vc, pbest);
            if (d.first && d.second != NEGATIVE_INFINITY) {
                if (allele->is_reference()) {
                    ref_quals.emplace_back(d.second);
                }
                else if (vc->has_allele(allele)) {
                    alt_quals.emplace_back(d.second);
                }
            }
        }
    }
}

bool RankSumTest::is_usable_read(pReadRecord read, [[maybe_unused]] pVariant vc)
{
    return read->mapping_quality() != 0 && read->mapping_quality() != QualityUtils::s_mapping_quality_unavailable;
}

OptionalDouble RankSumTest::get_element_for_read(pReadRecord read, pVariant vc, [[maybe_unused]] void* best_allele)
{
    return get_element_for_read(read, vc);
}

}  // namespace rovaca