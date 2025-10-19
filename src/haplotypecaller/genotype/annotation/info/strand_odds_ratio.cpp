#include "strand_odds_ratio.h"

#include "genotype_macors.h"
#include "info_data.hpp"

namespace rovaca
{

static constexpr int32_t s_min_count = 0;
static constexpr double s_pseudocount = 1;

double StrandOddsRatio::calculate_sor(const Int32Vector2D& table)
{
    double t00 = table[0][0] + s_pseudocount;
    double t01 = table[0][1] + s_pseudocount;
    double t11 = table[1][1] + s_pseudocount;
    double t10 = table[1][0] + s_pseudocount;

    double ratio = (t00 / t01) * (t11 / t10) + (t01 / t00) * (t10 / t11);

    double ref_ratio = std::min(t00, t01) / std::max(t00, t01);
    double alt_ratio = std::min(t10, t11) / std::max(t10, t11);

    return std::log(ratio) + std::log(ref_ratio) - std::log(alt_ratio);
}

void StrandOddsRatio::calculate_annotation_from_gt_field([[maybe_unused]] pGenotypesContext gc, [[maybe_unused]] pInfoData target,
                                                         [[maybe_unused]] pMemoryPool pool)
{}

void StrandOddsRatio::calculate_annotation_from_likelihoods(pVariant vc, pRALikelihoods likelihoods, pInfoData target, pMemoryPool pool)
{
    Int32Vector2D table = get_contingency_table(vc, likelihoods, s_min_count, pool);
    double sor = calculate_sor(table);
    target->set_sor(sor);
}

}  // namespace rovaca