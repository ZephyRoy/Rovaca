#include "qual_by_depth.h"

#include <numeric>

#include "allele_likelihoods.hpp"
#include "rovaca_logger.h"
#include "genotype.h"
#include "genotypes_context.hpp"
#include "info_data.hpp"
#include "variant.h"

namespace rovaca
{

static constexpr double MAX_QD_BEFORE_FIXING = 35.0;

void QualByDepth::annotate([[maybe_unused]] pRefFragment ref, pVariant vc, pRALikelihoods likelihoods, pInfoData target,
                           [[maybe_unused]] pMemoryPool pool)
{
    if (!vc->has_log10_error()) {
        return;
    }

    pGenotypesContext genotypes = vc->genotype();
    if (nullptr == genotypes || genotypes->empty()) {
        return;
    }

    int32_t depth = get_depth(genotypes, likelihoods);
    if (depth == 0) {
        return;
    }

    double qual = -10 * vc->log10_error();
    double qd = qual / (double)depth;

    /*!
     * @note 此处原逻辑需要对超过MAX_QD_BEFORE_FIXING的qd做一个随机的调整:fix_too_high_qd
     * 新架构中因为多线程无法保证一致，不再进行随机过程
     */
    target->set_qd(std::min(qd, MAX_QD_BEFORE_FIXING));
}

int32_t QualByDepth::get_depth(pGenotypesContext gc, pRALikelihoods likelihoods)
{
    int32_t depth = 0;
    int32_t adrestricted_depth = 0;

    pGenotype g;
    for (size_t i = 0, len = gc->size(); i < len; ++i) {
        g = gc->at(i);
        if (!g->is_het() && !g->is_hom_var()) {
            continue;
        }

        if (g->has_ad()) {
            const Int32Vector& ad = g->ad();
            int32_t total_addepth = std::accumulate(ad.begin(), ad.end(), 0);
            if (total_addepth != 0) {
                if (total_addepth - ad[0] > 1) {
                    adrestricted_depth += total_addepth;
                }
                depth += total_addepth;
                continue;
            }
        }

        if (nullptr != likelihoods) {
            depth += (int32_t)likelihoods->sample_evidence_count((size_t)g->sample_id());
        }
        else if (g->has_dp()) {
            depth += g->get_dp();
        }
    }

    // if the ad-restricted depth is a usable value (i.e. not zero), then we should use that one going forward
    if (adrestricted_depth > 0) {
        depth = adrestricted_depth;
    }
    return depth;
}

}  // namespace rovaca