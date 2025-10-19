#include "hom_ref_block.h"

#include <cstring>

#include "rovaca_logger.h"
#include "genotype.h"
#include "genotypes_context.hpp"
#include "utils/adapter_utils.h"
#include "utils/rovaca_variant_context_utils.h"
#include "variant.h"

#define HOM_PL_COUNT (3)

namespace rovaca
{

HomRefBlock::HomRefBlock(pVariant start_vc, int32_t lower_gq_bound, int32_t upper_gq_bound, int32_t default_ploidy)
    : tid_(start_vc->get_tid())
    , start_(start_vc->get_start())
    , end_(start_vc->get_stop())
    , ref_allele_(start_vc->ref_allele())
    , alt_allele_(start_vc->biallelic_alt())
    , lower_gq_bound_(lower_gq_bound)
    , upper_gq_bound_(upper_gq_bound)
    , ploidy_(start_vc->get_max_ploidy(default_ploidy))
    , min_gq_(INVALID_INT)
    , min_dp_(INT32_MAX)
    , dps_()
    , min_pls_(HOM_PL_COUNT)
{
    pGenotype g = start_vc->genotype()->at(0);
    if (g->has_dp()) {
        dps_.insert(std::max(0, g->get_dp()));
        min_dp_ = std::min(min_dp_, std::max(0, g->get_dp()));
    }
    if (g->has_likelihoods()) {
        min_pls_.clear();
        std::copy(g->pl().begin(), g->pl().end(), std::back_inserter(min_pls_));
    }
    if (g->has_gq()) {
        min_gq_ = g->get_gq();
    }
}

bool HomRefBlock::is_contiguous(pVariant vc) const { return vc->within_distance_of(*this, 1); }

void HomRefBlock::add(int64_t pos, int64_t end, pGenotype g)
{
    un_used(pos);
    // CHECK_CONDITION_EXIT(pos != end, "pos != end");
    CHECK_CONDITION_EXIT(g->get_ploidy() != ploidy_, "cannot add a genotype with a different ploidy");
    CHECK_CONDITION_EXIT(!within_bounds(std::min(g->get_gq(), MAX_GENOTYPE_QUAL)), "gq is not in this range");
    if (min_pls_.empty()) {
        std::copy(g->pl().begin(), g->pl().end(), std::back_inserter(min_pls_));
    }
    else {
        if (g->has_likelihoods()) {
            const Int32Vector &current_pls = g->pl();
            CHECK_CONDITION_EXIT(current_pls.size() != min_pls_.size(), "trying to merge different PL array sizes");
            for (size_t i = 0, len = current_pls.size(); i < len; ++i) {
                min_pls_[i] = std::min(min_pls_[i], current_pls[i]);
            }
        }
    }

    if (!min_pls_.empty()) {
        min_gq_ = ROVACAVariantContextUtils::calculate_gq_from_pls(min_pls_);
    }
    else {
        min_gq_ = min_gq_ == INVALID_INT ? g->get_gq() : std::min(min_gq_, g->get_gq());
    }

    end_ = end;
    if (g->has_dp()) {
        dps_.insert(std::max(0, g->get_dp()));
        min_dp_ = std::min(min_dp_, std::max(0, g->get_dp()));
    }
}

bcf1_t *HomRefBlock::to_variant(bcf_hdr_t *vcf_herder, bam_hdr_t *bam_herder, [[maybe_unused]] int32_t sample_id, bool floor_blocks,
                                bcf1_t *b)
{
    int32_t gq = !floor_blocks ? min_gq() : gq_lowerbound();
    int32_t dp = median_dp();
    return AdapterUtils::variant2bcf(bam_herder, vcf_herder, tid_, start_, ref_allele_, int32_t(end_), dp, gq, min_dp_, min_pls_, b);
}

int32_t HomRefBlock::median_dp() const
{
    auto it = dps_.begin();
    advance(it, dps_.size() / 2);
    if (dps_.size() & 1) {
        return *it;
    }
    else {
        auto it2 = it;
        it2--;
        return int32_t(std::round((*it + *it2) / 2.0));
    }
}

}  // namespace rovaca