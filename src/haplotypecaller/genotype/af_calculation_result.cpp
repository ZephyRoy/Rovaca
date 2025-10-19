#include "af_calculation_result.h"

#include <algorithm>

#include "rovaca_logger.h"
#include "math_utils.h"
#include "quality_utils.h"

namespace rovaca
{

static const double s_epsilon = 1.0e-10;

pAFCalculationResult AFCalculationResult::create(double log10posterior_of_no_variant, AlleleToDoubleMap&& log10p_ref_by_allele,
                                                 Int32Vector&& allele_counts_of_mle, const AlleleVector& alleles_used_in_genotyping,
                                                 pMemoryPool pool)
{
    AlleleVector a(alleles_used_in_genotyping, pool);
    return new ALLOC_TYPE_IN_POOL(pool, AFCalculationResult)
        AFCalculationResult{log10posterior_of_no_variant, std::forward<AlleleToDoubleMap>(log10p_ref_by_allele),
                            std::forward<Int32Vector>(allele_counts_of_mle), std::move(a)};
}

double AFCalculationResult::get_log10posterior_of_allele_absent(pAllele allele)
{
    CHECK_CONDITION_EXIT(!_log10p_ref_by_allele.count(allele), "unknown allele");
    return _log10p_ref_by_allele.at(allele);
}

bool AFCalculationResult::passes_threshold(pAllele allele, double phred_scale_qual_threshold)
{
    return get_log10posterior_of_allele_absent(allele) + s_epsilon < QualityUtils::qual_to_error_prob_log10(phred_scale_qual_threshold);
}

AFCalculationResult::AFCalculationResult(double log10posterior_of_no_variant, AlleleToDoubleMap&& log10p_ref_by_allele,
                                         Int32Vector&& allele_counts_of_mle, AlleleVector&& alleles_used_in_genotyping)
    : _log10posterior_of_no_variant(log10posterior_of_no_variant)
    , _log10p_ref_by_allele(std::move(log10p_ref_by_allele))
    , _allele_counts_of_mle(std::move(allele_counts_of_mle))
    , _alleles_used_in_genotyping(std::move(alleles_used_in_genotyping))
{
    CHECK_CONDITION_EXIT(log10posterior_of_no_variant > 0.0, "log10 posterior must be a valid log probability");
    CHECK_CONDITION_EXIT(_alleles_used_in_genotyping.empty(), "alleles_used_in_genotyping must be non-null list of at least 1 value");
    CHECK_CONDITION_EXIT(_allele_counts_of_mle.size() != _alleles_used_in_genotyping.size() - 1,
                         "_allele_counts_of_mle.size() != alleles_used_in_genotyping.size() - 1");
    CHECK_CONDITION_EXIT(_log10p_ref_by_allele.size() != _alleles_used_in_genotyping.size() - 1,
                         "_log10p_ref_by_allele.size() != _alleles_used_in_genotyping.size() - 1");

    bool contains;
    for (const auto& itr : _log10p_ref_by_allele) {
        contains = false;
        for (pAllele a : _alleles_used_in_genotyping) {
            if (a->equals(*itr.first)) {
                contains = true;
                break;
            }
        }
        CHECK_CONDITION_EXIT(!contains, "log10p_ref_by_allele doesn't contain all of the alleles used in genotyping");
    }
}

double AFCalculationResult::log10prob_variant_present() const { return math_utils::log10one_minus_pow10(_log10posterior_of_no_variant); }

int32_t AFCalculationResult::get_allele_count_at_mle(pAllele allele) const
{
    CHECK_CONDITION_EXIT(allele->is_reference(), "cannot get the alt allele index for reference allele");
    int32_t index_in_all_alleles_including_ref = INVALID_INT;
    for (int32_t i = 0, len = (int32_t)_alleles_used_in_genotyping.size(); i < len; ++i) {
        if (_alleles_used_in_genotyping.at(i)->equals(*allele)) {
            index_in_all_alleles_including_ref = i;
            break;
        }
    }
    CHECK_CONDITION_EXIT(INVALID_INT == index_in_all_alleles_including_ref, "could not find allele");
    int32_t index_in_alt_alleles = index_in_all_alleles_including_ref - 1;
    return _allele_counts_of_mle.at(index_in_alt_alleles);
}

}  // namespace rovaca