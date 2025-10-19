#include "allele_frequency_calculator.h"

#include <algorithm>
#include <cmath>

#include "af_calculation_result.h"
#include "allele.h"
#include "genotype.h"
#include "genotype_allele_counts.h"
#include "genotype_likelihood_calculator.h"
#include "genotypes_context.hpp"
#include "index_range.hpp"
#include "rovaca_logger.h"
#include "math_utils.h"
#include "utils/genotype_utils.h"
#include "variant.h"

namespace rovaca
{

static constexpr int s_hom_ref_genotype_index = 0;
static constexpr double s_threshold = 0.1;  // threshold_for_allele_count_convergence

pAlleleFrequencyCalculator AlleleFrequencyCalculator::create(double ref_pseudocount, double snp_pseudocount, double indel_pseudocount,
                                                             int32_t default_ploidy)
{
    return new AlleleFrequencyCalculator{ref_pseudocount, snp_pseudocount, indel_pseudocount, default_ploidy};
}

pAlleleFrequencyCalculator AlleleFrequencyCalculator::make_calculator(int32_t ploidy, double snp_heterozygosity,
                                                                      double indel_heterozygosity, double heterozygosity_standard_deviation)
{
    double ref_pseudocount = snp_heterozygosity / std::pow(heterozygosity_standard_deviation, 2);
    double snp_pseudocount = snp_heterozygosity * ref_pseudocount;
    double indel_pseudocount = indel_heterozygosity * ref_pseudocount;
    return new AlleleFrequencyCalculator{ref_pseudocount, snp_pseudocount, indel_pseudocount, ploidy};
}

pAFCalculationResult AlleleFrequencyCalculator::calculate(pVariant vc, int32_t default_ploidy) const
{
    pGenotypesContext genotype = vc->genotype();
    bool has_likelihoods = false;
    for (size_t i = 0, len = genotype->size(); i < len; ++i) {
        if (genotype->at(i)->has_likelihoods()) {
            has_likelihoods = true;
            break;
        }
    }
    CHECK_CONDITION_EXIT(!has_likelihoods, "Variant must contain at least one genotype with likelihoods");

    size_t num_alleles = vc->allele_num();
    const AlleleVector& alleles = vc->alleles();
    uint32_t ref_length = vc->ref_allele()->length();
    CHECK_CONDITION_EXIT(num_alleles <= 1,
                         "Variant at {}:{}-{} has only a single reference allele, but get_log10pnon_ref requires at least alternate allele",
                         vc->get_tid(), vc->get_start(), vc->get_stop());

    return calculate((int32_t)num_alleles, alleles, genotype, default_ploidy, ref_length);
}

AlleleFrequencyCalculator::AlleleFrequencyCalculator(double ref_pseudocount, double snp_pseudocount, double indel_pseudocount,
                                                     int32_t default_ploidy)
    : _ref_pseudocount(ref_pseudocount)
    , _snp_pseudocount(snp_pseudocount)
    , _indel_pseudocount(indel_pseudocount)
    , _default_ploidy(default_ploidy)
{}

pAFCalculationResult AlleleFrequencyCalculator::calculate(int32_t num_allele, const AlleleVector& alleles, pGenotypesContext genotypes,
                                                          int32_t default_ploidy, uint32_t ref_length) const
{
    pMemoryPool pool = alleles.get_allocator().resource();
    DoubleVector prior_pseudocounts(pool);
    prior_pseudocounts.reserve(num_allele);

    std::for_each(alleles.begin(), alleles.end(), [&](pAllele a) {
        double prior = a->is_reference() ? _ref_pseudocount : (a->length() == ref_length ? _snp_pseudocount : _indel_pseudocount);
        prior_pseudocounts.emplace_back(prior);
    });

    DoubleVector allele_counts(num_allele, alleles.get_allocator());
    double sum, flat_log10allele_frequency = -MathUtils::log10((int32_t)num_allele);  // log10(1 / num_alleles)
    DoubleVector log10allele_frequencies(num_allele, flat_log10allele_frequency, alleles.get_allocator());

    for (double diff = POSITIVE_INFINITY; diff > s_threshold;) {
        DoubleVector new_allele_counts = effective_allele_counts(genotypes, log10allele_frequencies);

        DoubleVector subtract_arr = math_utils::ebe_subtract(allele_counts, new_allele_counts);
        std::for_each(subtract_arr.begin(), subtract_arr.end(), [](double& d) { d = std::abs(d); });
        diff = *std::max_element(subtract_arr.begin(), subtract_arr.end());
        allele_counts = std::move(new_allele_counts);

        DoubleVector posterior_pseudocounts = math_utils::ebe_add(prior_pseudocounts, allele_counts);

        // first iteration uses flat prior in order to avoid local minimum where the prior + no pseudocounts gives such a low effective
        // allele frequency that it overwhelms the genotype likelihood of a real variant basically, we want a chance to get non-zero
        // pseudocounts before using a prior that's biased against a variant
        sum = std::accumulate(posterior_pseudocounts.begin(), posterior_pseudocounts.end(), 0.0);
        auto func = [&](double x) -> double { return std::log10(x / sum); };
        log10allele_frequencies = math_utils::apply_to_array(posterior_pseudocounts, func);
    }

    double log10pno_variant = 0;
    DoubleVector log10pof_zero_counts_by_allele(num_allele, pool);

    auto has_span_allele = [](pAllele a) { return StaticAllele::get_instance()->_span_del->equals(*a); };
    bool spanning_deletion_present = std::any_of(alleles.begin(), alleles.end(), has_span_allele);
    std::pmr::map<int32_t, std::pmr::vector<int32_t>> non_variant_indices_by_ploidy(pool);
    DoubleVector2D log10absent_posteriors(num_allele, DoubleVector(pool), pool);
    std::for_each(log10absent_posteriors.begin(), log10absent_posteriors.end(), [](DoubleVector& arr) { arr.reserve(20); });

    int32_t ploidy;
    pGenotype g;
    pGenotypeLikelihoodCalculator calculator;
    for (size_t i = 0, len = genotypes->size(); i < len; ++i) {
        g = genotypes->at(i);
        if (!GenotypeUtils::genotype_is_usable_for_afcalculation(g)) {
            continue;
        }
        ploidy = g->get_ploidy() == 0 ? default_ploidy : g->get_ploidy();
        calculator = GenotypeLikelihoodCalculator::create(ploidy, num_allele, pool);
        DoubleVector log10genotype_posteriors = log10normalized_genotype_posteriors(g, calculator, log10allele_frequencies);

        if (!spanning_deletion_present) {
            log10pno_variant += log10genotype_posteriors.at(s_hom_ref_genotype_index);
        }
        else {
            if (!non_variant_indices_by_ploidy.count(ploidy)) {
                int32_t span_del_index = -1;
                for (int32_t ai = 0, alen = (int32_t)alleles.size(); ai < alen; ++ai) {
                    if (alleles.at(ai)->equals(*StaticAllele::get_instance()->_span_del)) {
                        span_del_index = ai;
                        break;
                    }
                }
                auto func = [&](int32_t n) -> int32_t {
                    Int32Vector allele_count_array{{0, ploidy - n, span_del_index, n}, pool};
                    return calculator->allele_counts_to_index(allele_count_array);
                };
                Int32Vector indices = IndexRange(0, ploidy + 1).map_to_integer(func, pool);
                non_variant_indices_by_ploidy.insert({ploidy, std::move(indices)});
            }

            const Int32Vector& non_variant_indices = non_variant_indices_by_ploidy.at(ploidy);
            auto func = [&](int32_t n) -> double { return log10genotype_posteriors.at(n); };
            DoubleVector non_variant_log10posteriors = math_utils::apply_to_array(non_variant_indices, func);
            // when the only alt allele is the spanning deletion the probability that the site is non-variant may be so close to 1 that
            // finite precision error in log10sum_log10 yields a positive value, which is bogus.  thus we cap it at 0.
            log10pno_variant += std::min(0.0, math_utils::log10_sum_log10_1(non_variant_log10posteriors));
        }

        if (2 == num_allele && !spanning_deletion_present) {
            continue;
        }

        // for each allele, we collect the log10 probabilities of genotypes in which the allele is absent, then add (in log space) to get
        // the log10 probability that the allele is absent in this sample
        std::for_each(log10absent_posteriors.begin(), log10absent_posteriors.end(), [](DoubleVector& d) { d.clear(); });

        for (int32_t genotype = 0, count = calculator->genotype_count(); genotype < count; genotype++) {
            double log10genotype_posterior = log10genotype_posteriors.at(genotype);
            auto func = [&](int32_t a) { log10absent_posteriors.at(a).push_back(log10genotype_posterior); };
            calculator->genotype_allele_counts_at(genotype)->for_each_absent_allele_index(func, num_allele);
        }

        DoubleVector log10pno_allele(pool);
        log10pno_allele.reserve(log10absent_posteriors.size());
        std::for_each(log10absent_posteriors.begin(), log10absent_posteriors.end(), [&](DoubleVector& buffer) {
            // if prob of non hom ref > 1 due to finite precision, short-circuit to avoid NaN
            log10pno_allele.push_back(std::min(0.0, math_utils::log10_sum_log10_1(buffer)));
        });

        // multiply the cumulative probabilities of alleles being absent, which is addition of logs
        math_utils::add_to_array_in_place(log10pof_zero_counts_by_allele, log10pno_allele);
    }

    // for biallelic the allele-specific qual equals the variant qual, and we short-circuited the calculation above
    if (2 == num_allele && !spanning_deletion_present) {
        log10pof_zero_counts_by_allele.at(1) = log10pno_variant;
    }

    // unfortunately afcalculation_result expects integers for the mle.  we really should emit the em no-integer values
    // which are valuable (eg in combine_gvcfs) as the sufficient statistics of the dirichlet posterior on allele frequencies
    Int32Vector integer_alt_allele_counts(pool);
    integer_alt_allele_counts.reserve(num_allele - 1);
    std::for_each(allele_counts.begin() + 1, allele_counts.end(),
                  [&](double x) { integer_alt_allele_counts.push_back((int32_t)std::round(x)); });

    // skip the ref allele (index 0)
    AlleleToDoubleMap log10pref_by_allele{pool};
    for (int32_t i = 1; i < num_allele; ++i) {
        log10pref_by_allele.insert({alleles.at(i), log10pof_zero_counts_by_allele.at(i)});
    }

    return AFCalculationResult::create(log10pno_variant, std::move(log10pref_by_allele), std::move(integer_alt_allele_counts), alleles,
                                       pool);
}

DoubleVector AlleleFrequencyCalculator::effective_allele_counts(pGenotypesContext genotypes, const DoubleVector& log10allele_frequencies)
{
    pMemoryPool pool = log10allele_frequencies.get_allocator().resource();
    int32_t num_alleles = (int32_t)log10allele_frequencies.size();
    DoubleVector log10result(num_alleles, NEGATIVE_INFINITY, log10allele_frequencies.get_allocator());

    int32_t genotype_count;
    pGenotype g;
    pGenotypeLikelihoodCalculator calculator;
    for (size_t i = 0, len = genotypes->size(); i < len; ++i) {
        g = genotypes->at(i);
        if (!GenotypeUtils::genotype_is_usable_for_afcalculation(g)) {
            continue;
        }
        calculator = GenotypeLikelihoodCalculator::create(g->get_ploidy(), num_alleles, pool);

        DoubleVector log10genotypeposteriors = log10normalized_genotype_posteriors(g, calculator, log10allele_frequencies);

        genotype_count = calculator->genotype_count();
        IndexRange(0, genotype_count).for_each([&](int32_t genotype_index) {
            calculator->genotype_allele_counts_at(genotype_index)
                ->for_each_allele_index_and_count([&](int32_t allele_index, int32_t count) {
                    log10result.at(allele_index) = math_utils::log10_sum_log10_1(
                        log10result.at(allele_index), log10genotypeposteriors.at(genotype_index) + MathUtils::log10(count));
                });
        });
    }

    math_utils::apply_to_array_in_place(log10result, [](double x) { return std::pow(10.0, x); });

    return log10result;
}

DoubleVector AlleleFrequencyCalculator::log10normalized_genotype_posteriors(pGenotype g, pGenotypeLikelihoodCalculator calculator,
                                                                            const DoubleVector& log10allele_frequencies)
{
    pMemoryPool pool = log10allele_frequencies.get_allocator().resource();
    DoubleVector log10likelihoods(pool);
    if (g->has_likelihoods()) {
        log10likelihoods = std::move(g->get_likelihoods()->_log10likelihoods);
    }
    else if (g->is_hom_ref()) {
        CHECK_CONDITION_EXIT(2 != g->get_ploidy(), "likelihoods are required to calculate posteriors for hom-refs with ploidy != 2");
        if (g->has_gq()) {
            int32_t num_alleles = (int32_t)log10allele_frequencies.size();
            log10likelihoods = GenotypeUtils::make_approximate_diploid_log10likelihoods_from_gq(g, num_alleles, pool);
        }
        else {
            RovacaLogger::error("does not contain GQ necessary to calculate posteriors");
            exit(EXIT_FAILURE);
        }
    }
    else {
        RovacaLogger::error("does not contain likelihoods necessary to calculate posteriors");
        exit(EXIT_FAILURE);
    }
    int32_t genotype_count = calculator->genotype_count();
    DoubleVector log10posteriors = IndexRange(0, genotype_count)
                                       .map_to_double(
                                           [&](int32_t i) -> double {
                                               pGenotypeAlleleCounts gac = calculator->genotype_allele_counts_at(i);
                                               return gac->log10combination_count() + log10likelihoods.at(i) +
                                                      gac->sum_over_allele_indices_and_counts([&](int32_t idx, int32_t count) {
                                                          return count * log10allele_frequencies.at(idx);
                                                      });
                                           },
                                           pool);

    return math_utils::normalize_log10(log10posteriors);
}

}  // namespace rovaca