#include "allele_subsetting_utils.h"

#include <algorithm>

#include "allele.h"
#include "rovaca_logger.h"
#include "genotype.h"
#include "genotype_allele_counts.h"
#include "genotype_likelihood_calculator.h"
#include "genotype_likelihoods.h"
#include "genotypes_context.hpp"
#include "indexed_allele_list.hpp"
#include "rovaca_variant_context_utils.h"
#include "math_utils.h"
#include "variant.h"

namespace rovaca
{

#define PL_INDEX_OF_HOM_REF (0)

Int32Vector AlleleSubsettingUtils::subsetted_plindices(int32_t ploidy, const AlleleVector& old_alleles, const AlleleVector& new_alleles,
                                                       InterfaceAlleleListPermutation<pAllele>* allele_permutation, pMemoryPool pool)
{
    int32_t genotype_count = GenotypeLikelihoods::num_likelihoods((int32_t)new_alleles.size(), ploidy);
    Int32Vector result(genotype_count, pool);

    Int32Vector new_allele_counts(pool);
    new_allele_counts.reserve(new_alleles.size() * 2);

    bool contains_only_new_alleles;
    pGenotypeAlleleCounts old_allele_counts;
    pGenotypeLikelihoodCalculator calculator = GenotypeLikelihoodCalculator::create(ploidy, (int32_t)old_alleles.size(), pool);
    for (int32_t old_plindex = 0, pl_count = calculator->genotype_count(); old_plindex < pl_count; old_plindex++) {
        old_allele_counts = calculator->genotype_allele_counts_at(old_plindex);

        contains_only_new_alleles = true;
        for (int32_t j = 0, distinct_allele_count = old_allele_counts->distinct_allele_count(); j < distinct_allele_count; ++j) {
            if (!allele_permutation->is_kept(old_allele_counts->allele_index_at(j))) {
                contains_only_new_alleles = false;
                break;
            }
        }

        if (contains_only_new_alleles) {
            // make an array in the format described in {@link GenotypeAlleleCounts}:
            // [(new) index of first allele, count of first allele, (new) index of second allele, count of second allele. . .]
            new_allele_counts.clear();
            for (int32_t new_allele_index = 0, len = (int32_t)new_alleles.size(); new_allele_index < len; ++new_allele_index) {
                new_allele_counts.emplace_back(new_allele_index);
                new_allele_counts.emplace_back(old_allele_counts->allele_count_for(allele_permutation->from_index(new_allele_index)));
            }
            int32_t new_plindex = calculator->allele_counts_to_index(new_allele_counts);
            result.at(new_plindex) = old_plindex;
        }
    }

    return result;
}

pGenotypesContext AlleleSubsettingUtils::subset_alleles(pGenotypesContext old, int32_t default_ploidy, const AlleleVector& original_alleles,
                                                        const AlleleVector& alleles_to_keep, pGenotypePriorCalculator gpc,
                                                        GenotypeAssignmentMethod assignment_method,
                                                        const StringSet& allele_based_length_annots, pMemoryPool pool)
{
    un_used(allele_based_length_annots);

    CHECK_CONDITION_EXIT(nullptr == old, "original pGenotypesContext must not be nullptr");
    CHECK_CONDITION_EXIT(alleles_to_keep.empty(), "must keep at least one allele");
    CHECK_CONDITION_EXIT(!alleles_to_keep.front()->is_reference(), "first allele must be the reference allele");

    pGenotypesContext new_genotypes = GenotypesContext::create(pool);

    auto* old_list = IndexedAlleleList<pAllele>::create(original_alleles, pool);
    auto* new_list = IndexedAlleleList<pAllele>::create(alleles_to_keep, pool);
    auto* allele_permutation = old_list->permutation(new_list, pool);

    std::pmr::map<int32_t, std::pmr::vector<int32_t>> subsetted_likelihood_indices_by_ploidy{pool};
    pGenotype g;
    int32_t ploidy;
    for (size_t gi = 0, glen = old->size(); gi < glen; ++gi) {
        g = old->at(gi);
        ploidy = g->get_ploidy() > 0 ? g->get_ploidy() : default_ploidy;
        if (!subsetted_likelihood_indices_by_ploidy.count(ploidy)) {
            Int32Vector plindices = subsetted_plindices(ploidy, original_alleles, alleles_to_keep, allele_permutation, pool);
            subsetted_likelihood_indices_by_ploidy.insert({ploidy, std::move(plindices)});
        }
        const Int32Vector& subsetted_likelihood_indices = subsetted_likelihood_indices_by_ploidy.at(ploidy);
        int32_t expected_num_likelihoods = GenotypeLikelihoods::num_likelihoods((int32_t)original_alleles.size(), ploidy);
        // create the new likelihoods array from the alleles we are allowed to use
        DoubleVector new_likelihoods(pool);
        double new_log10gq = NEGATIVE_INFINITY;
        if (g->has_likelihoods()) {
            DoubleVector original_likelihoods = std::move(g->get_likelihoods()->_log10likelihoods);
            if (ROVACA_LIKELY((int32_t)original_likelihoods.size() == expected_num_likelihoods)) {
                new_likelihoods.reserve(subsetted_likelihood_indices.size());
                std::for_each(subsetted_likelihood_indices.begin(), subsetted_likelihood_indices.end(),
                              [&](int32_t idx) { new_likelihoods.emplace_back(original_likelihoods.at(idx)); });
                new_likelihoods = math_utils::scale_log_space_array_for_numerical_stability(new_likelihoods);
            }

            if (!new_likelihoods.empty()) {
                if (new_likelihoods.size() > 1) {
                    auto max_ele_itr = std::max_element(new_likelihoods.begin(), new_likelihoods.end());
                    int32_t pl_index = (int32_t)std::distance(new_likelihoods.begin(), max_ele_itr);
                    new_log10gq = GenotypeLikelihoods::get_gq_log10from_likelihoods(pl_index, new_likelihoods, pool);
                }
                else {
                    new_log10gq = g->get_gq() / -10.0;
                }
            }
        }
        else if (g->has_gq()) {
            new_log10gq = -0.1 * g->get_gq();
        }

        pGenotype new_g = Genotype::create(pool);
        if (new_log10gq != NEGATIVE_INFINITY && g->has_gq()) {
            new_g->set_log10perror(new_log10gq);
        }
        new_g->set_pl(GenotypeLikelihoods::gls_to_pls(new_likelihoods));
        new_g->set_id(g->sample_id());

        ROVACAVariantContextUtils::make_genotype_call(ploidy, new_g, assignment_method, new_likelihoods, alleles_to_keep, g->alleles(), gpc);

        // todo: 以下代码hc暂未使用
        // restrict ad to the new allele subset
        // if(g.has_ad()) {
        //     final int[] new_ad = get_new_allele_based_read_count_annotation(alleles_to_keep, allele_permutation, g.get_ad());
        //     gb.ad(new_ad);
        //     // if we have recalculated ad and the original genotype had af but was then removed, then recalculate af based on ad counts
        //     if (allele_based_length_annots.contains(gatkvcfconstants.allele_fraction_key) &&
        //     g.has_extended_attribute(gatkvcfconstants.allele_fraction_key)) {
        //         attributes_to_remove.remove(gatkvcfconstants.allele_fraction_key);
        //         final double[] new_afs = math_utils.normalize_sum_to_one(arrays.stream(new_ad).map_to_double(x -> x).to_array());
        //         gb.attribute(gatkvcfconstants.allele_fraction_key, arrays.copy_of_range(new_afs, 1, new_afs.length)); //omit the first
        //         entry of the array corresponding to the reference
        //     }
        // }
        // if (attributes_to_remove.size() > 0) {
        //     attributes_removed_one_shot_logger.warn(() -> "the following attributes have been removed at sites where alleles were subset:
        //     " + string.join(",", attributes_to_remove));
        // }
        new_genotypes->add(new_g);
    }
    return new_genotypes;
}

pGenotypesContext AlleleSubsettingUtils::subset_alleles(pGenotypesContext old, int32_t default_ploidy, const AlleleVector& original_alleles,
                                                        const AlleleVector& alleles_to_keep, pGenotypePriorCalculator gpc,
                                                        GenotypeAssignmentMethod assignment_method, pMemoryPool pool)
{
    return subset_alleles(old, default_ploidy, original_alleles, alleles_to_keep, gpc, assignment_method, {}, pool);
}

AlleleVector AlleleSubsettingUtils::calculate_most_likely_alleles(pVariant vc, int32_t default_ploidy, int32_t num_alt_alleles_to_keep,
                                                                  bool ensure_return_contains_alt, pMemoryPool pool)
{
    CHECK_CONDITION_EXIT(nullptr == vc, "nullptr == vc");

    bool has_symbolic_non_ref = vc->has_allele(StaticAllele::get_instance()->_non_ref_allele.get());
    int32_t number_of_alleles_that_arent_proper_alts = has_symbolic_non_ref ? 2 : 1;
    int32_t number_of_proper_alt_alleles = int32_t(vc->allele_num()) - number_of_alleles_that_arent_proper_alts;

    const AlleleVector& alleles = vc->alleles();
    if (num_alt_alleles_to_keep >= number_of_proper_alt_alleles) {
        return alleles;
    }

    DoubleVector likelihoods_sums = calculate_likelihood_sums(vc, default_ploidy, ensure_return_contains_alt, pool);

    return filter_to_max_number_of_alt_alleles_based_on_scores(num_alt_alleles_to_keep, alleles, likelihoods_sums, pool);
}

DoubleVector AlleleSubsettingUtils::calculate_likelihood_sums(pVariant vc, int32_t default_ploidy, bool count_alleles_without_hom_ref,
                                                              pMemoryPool pool)
{
    DoubleVector likelihood_sums{vc->allele_num(), pool};

    pGenotypesContext genotypes = vc->genotype();
    size_t begin = count_alleles_without_hom_ref ? 1 : 0;

    for (size_t i = 0, len = genotypes->size(); i < len; ++i) {
        pGenotype g = genotypes->at(i);
        pGenotypeLikelihoods gls = g->get_likelihoods();
        if (nullptr == gls) {
            continue;
        }
        DoubleVector& gls_vector = gls->_log10likelihoods;
        size_t index_of_most_likely_variant_genotype = math_utils::max_element_index(gls_vector, begin, gls_vector.size());

        double best_gl = gls_vector.at(index_of_most_likely_variant_genotype);
        double second_gl = gls_vector.at(PL_INDEX_OF_HOM_REF);
        double gl_diff_between_ref_and_best_variant_genotype = std::abs(best_gl - second_gl);
        int32_t ploidy = g->get_ploidy() > 0 ? g->get_ploidy() : default_ploidy;
        auto allele_count = int32_t(vc->allele_num());

        pGenotypeLikelihoodCalculator calculator = GenotypeLikelihoodCalculator::create(ploidy, allele_count, pool);
        Int32Vector allele_counts = calculator->genotype_allele_counts_at(int32_t(index_of_most_likely_variant_genotype))
                                        ->allele_counts_by_index(allele_count - 1, pool);
        for (size_t allele = 1, allele_counts_len = allele_counts.size(); allele < allele_counts_len; allele++) {
            if (allele_counts[allele] > 0) {
                likelihood_sums[allele] += gl_diff_between_ref_and_best_variant_genotype;
            }
        }
    }

    return likelihood_sums;
}

AlleleVector AlleleSubsettingUtils::filter_to_max_number_of_alt_alleles_based_on_scores(int32_t num_alt_alleles_to_keep,
                                                                                        const AlleleVector& alleles,
                                                                                        const DoubleVector& likelihoods, pMemoryPool pool)
{
    auto find_result = std::find(alleles.begin(), alleles.end(), StaticAllele::get_instance()->_non_ref_allele.get());
    size_t non_ref_alt_allele_index = find_result - alleles.begin();
    size_t num_alleles = alleles.size();

    Int32Vector indexes{pool};
    indexes.reserve(num_alleles - 2);
    for (size_t i = 1; i < num_alleles; ++i) {
        if (i == non_ref_alt_allele_index) continue;
        indexes.push_back(int32_t(i));
    }
    std::sort(indexes.begin(), indexes.end(), [&](int32_t l, int32_t r) { return likelihoods.at(l) > likelihoods.at(r); });
    Int32Vector proper_alt_indexes_to_keep{indexes.begin(), indexes.begin() + num_alt_alleles_to_keep, pool};

    AlleleVector result{pool};
    result.reserve(num_alt_alleles_to_keep + 2);
    for (size_t i = 0; i < num_alleles; i++) {
        if (i == 0 || i == non_ref_alt_allele_index ||
            std::find(proper_alt_indexes_to_keep.begin(), proper_alt_indexes_to_keep.end(), i) != proper_alt_indexes_to_keep.end()) {
            result.push_back(alleles[i]);
        }
    }

    return result;
}

}  // namespace rovaca