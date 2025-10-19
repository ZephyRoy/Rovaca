#include "rovaca_variant_context_utils.h"

#include <algorithm>
#include <cstring>
#include <memory_resource>
#include <numeric>

#include "allele.h"
#include "genotype_allele_counts.h"
#include "genotype_likelihood_calculator.h"
#include "genotype_macors.h"
#include "genotype_struct.h"
#include "genotypes_context.hpp"
#include "index_range.hpp"
#include "rovaca_logger.h"
#include "math_utils.h"
#include "utils/alignment_utils.h"
#include "variant.h"

namespace rovaca
{

enum RepeatUnits { k_first = 0, k_second, k_total };

static constexpr double s_sum_gl_thresh_nocall = -0.1;

static bool is_non_symbolic_extendable_allele(pAllele a)
{
    return !(a->is_reference() || a->is_symbolic() || a->equals(*StaticAllele::get_instance()->_span_del));
}

void ROVACAVariantContextUtils::create_allele_mapping(pAllele ref_allele, pAllele in_ref, const AlleleVector& in_alts, AlleleMap& mm)
{
    uint32_t ref_allele_len = ref_allele->length();
    uint32_t in_allele_len = in_ref->length();
    CHECK_CONDITION_EXIT(ref_allele_len < in_allele_len, "bug: ref_allele->length() < in_ref->length()");

    if (ref_allele_len == in_allele_len) {
        std::for_each(in_alts.begin(), in_alts.end(), [&](pAllele a) {
            if (!mm.old2new_map.count(a)) {
                mm.old2new_map.insert({a, a});
                mm.new_arr.emplace_back(a);
            }
        });
        return;
    }

    pMemoryPool pool = mm.new_arr.get_allocator().resource();
    pBases ref_bases = ref_allele->get_display_string();
    uint32_t current_len, extra_len = ref_allele_len - in_allele_len;
    for (pAllele a : in_alts) {
        if (is_non_symbolic_extendable_allele(a)) {
            current_len = a->length();
            auto new_base = new ALLOC_FLEXIBLE_IN_POOL(pool, Bases, current_len + extra_len, uint8_t) Bases{current_len + extra_len};
            memcpy(new_base->data, a->get_display_string()->data, current_len * sizeof(uint8_t));
            memcpy(new_base->data + current_len, ref_bases->data + in_allele_len, extra_len * sizeof(uint8_t));
            pAllele extend_allele = Allele::create_allele(new_base, a->is_reference(), pool);
            if (!mm.old2new_map.count(extend_allele)) {
                mm.old2new_map.insert({a, extend_allele});
                mm.new_arr.emplace_back(extend_allele);
            }
        }
        else if (a->equals(*StaticAllele::get_instance()->_span_del)) {
            if (!mm.old2new_map.count(a)) {
                mm.old2new_map.insert({a, a});
                mm.new_arr.emplace_back(a);
            }
        }
    }
}

void ROVACAVariantContextUtils::create_allele_mapping(pAllele ref_allele, pVariant vc, AlleleMap& mm)
{
    return create_allele_mapping(ref_allele, vc->ref_allele(), vc->alt_alleles(), mm);
}

void ROVACAVariantContextUtils::resolve_incompatible_alleles(pAllele ref_allele, pVariant vc, AlleleMap& mm)
{
    if (ref_allele->equals(*vc->ref_allele())) {
        const AlleleVector& alt_alleles = vc->alt_alleles();
        std::for_each(alt_alleles.begin(), alt_alleles.end(), [&](pAllele a) {
            if (!mm.old2new_map.count(a)) {
                mm.old2new_map.insert({a, a});
                mm.new_arr.emplace_back(a);
            }
        });
    }
    else {
        create_allele_mapping(ref_allele, vc, mm);
    }
}

pAllele ROVACAVariantContextUtils::determine_reference_allele(pAllele ref1, pAllele ref2)
{
    if (nullptr == ref1 || ref1->length() < ref2->length()) {
        return ref2;
    }
    else if (nullptr == ref2 || ref2->length() < ref1->length()) {
        return ref1;
    }
    else if (ref1->length() == ref2->length() && !ref1->equals(*ref2)) {
        RovacaLogger::error("the provided reference alleles do not appear to represent the same position");
        exit(EXIT_FAILURE);
    }
    else {
        return ref1;
    }
}

pAllele ROVACAVariantContextUtils::determine_reference_allele(const VariantVector& vcs, pInterfaceLocatable loc)
{
    pAllele ref{nullptr}, vc_ref{nullptr};
    for (pVariant vc : vcs) {
        if (nullptr == loc || loc->get_start() == vc->get_start()) {
            vc_ref = vc->ref_allele();
            ref = determine_reference_allele(ref, vc_ref);
        }
    }
    return ref;
}

pAllele ROVACAVariantContextUtils::determine_reference_allele(const VariantVector& vcs) { return determine_reference_allele(vcs, nullptr); }

VariantVector ROVACAVariantContextUtils::sort_variant_contexts_by_priority(const VariantVector& unsorted_vcs, GenotypeMergeType g)
{
    CHECK_CONDITION_EXIT(g == PRIORITIZE && unsorted_vcs.empty(), "cannot merge calls by priority with a null priority list");

    VariantVector result(unsorted_vcs.get_allocator());
    if (unsorted_vcs.empty() || g == UNSORTED) {
        result = unsorted_vcs;
    }
    else {
        result.reserve(unsorted_vcs.size());
        std::pmr::map<int32_t, pVariant> sorted;
        for (pVariant vc : unsorted_vcs) {
            sorted.insert({vc->source_id(), vc});
        }
        for (const auto& itr : sorted) {
            result.emplace_back(itr.second);
        }
    }
    return result;
}

pVariant ROVACAVariantContextUtils::simple_merge(const VariantVector& unsorted_vcs, FilteredRecordMergeType f, GenotypeMergeType g,
                                               bool filtered_are_uncalled)
{
    un_used(f);
    un_used(filtered_are_uncalled);

    if (unsorted_vcs.empty()) {
        return nullptr;
    }

    VariantVector pre_filtered_vcs = sort_variant_contexts_by_priority(unsorted_vcs, g);

    /*!
     * @note 此处原本有一个过滤操作
     * make sure all variant contexts are padded with reference base in case of indels if necessary
     */

    pVariant first = pre_filtered_vcs.front();
    int32_t source_id = first->source_id();
    pAllele ref_allele = determine_reference_allele(pre_filtered_vcs);

    pVariant longest_vc = first;

    pMemoryPool pool = unsorted_vcs.get_allocator().resource();

    size_t total_allele_num = 1;
    std::for_each(pre_filtered_vcs.begin(), pre_filtered_vcs.end(),
                  [&total_allele_num](pVariant vc) { total_allele_num += vc->allele_num() - 1; });
    AlleleVector new_alleles{pool};
    new_alleles.reserve(total_allele_num);
    new_alleles.push_back(ref_allele);

    AlleleMap sort_allele_map{pool};  // 既要无重复，又要不变更相对顺序，又要可以满足old_allele -> new_allele
    std::pmr::unordered_set<pAllele, AlleleHash, AlleleEqual> unique_set{pool};

    for (pVariant vc : pre_filtered_vcs) {
        CHECK_CONDITION_EXIT(vc->get_start() != longest_vc->get_start(), "bug: attempting to merge Variant with different start sites");
        if (vc->get_length_on_reference() > longest_vc->get_length_on_reference()) {
            longest_vc = vc;
        }

        resolve_incompatible_alleles(ref_allele, vc, sort_allele_map);
        std::for_each(sort_allele_map.new_arr.begin(), sort_allele_map.new_arr.end(), [&](pAllele a) {
            if (!unique_set.count(a)) {
                unique_set.insert(a);
                new_alleles.push_back(a);
            }
        });
        sort_allele_map.clear();

        /*!
         * @note 此处还有一些关于 Variant 属性的操作，HC未使用，暂不实现
         */
    }

    /*!
     * @note 此处还有一些关于 Variant 属性的操作，HC未使用，暂不实现
     */

    pVariant new_vc = Variant::create(pool);
    new_vc->set_source_id(source_id);
    new_vc->set_tid(longest_vc->get_tid());
    new_vc->set_start(longest_vc->get_start());
    new_vc->set_stop(longest_vc->get_stop());
    new_vc->set_alleles(new_alleles);
    return new_vc;
}

pGenotypesContext ROVACAVariantContextUtils::subset_to_ref_only(pVariant vc, int32_t default_ploidy, pMemoryPool pool)
{
    CHECK_CONDITION_EXIT(nullptr == vc, "vc cannot be nullptr");
    CHECK_CONDITION_EXIT(default_ploidy < 1, "default_ploidy must be >= 1");

    pAllele ref_allele = vc->ref_allele();
    AlleleVector diploid_ref_alleles({ref_allele, ref_allele}, pool);
    pGenotypesContext old_genotypes = vc->genotype();
    pGenotype g;
    int32_t ploidy;
    for (size_t gi = 0, g_len = old_genotypes->size(); gi < g_len; ++gi) {
        g = old_genotypes->at(gi);
        ploidy = g->get_ploidy() == 0 ? default_ploidy : g->get_ploidy();
        g->set_alleles({(size_t)ploidy, ref_allele, pool});
    }

    return old_genotypes;
}

void ROVACAVariantContextUtils::make_genotype_call(int32_t ploidy, pGenotype g, GenotypeAssignmentMethod method,
                                                 const DoubleVector& likelihoods, const AlleleVector& alleles_to_use,
                                                 const AlleleVector& original_gt, [[maybe_unused]] pGenotypePriorCalculator gpc)
{
    CHECK_CONDITION_EXIT(original_gt.empty() && method == BEST_MATCH_TO_ORIGINAL,
                         "original_gt cannot be empty if method is BEST_MATCH_TO_ORIGINAL");
    pMemoryPool pool = alleles_to_use.get_allocator().resource();

    switch (method) {
        case SET_TO_NO_CALL:
        case SET_TO_NO_CALL_NO_ANNOTATIONS: {
            g->set_alleles({(size_t)ploidy, StaticAllele::get_instance()->_no_call.get(), likelihoods.get_allocator()});
            break;
        }
        case USE_PLS_TO_ASSIGN:
        case PREFER_PLS: {
            if (likelihoods.empty() || !is_informative(likelihoods)) {
                if (method == PREFER_PLS) {
                    CHECK_CONDITION_EXIT(original_gt.empty(), "original_gt cannot be empty if method is BEST_MATCH_TO_ORIGINAL");
                    g->set_alleles(best_match_to_original_gt(alleles_to_use, original_gt));
                }
                else {
                    g->set_alleles({(size_t)ploidy, StaticAllele::get_instance()->_no_call.get(), likelihoods.get_allocator()});
                }
            }
            else {
                int32_t max_likelihood_index = (int32_t)math_utils::max_element_index(likelihoods);
                pGenotypeLikelihoodCalculator calculator =
                    GenotypeLikelihoodCalculator::create(ploidy, (int32_t)alleles_to_use.size(), pool);
                pGenotypeAlleleCounts gac = calculator->genotype_allele_counts_at(max_likelihood_index);

                AlleleVector final_alleles = gac->as_allele_list(alleles_to_use);
                auto find_itr = std::find_if(final_alleles.begin(), final_alleles.end(),
                                             [](pAllele a) { return a->equals(*StaticAllele::get_instance()->_non_ref_allele); });
                if (find_itr != std::end(final_alleles)) {
                    pAllele ref = *std::find_if(alleles_to_use.begin(), alleles_to_use.end(), [](pAllele a) { return a->is_reference(); });
                    AlleleVector new_alleles(ploidy, ref, alleles_to_use.get_allocator());
                    Int32Vector new_pl(likelihoods.size(), alleles_to_use.get_allocator());
                    g->set_alleles(std::move(new_alleles));
                    g->set_pl(std::move(new_pl));
                    g->set_log10perror(0.0);
                }
                else {
                    g->set_alleles(std::move(final_alleles));
                }
                int32_t num_alt_alleles = (int32_t)alleles_to_use.size() - 1;
                if (num_alt_alleles > 0) {
                    g->set_log10perror(GenotypeLikelihoods::get_gq_log10from_likelihoods(max_likelihood_index, likelihoods, pool));
                }
            }
            break;
        }
        case BEST_MATCH_TO_ORIGINAL: {
            g->set_alleles(best_match_to_original_gt(alleles_to_use, original_gt));
            break;
        }
        case USE_POSTERIORS_ANNOTATION:
        case DO_NOT_ASSIGN_GENOTYPES:
        case USE_POSTERIOR_PROBABILITIES:
        default: {
            // todo: hc未使用，暂未实现
            RovacaLogger::error("invalid code");
            exit(EXIT_FAILURE);
        }
    }
}

bool ROVACAVariantContextUtils::is_informative(const DoubleVector& gls)
{
    return std::accumulate(gls.begin(), gls.end(), 0.0) < s_sum_gl_thresh_nocall;
}

AlleleVector ROVACAVariantContextUtils::best_match_to_original_gt(const AlleleVector& alleles_to_use, const AlleleVector& original_gt)
{
    AlleleVector best(alleles_to_use.get_allocator());
    best.reserve(original_gt.size());
    pAllele ref = alleles_to_use.front();
    for (pAllele original_allele : original_gt) {
        auto find_itr = std::find(alleles_to_use.begin(), alleles_to_use.end(), original_allele);
        best.emplace_back((find_itr != std::end(alleles_to_use) || !original_allele->is_called()) ? original_allele : ref);
    }
    return best;
}

pVariant ROVACAVariantContextUtils::reverse_trim_alleles(pVariant untrimmed_result) { return trim_alleles(untrimmed_result, false, true); }

pVariant ROVACAVariantContextUtils::trim_alleles(pVariant input_vc, bool trim_forward, bool trim_reverse)
{
    const AlleleVector& alleles = input_vc->alleles();
    if (alleles.size() <= 1 || std::any_of(alleles.begin(), alleles.end(), [](pAllele a) {
            return a->length() == 1 && !a->equals(*StaticAllele::get_instance()->_span_del);
        })) {
        return input_vc;
    }

    pMemoryPool pool = input_vc->alleles().get_allocator().resource();
    std::pmr::vector<uint8_t*> sequences{pool};
    sequences.reserve(alleles.size());
    std::pmr::vector<IndexRange> bounds{pool};
    bounds.reserve(alleles.size());

    std::for_each(alleles.begin(), alleles.end(), [&](pAllele a) {
        if (!a->is_symbolic() && !a->equals(*StaticAllele::get_instance()->_span_del)) {
            sequences.push_back(a->get_bases()->data);
            bounds.emplace_back(0, (int32_t)a->length());
        }
    });

    std::pair<int32_t, int32_t> shifts = AlignmentUtils::normalize_alleles(sequences, bounds, 0, true);
    int32_t end_trim = shifts.second;
    int32_t start_trim = -shifts.first;

    bool empty_allele = std::any_of(bounds.begin(), bounds.end(), [](const IndexRange& i) { return i.size() == 0; });
    bool restore_one_base_at_end = empty_allele && start_trim == 0;
    bool restore_one_base_at_start = empty_allele && start_trim > 0;

    // if the end trimming consumed all the bases, leave one base
    int32_t end_bases_to_clip = restore_one_base_at_end ? end_trim - 1 : end_trim;
    int32_t start_bases_to_clip = restore_one_base_at_start ? start_trim - 1 : start_trim;

    return trim_alleles(input_vc, (trim_forward ? start_bases_to_clip : 0) - 1, trim_reverse ? end_bases_to_clip : 0);
}

pVariant ROVACAVariantContextUtils::trim_alleles(pVariant input_vc, int32_t fwd_trim_end, int32_t rev_trim)
{
    // nothing to do, so just return input_vc unmodified
    if (fwd_trim_end == -1 && rev_trim == 0) {
        return input_vc;
    }

    const AlleleVector& alleles = input_vc->alleles();
    pMemoryPool pool = alleles.get_allocator().resource();

    AlleleMap mm{pool};

    pAllele span_del = StaticAllele ::get_instance()->_span_del.get();
    for (const auto& item : alleles) {
        if (item->is_symbolic() || item->equals(*span_del)) {
            mm.new_arr.push_back(item);
            mm.old2new_map.insert({item, item});
        }
        else {
            uint32_t new_length = item->length() - rev_trim - fwd_trim_end - 1;
            pBases old_bases = item->get_display_string();
            pAllele new_allele;
            if (1 == new_length) {
                new_allele = Allele::create_allele(old_bases->data[fwd_trim_end + 1], item->is_reference());
            }
            else {
                auto new_bases = new ALLOC_FLEXIBLE_IN_POOL(pool, Bases, new_length, uint8_t) Bases{new_length};
                memcpy(new_bases->data, old_bases->data + fwd_trim_end + 1, new_length * sizeof(uint8_t));
                new_allele = Allele::create_allele(new_bases, item->is_reference(), pool);
            }
            mm.new_arr.push_back(new_allele);
            mm.old2new_map.insert({item, new_allele});
        }
    }

    pGenotypesContext genotypes = update_genotypes_with_mapped_alleles(input_vc->genotype(), mm);

    int64_t start = input_vc->get_start() + fwd_trim_end + 1;
    input_vc->set_start(INVALID_INT);
    input_vc->set_stop(start + mm.new_arr.at(0)->length() - 1);
    input_vc->set_alleles(mm.new_arr);
    input_vc->set_start(start);  // start、stop、alleles都被设置时会触发校验
    input_vc->set_genotype(genotypes);
    return input_vc;
}

pGenotypesContext ROVACAVariantContextUtils::update_genotypes_with_mapped_alleles(pGenotypesContext genotypes, const AlleleMap& mm)
{
    if (nullptr == genotypes) {
        return nullptr;
    }

    pGenotype g;
    for (size_t i = 0, len = genotypes->size(); i < len; ++i) {
        g = genotypes->at(i);
        g->set_alleles(mm.remap(g->alleles()));
    }
    return genotypes;
}

bool ROVACAVariantContextUtils::get_num_tandem_repeat_units(pVariant vc, pRefFragment ref_bases_starting_at_vc_without_pad,
                                                          pRepeatUnitsResult result, pMemoryPool pool)
{
    CHECK_CONDITION_EXIT(vc == nullptr || ref_bases_starting_at_vc_without_pad == nullptr, "nullptr");
    if (!vc->is_indel1()) {
        return false;
    }

    RefFragment ref_allele_bases;
    ref_allele_bases.len = vc->ref_allele()->length() - 1;
    ref_allele_bases.data = &(vc->ref_allele()->get_bases()->data[1]);

    RepeatUnitsResult loop_result{pool};
    RefFragment alt_allele_bases;
    const AlleleVector& alt_alleles = vc->alt_alleles();
    for (pAllele allele : alt_alleles) {
        alt_allele_bases.len = allele->length() - 1;
        alt_allele_bases.data = &(allele->get_bases()->data[1]);
        get_num_tandem_repeat_units(&ref_allele_bases, &alt_allele_bases, ref_bases_starting_at_vc_without_pad, &loop_result, pool);

        const Int32Vector& repetition_count = loop_result.lengths;
        // repetition count = 0 means allele is not a tandem expansion of context
        if (repetition_count.at(0) == 0 || repetition_count.at(1) == 0) {
            return false;
        }

        if (result->lengths.empty()) {
            result->lengths.push_back(repetition_count.at(0));
        }
        result->lengths.push_back(repetition_count.at(1));
        result->bases = loop_result.bases;
    }
    return true;
}

bool ROVACAVariantContextUtils::get_num_tandem_repeat_units(pRefFragment ref_bases, pRefFragment alt_bases,
                                                          pRefFragment remaining_ref_context, pRepeatUnitsResult result, pMemoryPool pool)
{
    // find first repeat unit based on either ref or alt, whichever is longer
    pRefFragment long_b = alt_bases->len > ref_bases->len ? alt_bases : ref_bases;

    result->bases.len = find_repeated_substring(long_b);
    result->bases.data = long_b->data;

    if (result->lengths.size() < 2) {
        result->lengths.resize(2);
    }

    // look for repetitions forward on the ref bases (i.e. starting at beginning of ref bases)
    int32_t repetitions_in_ref = find_number_of_repetitions(&result->bases, ref_bases, true);

    uint32_t max_len = (alt_bases->len > ref_bases->len ? alt_bases->len : ref_bases->len) + remaining_ref_context->len;
    RefFragment test_string{0, nullptr};
    test_string.data = (uint8_t*)pool->allocate(max_len);

    memset(test_string.data, 0, max_len);
    test_string.len = ref_bases->len + remaining_ref_context->len;
    memcpy(test_string.data, ref_bases->data, ref_bases->len * sizeof(uint8_t));
    memcpy(test_string.data + ref_bases->len, remaining_ref_context->data, remaining_ref_context->len * sizeof(uint8_t));
    result->lengths[0] = find_number_of_repetitions(&result->bases, &test_string, true) - repetitions_in_ref;

    memset(test_string.data, 0, max_len);
    test_string.len = alt_bases->len + remaining_ref_context->len;
    memcpy(test_string.data, alt_bases->data, alt_bases->len * sizeof(uint8_t));
    memcpy(test_string.data + alt_bases->len, remaining_ref_context->data, remaining_ref_context->len * sizeof(uint8_t));
    result->lengths[1] = find_number_of_repetitions(&result->bases, &test_string, true) - repetitions_in_ref;

    return true;
}

uint32_t ROVACAVariantContextUtils::find_repeated_substring(pRefFragment bases)
{
    uint32_t rep_length;
    bool all_bases_match;
    RefFragment candidate_repeat_unit, base_piece;
    for (rep_length = 1; rep_length <= bases->len; ++rep_length) {
        candidate_repeat_unit.data = bases->data;
        candidate_repeat_unit.len = rep_length;
        all_bases_match = true;
        for (uint32_t start = rep_length; start < bases->len; start += rep_length) {
            base_piece.len = candidate_repeat_unit.len;
            base_piece.data = bases->data + start;
            if (!ref_fragment_equal(&candidate_repeat_unit, &base_piece)) {
                all_bases_match = false;
                break;
            }
        }
        if (all_bases_match) {
            return rep_length;
        }
    }
    return rep_length;
}

int32_t ROVACAVariantContextUtils::find_number_of_repetitions(pRefFragment repeat_unit, pRefFragment test_string, bool leading_repeats)
{
    CHECK_CONDITION_EXIT(repeat_unit == nullptr, "repeat_unit");
    CHECK_CONDITION_EXIT(test_string == nullptr, "test_string");
    CHECK_CONDITION_EXIT(repeat_unit->len == 0, "repeat_unit->len == 0");
    if (test_string->len == 0) {
        return 0;
    }
    return find_number_of_repetitions(repeat_unit, 0, int32_t(repeat_unit->len), test_string, 0, int32_t(test_string->len),
                                      leading_repeats);
}

int32_t ROVACAVariantContextUtils::find_number_of_repetitions(pRefFragment repeat_unit_full, int32_t offset_in_repeat_unit_full,
                                                            int32_t repeat_unit_length, pRefFragment test_string_full,
                                                            int32_t offset_in_test_string_full, int32_t test_string_length,
                                                            bool leading_repeats)
{
    if (test_string_length == 0) {
        return 0;
    }
    int32_t length_difference = test_string_length - repeat_unit_length;
    if (leading_repeats) {
        int32_t num_repeats = 0;
        // look forward on the test string
        for (int32_t start = 0; start <= length_difference; start += repeat_unit_length) {
            if (equal_range(test_string_full, start + offset_in_test_string_full, repeat_unit_full, offset_in_repeat_unit_full,
                            repeat_unit_length)) {
                num_repeats++;
            }
            else {
                return num_repeats;
            }
        }
        return num_repeats;
    }
    else {
        // look backward. for example, if repeat_unit = at and test_string = gatat, number of repeat units is still 2
        int32_t num_repeats = 0;
        // look backward on the test string
        for (int32_t start = length_difference; start >= 0; start -= repeat_unit_length) {
            if (equal_range(test_string_full, start + offset_in_test_string_full, repeat_unit_full, offset_in_repeat_unit_full,
                            repeat_unit_length)) {
                num_repeats++;
            }
            else {
                return num_repeats;
            }
        }
        return num_repeats;
    }
}

bool ROVACAVariantContextUtils::ref_fragment_equal(pRefFragment first, pRefFragment second)
{
    CHECK_CONDITION_EXIT(first == nullptr || second == nullptr, "nullptr");
    if (first->len != second->len) {
        return false;
    }
    return std::memcmp(first->data, second->data, first->len) == 0;
}

bool ROVACAVariantContextUtils::equal_range(pRefFragment left, int32_t left_offset, pRefFragment right, int32_t right_offset, int32_t length)
{
    CHECK_CONDITION_EXIT(nullptr == left || nullptr == right, "nullptr");
    CHECK_CONDITION_EXIT(left_offset + length > int32_t(left->len), "left_offset + length > int32_t(left->len)");
    CHECK_CONDITION_EXIT(right_offset + length > int32_t(right->len), "right_offset + length > int32_t(right->len)");

    uint8_t* left_data = left->data + left_offset;
    uint8_t* right_data = right->data + right_offset;
    return std::memcmp(left_data, right_data, length) == 0;
}

pVariant ROVACAVariantContextUtils::get_overlapping_variant_context(InterfaceLocatable* cur_pos, const VariantVector& maybe_overlapping)
{
    pVariant ret = nullptr;
    for (pVariant vc : maybe_overlapping) {
        if (cur_pos->overlaps(*vc)) {
            if (nullptr == ret || vc->get_start() > ret->get_start()) {
                ret = vc;
            }
        }
    }
    return ret;
}

AlleleVector ROVACAVariantContextUtils::homozygous_allele_list(pAllele a, int32_t ploidy, pMemoryPool pool)
{
    CHECK_CONDITION_EXIT(nullptr == a, "nullptr == a");
    return {(size_t)ploidy, a, pool};
}

int32_t ROVACAVariantContextUtils::calculate_gq_from_pls(const Int32Vector& pls)
{
    CHECK_CONDITION_EXIT(pls.size() < 2, "array of pl values must contain at least two elements");
    int32_t first = pls[0];
    int32_t second = pls[1];
    if (first > second) {
        second = first;
        first = pls[1];
    }
    for (int32_t i = 2, len = int32_t(pls.size()); i < len; i++) {
        int32_t candidate = pls[i];
        if (candidate >= second) {
            continue;
        }
        if (candidate <= first) {
            second = first;
            first = candidate;
        }
        else {
            second = candidate;
        }
    }
    return second - first;
}

int32_t ROVACAVariantContextUtils::calculate_gq_from_pls(const std::vector<int32_t>& pls)
{
    CHECK_CONDITION_EXIT(pls.size() < 2, "array of pl values must contain at least two elements");
    int32_t first = pls[0];
    int32_t second = pls[1];
    if (first > second) {
        second = first;
        first = pls[1];
    }
    for (int32_t i = 2, len = int32_t(pls.size()); i < len; i++) {
        int32_t candidate = pls[i];
        if (candidate >= second) {
            continue;
        }
        if (candidate <= first) {
            second = first;
            first = candidate;
        }
        else {
            second = candidate;
        }
    }
    return second - first;
}

std::pmr::vector<Event> ROVACAVariantContextUtils::split_variant_context_to_biallelics_event(pVariant vc, bool trim_left, pMemoryPool pool)
{
    VariantVector vcs = ROVACAVariantContextUtils::split_variant_context_to_biallelics(vc, trim_left, pool);
    std::pmr::vector<Event> result{pool};
    result.reserve(vcs.size());

    std::for_each(vcs.begin(), vcs.end(), [&](pVariant v) {
        result.push_back({v->get_tid(), v->get_start(), v->ref_allele(), v->biallelic_alt(), pool});
    });

    return result;
}

VariantVector ROVACAVariantContextUtils::split_variant_context_to_biallelics(pVariant vc, bool trim_left, pMemoryPool pool)
{
    VariantVector result{pool};
    if (!vc->is_variant()) {
        return result;
    }
    else if (vc->is_biallelic()) {
        result.push_back(vc);
    }
    else {
        pVariant new_vc = nullptr;
        int32_t tid = vc->get_tid();
        int64_t start = vc->get_start();
        int64_t stop = vc->get_stop();
        pAllele ref = vc->ref_allele();
        pBases rsid = vc->db_id();
        const AlleleVector& alt_alleles = vc->alt_alleles();
        for (pAllele alt : alt_alleles) {
            AlleleVector new_alleles{{ref, alt}, pool};
            new_vc = Variant::create(pool);

            new_vc->set_tid(tid);
            new_vc->set_start(start);
            new_vc->set_stop(stop);
            new_vc->set_alleles(new_alleles);
            new_vc->set_id(rsid);

            result.push_back(trim_alleles(new_vc, trim_left, true));
        }
    }

    return result;
}

static bool different_last_base(pBases ref, pBases alt)
{
    return 0 == ref->num || 0 == alt->num || ref->data[ref->num - 1] != alt->data[alt->num - 1];
}

void Event::make_minimal_representation(pAllele ref, pAllele alt, pMemoryPool pool)
{
    if (1 == ref->length() || 1 == alt->length() || different_last_base(ref->get_bases(), alt->get_bases())) {
        ref_ = ref;
        alt_ = alt;
    }
    else {
        pBases ref_base = ref->get_bases();
        pBases alt_base = alt->get_bases();
        if (ref_base->num == alt_base->num && 0 == strncmp((char*)ref_base->data, (char*)alt_base->data, ref_base->num)) {
            RovacaLogger::error("ref and alt alleles are identical");
            exit(0);
        }

        uint32_t overlap_count = 0;
        uint32_t min_len = std::min(ref_base->num, alt_base->num);
        while (overlap_count < min_len &&
               ref_base->data[ref_base->num - 1 - overlap_count] == alt_base->data[alt_base->num - 1 - overlap_count]) {
            overlap_count++;
        }

        pBases nref_base, nalt_base;
        nref_base = new ALLOC_FLEXIBLE_IN_POOL(pool, Bases, ref_base->num - overlap_count, uint8_t) Bases{ref_base->num - overlap_count};
        nref_base->data[nref_base->num] = 0;
        memcpy((char*)nref_base->data, (char*)ref_base->data, nref_base->num);
        nalt_base = new ALLOC_FLEXIBLE_IN_POOL(pool, Bases, alt_base->num - overlap_count, uint8_t) Bases{alt_base->num - overlap_count};
        nalt_base->data[nalt_base->num] = 0;
        memcpy((char*)nalt_base->data, (char*)alt_base->data, nalt_base->num);

        ref_ = Allele::create_allele(nref_base, true, pool);
        alt_ = Allele::create_allele(nalt_base, false, pool);
    }

    stop_ = start_ + ref_->length() - 1;
}

bool Event::operator==(const Event& o) const { return start_ == o.start_ && ref_->equals(*o.ref_) && alt_->equals(*o.alt_); }

}  // namespace rovaca