#include "germline_genotying_engine.h"

#include "af_calculation_result.h"
#include "allele.h"
#include "allele_frequency_calculator.h"
#include "allele_likelihoods.hpp"
#include "block_combiner.h"
#include "event_map.h"
#include "forward.h"
#include "genotype.h"
#include "genotype_argument.h"
#include "genotype_likelihoods.h"
#include "genotype_macors.h"
#include "genotype_prior_calculator.h"
#include "genotypes_context.hpp"
#include "homogeneous_ploidy_model.hpp"
#include "independent_sample_genotypes_model.h"
#include "indexed_allele_list.hpp"
#include "indexed_sample_list.hpp"
#include "info_data.hpp"
#include "rovaca_logger.h"
#include "read_record.h"
#include "reference_confidence_model.h"
#include "simple_interval.h"
#include "utils/adapter_utils.h"
#include "utils/alignment_utils.h"
#include "utils/allele_subsetting_utils.h"
#include "utils/assembly_based_caller_utils.h"
#include "utils/rovaca_variant_context_utils.h"
#include "utils/read_record_utils.h"
#include "variant.h"
#include "variant_annotator_engine.h"

#define PRACTICAL_ALLELE_COUNT_TWO_PLOIDY (44)

namespace rovaca
{

static constexpr int32_t s_too_long_pl = 100000;

struct OutputAlleleSubset
{
    AlleleVector alleles;
    Int32Vector mle_counts;
    bool site_is_monomorphic;
    explicit OutputAlleleSubset(pMemoryPool pool)
        : alleles(pool)
        , mle_counts(pool)
        , site_is_monomorphic(true)
    {
        alleles.reserve(5);
        mle_counts.reserve(5);
    }
};

GermlineGenotyingEngine::GermlineGenotyingEngine()
    : sw_(sw_avx_init(false))
    , _upstream_deletions_loc()
    , _deletions_locs()
{}

GermlineGenotyingEngine::~GermlineGenotyingEngine()
{
    sw_avx_finit(sw_);
    clear_upstream_deletions_loc();
    delete _block_combiner;
    delete _annotator_engine;
    delete _prior_calculator;
    delete _genotyping_model;
    delete _allele_frequency_calculator;
}

void GermlineGenotyingEngine::init_engine_per_loop(pHCArgs args, pMemoryPool pool, pBamDataPool bpool, bam_hdr_t* header,
                                                   bcf_hdr_t* bcf_hdr, pInterfaceSampleList samples, pInterfacePloidyModel pm)
{
    _args = args;
    _pool = pool;
    _bam_pool = bpool;
    _header = header;
    _sample_list = samples;
    _ploidy_model = pm;
    _upstream_deletions_loc.clear();

    if (ROVACA_UNLIKELY(!_block_combiner)) {
        _block_combiner = new BlockCombiner{args->gvcf_gq_bands, bcf_hdr, header, false,
                                            args->reference_confidence_mode == ReferenceConfidenceMode::GVCF};
    }
    if (ROVACA_UNLIKELY(!_annotator_engine)) {
        init_engine();
    }
}

CalledHaplotypes GermlineGenotyingEngine::assign_genotype_likelihoods(pRHLikelihoods rh_likelihoods, pRefFragment ref,
                                                                      pSimpleInterval ref_loc, pSimpleInterval active_region,
                                                                      const Int32ToReadVectorMap& per_sample_filtered_read)
{
    CHECK_CONDITION_EXIT(ref_loc->size() != int64_t(ref->len), "ref_loc length must match ref bytes");
    CHECK_CONDITION_EXIT(!ref_loc->contains(*active_region), "ref_loc must contain active_region");

    auto& haplotypes = const_cast<HaplotypeVector&>(rh_likelihoods->get_alleles());
    Int64Set start_pos_key_set = EventMap::build_event_maps_for_haplotypes(haplotypes, ref, ref_loc, _args->max_mnp_distance, _pool);

    VariantVector return_calls{_pool};
    return_calls.reserve(20);
    AlleleVector given_alleles{_pool};
    HaplotypeList called_haplotypes{_pool};

    size_t merged_alleles_list_size_before_possible_trimming;

    int32_t ploidy = _args->sample_ploidy;
    bool esd = _args->disable_spanning_event_genotyping;
    bool mnd = _args->max_mnp_distance;
    bool erc = emit_reference_confidence();
    int64_t irom = _args->informative_read_overlap_margin;
    int64_t contig_length;
    pAllele ref_allele;
    pSimpleInterval variant_calling_relevant_overlap;

    for (int64_t loc : start_pos_key_set) {
        if (loc < active_region->get_start() || loc > active_region->get_stop()) {
            continue;
        }

        VariantVector events_at_this_loc = AssemblyBasedCallerUtils::get_variant_contexts_from_active_haplotypes(loc, haplotypes, !mnd);

        ref_allele = Allele::create_allele(ref->data[loc - ref_loc->get_start()], true);
        VariantVector events_span_dels_replaced = replace_span_dels(events_at_this_loc, ref_allele, loc);

        pVariant merged_vc = AssemblyBasedCallerUtils::make_merged_variant_context(events_span_dels_replaced);
        if (nullptr == merged_vc) {
            continue;
        }

        merged_alleles_list_size_before_possible_trimming = merged_vc->allele_num();
        AlleleToHaplotypeVectorMap allele_mapper = AssemblyBasedCallerUtils::create_allele_mapper(haplotypes, merged_vc, loc, !esd);

        if (ROVACA_UNLIKELY(allele_mapper.size() > PRACTICAL_ALLELE_COUNT_TWO_PLOIDY)) {
            merged_vc = remove_alt_alleles_if_too_many_genotypes(ploidy, allele_mapper, merged_vc);
            if (nullptr == merged_vc) continue;
        }

        const AlleleVector& alleles = merged_vc->alleles();
        pRALikelihoods ra_likelihoods = rh_likelihoods->marginalize(alleles, allele_mapper);

        contig_length = int64_t(_header->target_len[merged_vc->get_tid()]);
        variant_calling_relevant_overlap = SimpleInterval::create(*merged_vc, _pool)->expand_within_contig(irom, contig_length, _pool);
        ra_likelihoods->retain_evidence(variant_calling_relevant_overlap);
        ra_likelihoods->set_variant_calling_subset_used(variant_calling_relevant_overlap);

        if (erc) {
            merged_vc->add_non_ref_symbolic_allele();
            ra_likelihoods->add_non_reference_allele();
            ++merged_alleles_list_size_before_possible_trimming;
        }

        pGenotypesContext genotypes = calculate_gls_for_this_event(ra_likelihoods, merged_vc, ref, loc - ref_loc->get_start());
        merged_vc->set_genotype(genotypes);
        pVariant call = calculate_genotypes(merged_vc, given_alleles);
        if (nullptr != call) {
            ra_likelihoods = prepare_read_allele_likelihoods_for_annotation(rh_likelihoods, per_sample_filtered_read, erc, allele_mapper,
                                                                            ra_likelihoods, call, variant_calling_relevant_overlap);

            std::for_each(call->alleles().begin(), call->alleles().end(), [&](pAllele a) {
                if (allele_mapper.count(a)) {
                    const HaplotypeVector& hs = allele_mapper.at(a);
                    called_haplotypes.insert(called_haplotypes.end(), hs.begin(), hs.end());
                }
            });

            pVariant annotated_call = make_annotated_call(ref, call, ra_likelihoods, merged_alleles_list_size_before_possible_trimming);

            return_calls.push_back(annotated_call);
        }
    }

    if (!_args->do_not_run_physical_phasing) {
        VariantVector phasing_vc = AssemblyBasedCallerUtils::phase_calls(return_calls, called_haplotypes);
        return {phasing_vc, called_haplotypes};
    }
    else {
        return {return_calls, {}};
    }
}

VariantVector GermlineGenotyingEngine::call_non_active_site(pHaplotype ref_h, pRefFragment ref, pSimpleInterval ref_loc,
                                                            pSimpleInterval original, pSimpleInterval original_padded,
                                                            pSimpleInterval variant, pSimpleInterval variant_padded,
                                                            const CalledHaplotypes& calls, const ReadHashSet& original_reads,
                                                            const ReadHashSet& genotype_reads, std::pmr::list<bam1_t*>& extra_memory_reads)
{
    if (!ROVACA_LIKELY(emit_reference_confidence())) {
        return {calls.first, _pool};
    }

    int32_t ploidy = _ploidy_model->sample_ploidy(0);
    int64_t padding = _args->assembly_region_padding;
    int64_t chr_len = _header->target_len[original->get_tid()];

    VariantVector left_result{}, middle_result{}, right_result{};

    // 计算左侧
    IntervalPair lrp = AdapterUtils::non_variant_left_flank_region(original, original_padded, variant, padding, chr_len, _pool);
    if (lrp.first) {
        ReadHashSet left_overlaps_reads =
            AdapterUtils::trim_reads_by_region(original_reads, lrp.second, _pool, _bam_pool, extra_memory_reads);
        left_result = reference_model_for_no_variation(ref, ref_loc, lrp.first, lrp.second, ploidy, left_overlaps_reads);
    }

    // 计算中间
    ReferenceConfidenceModel r{_pool};
    middle_result = r.calculate_ref_confidence(calls.first, ref_h, variant, variant_padded, ploidy, genotype_reads);

    // 计算右侧
    IntervalPair rrp = AdapterUtils::non_variant_right_flank_region(original, original_padded, variant, padding, chr_len, _pool);
    if (rrp.first) {
        ReadHashSet right_overlaps_reads =
            AdapterUtils::trim_reads_by_region(original_reads, rrp.second, _pool, _bam_pool, extra_memory_reads);
        right_result = reference_model_for_no_variation(ref, ref_loc, rrp.first, rrp.second, ploidy, right_overlaps_reads);
    }

    VariantVector result{_pool};
    result.reserve(left_result.size() + middle_result.size() + right_result.size());
    std::copy(left_result.begin(), left_result.end(), std::back_inserter(result));
    std::copy(middle_result.begin(), middle_result.end(), std::back_inserter(result));
    std::copy(right_result.begin(), right_result.end(), std::back_inserter(result));
    return result;
}

VariantVector GermlineGenotyingEngine::reference_model_for_no_variation(pRefFragment ref, pSimpleInterval ref_loc, pSimpleInterval active,
                                                                        pSimpleInterval padded, int32_t ploidy, ReadHashSet& reads)
{
    AdapterUtils::filter_non_passing_reads2(reads, _args->mapping_quality_threshold, _pool);
    pHaplotype ref_haplotype = AssemblyBasedCallerUtils::create_reference_haplotype(ref, ref_loc, padded, _pool);
    return ReferenceConfidenceModel{_pool}.calculate_ref_confidence({}, ref_haplotype, active, padded, ploidy, reads);
}

VariantVector GermlineGenotyingEngine::replace_span_dels(const VariantVector& events_at_this_loc, pAllele ref_allele, int64_t loc)
{
    VariantVector result(_pool);
    result.reserve(events_at_this_loc.size());
    for (pVariant vc : events_at_this_loc) {
        if (ROVACA_LIKELY(vc->get_start() == loc)) {
            result.emplace_back(vc);
        }
        else {
            AlleleVector new_alleles{ref_allele, StaticAllele::get_instance()->_span_del.get()};
            pVariant replace_vc = Variant::create(_pool);
            replace_vc->set_tid(vc->get_tid());
            replace_vc->set_start(loc);
            replace_vc->set_stop(loc);
            replace_vc->set_alleles(new_alleles);
            replace_vc->set_source_id(vc->source_id());
            result.emplace_back(replace_vc);
        }
    }

    return result;
}

pVariant GermlineGenotyingEngine::calculate_genotypes(pVariant vc, const AlleleVector& given_alleles)
{
    if (cannot_be_genotyped(vc) || vc->sample_count() == 0) {
        return nullptr;
    }

    int32_t default_ploidy = _args->sample_ploidy;
    int32_t max_alt_alleles = _args->max_alternate_alleles;

    pVariant reduced_vc = vc;
    if (max_alt_alleles < (int32_t)reduced_vc->alt_alleles().size()) {
        AlleleVector alleles_to_keep =
            AlleleSubsettingUtils::calculate_most_likely_alleles(vc, default_ploidy, max_alt_alleles, false, _pool);
        pGenotypesContext reduced_genotypes =
            alleles_to_keep.size() == 1
                ? ROVACAVariantContextUtils::subset_to_ref_only(vc, default_ploidy, _pool)
                : AlleleSubsettingUtils::subset_alleles(vc->genotype(), default_ploidy, vc->alleles(), alleles_to_keep, _prior_calculator,
                                                        BEST_MATCH_TO_ORIGINAL, _pool);
        reduced_vc = Variant::create(_pool);
        reduced_vc->set_tid(vc->get_tid());
        reduced_vc->set_start(vc->get_start());
        reduced_vc->set_stop(vc->get_stop());
        reduced_vc->set_alleles(alleles_to_keep);
        reduced_vc->set_source_id(vc->source_id());
        reduced_vc->set_genotype(reduced_genotypes);
    }

    int32_t allele_count = (int32_t)reduced_vc->allele_num();
    int32_t ploidy = reduced_vc->get_max_ploidy(default_ploidy);
    int32_t max_pllength = GenotypeLikelihoods::num_likelihoods(allele_count, ploidy);
    if (ROVACA_UNLIKELY(max_pllength >= s_too_long_pl)) {
        RovacaLogger::warn("length of pl arrays is likely to reach {}", max_pllength);
    }

    pAFCalculationResult af_result = _allele_frequency_calculator->calculate(reduced_vc, default_ploidy);
    AlleleSet forced_alleles = AssemblyBasedCallerUtils::get_alleles_consistent_with_given_alleles(given_alleles, reduced_vc);
    pOutputAlleleSubset output_alternative_alleles = calculate_output_allele_subset(af_result, reduced_vc, forced_alleles);

    // note the math.abs is necessary because -10 * 0.0 => -0.0 which isn't nice
    double log10confidence = !output_alternative_alleles->site_is_monomorphic || _args->annotate_all_sites_with_pls
                                 ? af_result->log10prob_only_ref_allele_exists() + 0.00
                                 : af_result->log10prob_variant_present() + 0.00;

    // add 0.0 removes -0.0 occurrences.
    double phred_scaled_confidence = (-10.0 * log10confidence) + 0.0;

    // return a null call if we don't pass the confidence cutoff or the most likely allele frequency is zero
    // skip this if we are already looking at a vc with NON_REF as the first alt allele i.e. if we are in GenotypeGVCFs
    if (!passes_emit_threshold(phred_scaled_confidence, output_alternative_alleles->site_is_monomorphic) && !emit_all_active_sites() &&
        no_alleles_or_first_allele_is_not_non_ref(output_alternative_alleles->alleles) && forced_alleles.empty()) {
        return nullptr;
    }

    // return a null call if we aren't forcing site emission and the only alt allele is a spanning deletion
    if (!emit_all_active_sites() && output_alternative_alleles->alleles.size() == 1 &&
        output_alternative_alleles->alleles.front()->equals(*StaticAllele::get_instance()->_span_del)) {
        return nullptr;
    }

    // start constructing the resulting VC
    AlleleVector output_alleles(_pool);
    output_alleles.reserve(1 + output_alternative_alleles->alleles.size());
    output_alleles.emplace_back(reduced_vc->ref_allele());
    std::copy(output_alternative_alleles->alleles.begin(), output_alternative_alleles->alleles.end(), std::back_inserter(output_alleles));
    record_deletions(reduced_vc, output_alleles);

    pVariant new_vc = Variant::create(_pool);
    new_vc->set_tid(reduced_vc->get_tid());
    new_vc->set_start(reduced_vc->get_start());
    new_vc->set_stop(reduced_vc->get_stop());
    new_vc->set_alleles(output_alleles);
    new_vc->set_log10_error(log10confidence);

    /*! @note 以下代码 hc 未使用，暂未实现 */
#if 0
    if (!passes_call_threshold(phred_scaled_confidence)) {
        new_vc->filter(GATKVCFConstants.LOW_QUAL_FILTER_NAME);
    }
#endif

    pGenotypesContext genotypes =
        1 == output_alleles.size()
            ? ROVACAVariantContextUtils::subset_to_ref_only(reduced_vc, default_ploidy, _pool)
            : AlleleSubsettingUtils::subset_alleles(reduced_vc->genotype(), default_ploidy, reduced_vc->alleles(), output_alleles,
                                                    _prior_calculator, _args->genotype_assignment_method, _pool);

    /*! @note 以下代码 hc 未使用，暂未实现 */
#if 0
    if (configuration.genotype_args.use_posterior_probabilities_to_calculate_qual && has_posteriors(genotypes)) {
        final double log10no_variant_posterior = phred_no_variant_posterior_probability(output_alleles, genotypes) * -.1;
        final double qual_update = !output_alternative_alleles.site_is_monomorphic || configuration.annotate_all_sites_with_pls
                                       ? log10no_variant_posterior + 0.0
                                       : math_utils.log10one_minus_pow10(log10no_variant_posterior) + 0.0;
        if (!double.is_na_n(qual_update)) {
            builder.log10perror(qual_update);
        }
    }
#endif

    pInfoData info = InfoData::create(_pool);
    compose_call_attributes(reduced_vc, output_alternative_alleles->mle_counts, af_result, output_alleles, genotypes, info);

    new_vc->set_genotype(genotypes);
    new_vc->set_info(info);

    return new_vc;
}

pGenotypesContext GermlineGenotyingEngine::calculate_gls_for_this_event(pRALikelihoods likelihoods, pVariant merged_vc, pRefFragment ref,
                                                                        int64_t offset_for_ref_into_event)
{
    CHECK_CONDITION_EXIT(nullptr == likelihoods, "nullptr == likelihoods");
    CHECK_CONDITION_EXIT(nullptr == merged_vc, "nullptr == merged_vc");

    const AlleleVector& variant_alleles = merged_vc->alleles();
    InterfaceAlleleList<pAllele>* allele_list;
    if (ROVACA_LIKELY(likelihoods->number_of_alleles() == variant_alleles.size())) {
        allele_list = likelihoods;
    }
    else {
        allele_list = IndexedAlleleList<pAllele>::create(variant_alleles, _pool);
    }

    auto* l = _genotyping_model->calculate_likelihoods(allele_list, _ploidy_model, likelihoods, ref, offset_for_ref_into_event, _pool);

    pGenotypesContext result = GenotypesContext::create(_pool);
    size_t sample_count = _sample_list->number_of_samples();
    for (size_t i = 0; i < sample_count; ++i) {
        pGenotype g = Genotype::create(_pool);
        g->set_alleles({size_t(_args->sample_ploidy), StaticAllele::get_instance()->_no_call.get(), _pool});
        g->set_pl(std::move(l->sample_likelihoods((int32_t)i)->_pls));
        g->set_id(int32_t(i));
        result->add(g);
    }
    return result;
}

pGenotypePriorCalculator GermlineGenotyingEngine::resolve_genotype_prior_calculator(pHCArgs args)
{
    if (nullptr == args->dragstr_params || args->dont_use_dragstr_priors) {
        return GenotypePriorCalculator::assuming_hw(args->snp_heterozygosity, args->indel_heterozygosity);
    }
    else {
        /*! @note: 暂未实现 */
        return nullptr;
    }
}

bool GermlineGenotyingEngine::cannot_be_genotyped(pVariant vc)
{
    if (vc->allele_num() <= GenotypeLikelihoods::s_max_diploid_alt_alleles_that_can_be_genotyped) {
        pGenotypesContext gc = vc->genotype();
        for (size_t i = 0, len = gc->size(); i < len; ++i) {
            if (gc->at(i)->has_likelihoods()) {
                return false;
            }
        }
    }

    if (vc->allele_num() > GenotypeLikelihoods::s_max_diploid_alt_alleles_that_can_be_genotyped) {
        RovacaLogger::warn("attempting to genotype more than 50 alleles. site will be skipped at location {}:{}-{}", vc->get_tid(),
                         vc->get_start(), vc->get_stop());
    }
    else {
        RovacaLogger::warn(
            "no genotype contained sufficient data to recalculate site and allele qualities. site will be skipped at location {}:{}-{}",
            vc->get_tid(), vc->get_start(), vc->get_stop());
    }
    return true;
}

bool GermlineGenotyingEngine::is_spanning_deletion(pAllele allele)
{
    return allele->equals(*StaticAllele::get_instance()->_span_del) ||
           allele->equals(*StaticAllele::get_instance()->_spanning_deletion_symbolic_allele_deprecated);
}

bool GermlineGenotyingEngine::no_alleles_or_first_allele_is_not_non_ref(const AlleleVector& alleles)
{
    return alleles.empty() || !alleles.front()->equals(*StaticAllele::get_instance()->_non_ref_allele);
}

bool GermlineGenotyingEngine::force_keep_allele(pAllele allele)
{
    return allele->equals(*StaticAllele::get_instance()->_non_ref_allele, false) ||
           _args->reference_confidence_mode != ReferenceConfidenceMode::NONE;
}

bool GermlineGenotyingEngine::emit_all_active_sites() const { return _args->output_mode == EMIT_ALL_ACTIVE_SITES; }

bool GermlineGenotyingEngine::passes_call_threshold(double conf) const { return conf >= _args->standard_confidence_for_calling; }

bool GermlineGenotyingEngine::passes_emit_threshold(double conf, bool best_guess_is_ref) const
{
    return (_args->output_mode == EMIT_ALL_CONFIDENT_SITES || !best_guess_is_ref) && passes_call_threshold(conf);
}

pVariant GermlineGenotyingEngine::make_annotated_call(pRefFragment ref, pVariant call, pRALikelihoods ra_likelihoods,
                                                      size_t num_before_remove)
{
    pVariant untrimmed_result = _annotator_engine->annotate_context(ref, call, ra_likelihoods, db_offset_, db_data_, _pool);

    if (untrimmed_result->allele_num() == num_before_remove) {
        return untrimmed_result;
    }

    return ROVACAVariantContextUtils::reverse_trim_alleles(untrimmed_result);
}

pOutputAlleleSubset GermlineGenotyingEngine::calculate_output_allele_subset(pAFCalculationResult af_calculation_result, pVariant merged_vc,
                                                                            const AlleleSet& forced_alleles)
{
    pOutputAlleleSubset result = new ALLOC_TYPE_IN_POOL(_pool, OutputAlleleSubset) OutputAlleleSubset{_pool};

    const AlleleVector& alleles = af_calculation_result->get_alleles_used_in_genotyping();
    size_t alternative_allele_count = alleles.size() - 1;

    for (pAllele a : alleles) {
        if (ROVACA_LIKELY(!a->is_reference())) {
            // we want to keep the NON_REF symbolic allele but only in the absence of a non-symbolic allele, e.g.
            // if we combined a ref / NON_REF gVCF with a ref / alt gVCF
            bool is_non_ref_which_is_lone_alt_allele =
                1 == alternative_allele_count && a->equals(*StaticAllele::get_instance()->_non_ref_allele);
            bool is_plausible = af_calculation_result->passes_threshold(a, _args->standard_confidence_for_calling);

            // it's possible that the upstream deletion that spanned this site was not emitted, mooting the symbolic spanning deletion
            // allele
            bool is_spurious_spanning_deletion = is_spanning_deletion(a) && !is_vc_covered_by_deletion(merged_vc);

            bool to_output = (is_plausible || force_keep_allele(a) || is_non_ref_which_is_lone_alt_allele || forced_alleles.count(a)) &&
                             !is_spurious_spanning_deletion;

            result->site_is_monomorphic &= !(is_plausible && !is_spurious_spanning_deletion);

            if (to_output) {
                result->alleles.emplace_back(a);
                result->mle_counts.emplace_back(af_calculation_result->get_allele_count_at_mle(a));
            }
        }
    }

    return result;
}

bool GermlineGenotyingEngine::contains_calls(const VariantVector& calls)
{
    pGenotypesContext genotypes;
    for (pVariant vc : calls) {
        genotypes = vc->genotype();
        for (size_t i = 0, glen = genotypes->size(); i < glen; ++i) {
            if (genotypes->at(i)->is_called()) {
                return true;
            }
        }
    }
    return false;
}

void GermlineGenotyingEngine::clear_upstream_deletions_loc()
{
    for (auto i : _deletions_locs) {
        delete i;
    }
    _deletions_locs.clear();
    _upstream_deletions_loc.clear();
}

void GermlineGenotyingEngine::record_deletions(pVariant vc, const AlleleVector& alleles)
{
    while (!_upstream_deletions_loc.empty() &&
           (!(*_upstream_deletions_loc.begin())->contigs_match(*vc) || (*_upstream_deletions_loc.begin())->get_stop() < vc->get_start())) {
        _upstream_deletions_loc.erase(_upstream_deletions_loc.begin());
    }

    uint32_t ref_len, alt_len;
    ref_len = vc->ref_allele()->length();
    for (pAllele allele : alleles) {
        alt_len = allele->length();
        // in a deletion
        if (ref_len > alt_len) {
            pInterfaceLocatable loc = SimpleInterval::create(vc->get_tid(), vc->get_start(), vc->get_start() + ref_len - alt_len);
            _upstream_deletions_loc.insert(loc);
            _deletions_locs.push_back(loc);
        }
    }
}

bool GermlineGenotyingEngine::is_vc_covered_by_deletion(pVariant vc)
{
    return !_upstream_deletions_loc.empty() &&
           std::any_of(_upstream_deletions_loc.begin(), _upstream_deletions_loc.end(), [&](pInterfaceLocatable loc) {
               return loc->contigs_match(*vc) && loc->get_start() < vc->get_start() && vc->get_start() <= loc->get_stop();
           });
}

bool GermlineGenotyingEngine::emit_reference_confidence() const
{
    return _args->reference_confidence_mode != ReferenceConfidenceMode::NONE;
}

void GermlineGenotyingEngine::compose_call_attributes(pVariant vc, Int32Vector& mleac, pAFCalculationResult af_result,
                                                      const AlleleVector& all_alleles_to_use, pGenotypesContext genotypes, pInfoData info)
{
    un_used(vc);
    un_used(af_result);
    un_used(all_alleles_to_use);

    if (!mleac.empty()) {
        pGenotype g;
        int32_t an = 0;
        for (size_t i = 0, len = genotypes->size(); i < len; ++i) {
            g = genotypes->at(size_t(i));
            std::for_each(g->alleles().begin(), g->alleles().end(), [&](pAllele a) { an += a->is_called(); });
        }

        DoubleVector mleaf{_pool};
        std::for_each(mleac.begin(), mleac.end(), [&](int32_t ac) { mleaf.push_back(an == 0 ? NAN : std::min(1.0, (double)ac / an)); });
        info->set_mleaf(std::move(mleaf));
        info->set_mleac(std::move(mleac));
    }
    // 以下代码 hc 未使用，暂未实现
#if 0
    if (do_allele_specific_calcs) {
        list<integer> per_allele_quals = new array_list<>();
        // per-allele quals are not calculated for biallelic sites
        if (afresult.get_alleles_used_in_genotyping().size() > 2) {
            for (final allele a : all_alleles_to_use) {
                if (a.is_non_reference()) {
                    //*-10 to convert from log10-scale to phred-scale, as quals are typically represented
                    per_allele_quals.add((int)math.round(afresult.get_log10posterior_of_allele_absent(a) * -10));
                }
            }
        }
        else {
            //*-10 to convert from log10-scale to phred-scale, as quals are typically represented
            per_allele_quals.add((int)math.round(afresult.log10prob_only_ref_allele_exists() * -10));
        }
        attributes.put(gatkvcfconstants.as_qual_key,
                       per_allele_quals.stream().map_to_int(q->math.round(q)).boxed().collect(collectors.to_list()));
    }
    if (configuration.genotype_args.annotate_number_of_alleles_discovered) {
        attributes.put(gatkvcfconstants.number_of_discovered_alleles_key, vc.get_alternate_alleles().size());
    }
#endif
}

pRALikelihoods GermlineGenotyingEngine::prepare_read_allele_likelihoods_for_annotation(pRHLikelihoods rh_likelihoods,
                                                                                       const Int32ToReadVectorMap& filtered_read, bool erc,
                                                                                       const AlleleToHaplotypeVectorMap& allele_mapper,
                                                                                       pRALikelihoods ra_likelihoods, pVariant call,
                                                                                       pSimpleInterval loc)
{
    un_used(rh_likelihoods);
    un_used(erc);
    un_used(allele_mapper);

    pRALikelihoods ra_likelihoods_for_annotations;

    // 根据是否进行了污染过滤来决定是否需要重新计算基因型似然值
    if (ROVACA_LIKELY(_args->use_filtered_read_map_for_annotations || !_args->is_sample_contamination_present())) {
        ra_likelihoods_for_annotations = ra_likelihoods;
    }
    else {
        // todo: 以下代码 hc 未使用，暂未实现
        ra_likelihoods_for_annotations = ra_likelihoods;
        RovacaLogger::error("invalid site, hc not to do");
    }

    if (call->allele_num() != ra_likelihoods_for_annotations->number_of_alleles()) {
        auto* new_allele_list = IndexedAlleleList<pAllele>::create(call->alleles(), _pool);
        ra_likelihoods_for_annotations->update_non_ref_allele_likelihoods(new_allele_list);
    }

    /*!
     * @brief 为了减少下一步过滤操作的计算量，只对与当前变异位点重叠的reads进行过滤.
     * 从每个样本的过滤后的reads列表中，筛选出与当前变异位点重叠的reads.
     */
    Int32ToReadVectorMap filters = overlapping_filtered_reads(filtered_read, loc);
    ra_likelihoods_for_annotations->add_evidence(filters, 0.0);

    return ra_likelihoods_for_annotations;
}

Int32ToReadVectorMap GermlineGenotyingEngine::overlapping_filtered_reads(const Int32ToReadVectorMap& filtered_read, pSimpleInterval loc)
{
    Int32ToReadVectorMap sample_reads;
    for (const auto& tup : filtered_read) {
        int32_t sample_index = tup.first;
        const ReadVector& reads = tup.second;
        if (reads.empty()) {
            continue;
        }

        ReadVector filters(_pool);
        std::for_each(reads.begin(), reads.end(), [&](pReadRecord r) {
            if (r->overlaps(*loc)) {
                filters.emplace_back(r);
            }
        });

        if (!filters.empty()) {
            sample_reads.insert({sample_index, std::move(filters)});
        }
    }
    return sample_reads;
}

pVariant GermlineGenotyingEngine::remove_alt_alleles_if_too_many_genotypes(int32_t ploidy, const AlleleToHaplotypeVectorMap& allele_mapper,
                                                                           pVariant merged_vc)
{
    un_used(ploidy);
    un_used(allele_mapper);
    un_used(merged_vc);
    // RovacaLogger::warn("too many alt allele at {}:{}-{}", merged_vc->get_tid(), merged_vc->get_start(), merged_vc->get_stop());
    return nullptr;
}

void GermlineGenotyingEngine::init_engine()
{
    _annotator_engine = new VariantAnnotatorEngine{_args->reference_confidence_mode == ReferenceConfidenceMode::GVCF};
    _prior_calculator = GenotypePriorCalculator::assuming_hw(_args->snp_heterozygosity, _args->indel_heterozygosity);
    _genotyping_model = new IndependentSampleGenotypesModel{};
    _allele_frequency_calculator = AlleleFrequencyCalculator::make_calculator(
        _args->sample_ploidy, _args->snp_heterozygosity, _args->indel_heterozygosity, _args->heterozygosity_standard_deviation);
}

}  // namespace rovaca