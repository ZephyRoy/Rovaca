#include "assembly_based_caller_utils.h"

#include <algorithm>

#include "alignment_utils.h"
#include "allele_likelihoods.hpp"
#include "cigar_utils.h"
#include "rovaca_logger.h"
#include "event_map.h"
#include "genotype.h"
#include "genotypes_context.hpp"
#include "haplotype.h"
#include "rovaca_variant_context_utils.h"
#include "pileup_element.h"
#include "read_pileup.h"
#include "read_record.h"
#include "read_record_utils.h"
#include "variant.h"

namespace rovaca
{

const char* PHASE_01_DESCRIPTION = "0|1";
const char* PHASE_10_DESCRIPTION = "1|0";
static const size_t PHASESTR_LEN = strlen(PHASE_10_DESCRIPTION);

static constexpr int64_t s_default_adaptor_size = 100;

#define GET_DESCRIPTION(idx) (idx == PhaseGroup::PHASE_10 ? PHASE_10_DESCRIPTION : PHASE_01_DESCRIPTION)

AlleleToHaplotypeVectorMap AssemblyBasedCallerUtils::create_allele_mapper(const HaplotypeVector& haplotypes, pVariant merged, int64_t loc,
                                                                          bool emit_spanning_dels)
{
    AlleleToHaplotypeVectorMap result(haplotypes.get_allocator());

    pAllele ref = merged->ref_allele();
    uint32_t span_ref_len, merged_ref_len = ref->length();
    const AlleleVector& alleles = merged->alleles();
    std::for_each(alleles.begin(), alleles.end(), [&result](pAllele a) { result.insert({a, {}}); });

    for (pHaplotype h : haplotypes) {
        VariantVector spanning_events = h->event_map()->get_overlapping_events(loc);

        // no events --> this haplotype supports the reference at this locus
        if (spanning_events.empty()) {
            result.at(ref).emplace_back(h);
            continue;
        }

        for (pVariant spanning_event : spanning_events) {
            if (spanning_event->get_start() == loc) {
                // the event starts at the current location
                span_ref_len = spanning_event->ref_allele()->length();
                if (span_ref_len == merged_ref_len) {
                    // reference allele lengths are equal; we can just use the spanning event's alt allele in the case of GGA mode the
                    // spanning event might not match an allele in the mergedVC
                    if (result.count(spanning_event->alternate_allele_at(0))) {
                        // variant contexts derived from the event map have only one alt allele each, so we can just grab the first one
                        // (we're not assuming that the sample is biallelic)
                        result.at(spanning_event->alternate_allele_at(0)).emplace_back(h);
                    }
                }
                else if (span_ref_len < merged_ref_len) {
                    // spanning event has shorter ref allele than merged VC; we need to pad out its alt allele
                    AlleleMap allele_map(haplotypes.get_allocator().resource());
                    ROVACAVariantContextUtils::create_allele_mapping(ref, spanning_event, allele_map);

                    pAllele remapped_spanning_event_alt_allele = allele_map.old2new_map.at(spanning_event->alternate_allele_at(0));

                    // in the case of GGA mode the spanning event might not match an allele in the mergedVC
                    if (result.count(remapped_spanning_event_alt_allele)) {
                        result.at(remapped_spanning_event_alt_allele).emplace_back(h);
                    }
                }
                else {
                    // the process of creating the merged VC in AssemblyBasedCallerUtils::makeMergedVariantContext should have already
                    // padded out the reference allele, therefore this spanning VC must not be in events at this site because we're in GGA
                    // mode and it's not an allele we want
                    continue;
                }
            }
            else {
                if (emit_spanning_dels) {
                    // the event starts prior to the current location, so it's a spanning deletion
                    if (!result.count(StaticAllele::get_instance()->_span_del.get())) {
                        result.insert({StaticAllele::get_instance()->_span_del.get(), {}});
                    }
                    result.at(StaticAllele::get_instance()->_span_del.get()).emplace_back(h);
                    // there might be a del+ins at the site in question and this would miss one of them unless its a continue
                    break;  // Why is there a break here? Shouldn't this be a continue? Why should the first spanning event overlap?A
                }
                else {
                    result.at(ref).emplace_back(h);
                    break;
                }
            }
        }
    }

    return result;
}

VariantVector AssemblyBasedCallerUtils::get_variant_contexts_from_active_haplotypes(int64_t loc, const HaplotypeVector& haplotypes,
                                                                                    bool include_spanning_events)
{
    std::pmr::list<pVariant> result(haplotypes.get_allocator());
    std::pmr::unordered_set<pVariant, VariantHash, VariantEqual> unique_loc(haplotypes.get_allocator());

    for (pHaplotype h : haplotypes) {
        VariantVector variants_in_map = h->event_map()->get_overlapping_events(loc);
        for (pVariant vc : variants_in_map) {
            if (include_spanning_events || vc->get_start() == loc) {
                if (!unique_loc.count(vc)) {
                    unique_loc.insert(vc);
                    result.emplace_back(vc);
                }
            }
        }
    }

    return {{result.begin(), result.end()}, haplotypes.get_allocator()};
}

pVariant AssemblyBasedCallerUtils::make_merged_variant_context(const VariantVector& vcs)
{
    if (vcs.empty()) {
        return nullptr;
    }

    return ROVACAVariantContextUtils::simple_merge(vcs, KEEP_IF_ANY_UNFILTERED, PRIORITIZE, false);
}

AlleleSet AssemblyBasedCallerUtils::get_alleles_consistent_with_given_alleles(const AlleleVector& given_alleles, pVariant merged_vc)
{
    un_used(merged_vc);

    if (given_alleles.empty()) {
        return {};
    }

    // hc 中 given_alleles 一直为空，下面代码暂未实现
    RovacaLogger::error("invalid code");
    exit(EXIT_FAILURE);
}

void AssemblyBasedCallerUtils::realign_reads_to_their_best_haplotype(pRHLikelihoods read_likelihoods, pHaplotype reference,
                                                                     int64_t padded_reference_loc_start, p_lib_sw_avx sw, pMemoryPool pool,
                                                                     pBamDataPool bam_pool)
{
    auto lambda_tmp = [](void* h) -> double {
        pCigar cigar = static_cast<pHaplotype>(h)->cigar();
        int32_t reference_term = (static_cast<pHaplotype>(h)->is_reference() ? 1 : 0);
        int32_t cigar_term = cigar == nullptr ? 0 : (1 - int32_t(cigar->num));
        return (double)(reference_term + cigar_term);
    };

    size_t len = read_likelihoods->number_of_samples();

    for (size_t s = 0; s < len; s++) {
        auto best_alleles = read_likelihoods->best_alleles_breaking_ties((int32_t)s, lambda_tmp);
        for (auto ba : best_alleles) {
            pReadRecord original_read = ba->evidence;
            pHaplotype best_haplotype = (pHaplotype)ba->best_allele;
            bool is_informative = ba->is_informative();
            AlignmentUtils::create_read_aligned_to_ref(original_read, best_haplotype, reference, padded_reference_loc_start, is_informative,
                                                       sw, pool, bam_pool);
        }
    }
}

ReadPileupVector AssemblyBasedCallerUtils::get_pileups_over_reference(const ReadHashSet& un_sorted_reads, pSimpleInterval active_loc,
                                                                      pBases ref, int64_t bases_offset, pMemoryPool pool)
{
    int32_t tid = active_loc->get_tid();
    int64_t active_start_pos = active_loc->get_start();
    int64_t active_stop_pos = active_loc->get_stop();
    int64_t ref_pileups_size = active_loc->get_length();

    ReadPileupVector pileups{size_t(ref_pileups_size), nullptr, pool};
    for (int64_t i = 0, pos = active_start_pos; pos <= active_stop_pos; ++pos, ++i) {
        pileups.at(i) = ReadPileup::create(tid, pos, pos, pool);
    }

    // 使用set对read进行排序
    std::pmr::set<pReadRecord> sorted_reads{pool};
    std::for_each(un_sorted_reads.begin(), un_sorted_reads.end(), [&](pReadRecord r) { sorted_reads.insert(r); });

    pReadPileup current_pileup;
    int64_t read_start_pos, ref_offset;
    uint32_t num_of_cigar, ce, ce_op, ce_length;
    int32_t adaptor_boundary;
    int64_t insert_size;
    bool is_reverse_strand;

    for (pReadRecord read : sorted_reads) {
        read_start_pos = read->get_start();
        if (read_start_pos > active_stop_pos) {
            continue;
        }

        adaptor_boundary = read->get_adaptor_bounder();
        insert_size = read->insert_size();
        is_reverse_strand = read->is_reverse_strand();
        num_of_cigar = read->cigar_length();
        ref_offset = read_start_pos - active_start_pos;

        for (uint32_t ce_idx = 0, query_offset = 0; ce_idx < num_of_cigar && ref_offset < ref_pileups_size; ++ce_idx) {
            ce = read->cigar_i(ce_idx);
            ce_op = bam_cigar_op(ce);
            ce_length = bam_cigar_oplen(ce);

            if (cigar_op_is_hard_clip(ce_op) || cigar_op_is_ins(ce_op)) {
                query_offset += cigar_op_is_ins(ce_op) ? ce_length : 0;
                continue;
            }

            if (cigar_op_is_ref_skip(ce_op)) {
                ref_offset += ce_length;
                query_offset += ce_length;
                continue;
            }

            // D M = X 仅此4种情况, 但因为 =和X 常用 M 表示，故累加 qual 时需要比较 ref_base 与 read_base 来区分
            if (consumes_ref_bases(ce_op)) {
                for (uint32_t i = 0; i < ce_length; i++) {
                    if (dont_include_read_in_pileup(adaptor_boundary, insert_size, is_reverse_strand, active_start_pos + ref_offset + i)) {
                        continue;
                    }
                    if (ref_offset + i >= 0 && ref_offset + i < ref_pileups_size) {
                        current_pileup = pileups.at(ref_offset + i);
                        current_pileup->add(PileupElement::create(read, query_offset + i, ce, ce_idx, i, pool));

                        // logic like active Region
                        uint8_t qual = cigar_op_is_del(ce_op) ? REF_MODEL_DELETION_QUAL : read->qual_i(query_offset + i);
                        if (cigar_op_is_del(ce_op) || qual > BASE_QUAL_THRESHOLD) {
                            //  accumulate qual to his
                            if (cigar_op_is_del(ce_op) || ref->data[bases_offset + ref_offset + i] != read->seq_i(query_offset + i)) {
                                current_pileup->qual_his[FS_NON_REF][qual]++;
                                current_pileup->qual_max[FS_NON_REF] = std::max(current_pileup->qual_max[FS_NON_REF], qual);
                            }
                            else {
                                current_pileup->qual_his[FS_REF][qual]++;
                                current_pileup->qual_max[FS_REF] = std::max(current_pileup->qual_max[FS_REF], qual);
                            }
                        }
                    }
                }
            }
            ref_offset += consumes_ref_bases(ce_op) ? ce_length : 0;
            query_offset += consumes_read_bases(ce_op) ? ce_length : 0;
        }
    }

    return pileups;
}

bool AssemblyBasedCallerUtils::dont_include_read_in_pileup(int32_t adaptor_boundary, int64_t insert_size, bool is_reverse_strand,
                                                           int64_t pos)
{
    if (s_cannot_compute_adaptor_boundary == adaptor_boundary || insert_size > s_default_adaptor_size) {
        return false;
    }
    auto pos_ = static_cast<int64_t>(pos);
    return is_reverse_strand ? pos_ <= adaptor_boundary : pos_ >= adaptor_boundary;
}

VariantVector AssemblyBasedCallerUtils::phase_calls(const VariantVector& calls, const HaplotypeList& called_haplotypes)
{
    if (calls.empty()) {
        return {{}, calls.get_allocator()};
    }

    VariantToHaplotypeSetMap haplotype_map = construct_haplotype_mapping(calls, called_haplotypes);

    VariantToInt32PhaseGroupPairMap phase_set_map = construct_phase_set_mapping(calls, haplotype_map);

    std::pmr::set<int32_t> distinct_set(calls.get_allocator());
    for (const auto& item : phase_set_map) {
        distinct_set.insert(item.second.first);
    }
    int32_t unique_counter_end_value = (int32_t)distinct_set.size();

    return construct_phase_groups(calls, phase_set_map, unique_counter_end_value);
}

VariantToHaplotypeSetMap AssemblyBasedCallerUtils::construct_haplotype_mapping(const VariantVector& calls,
                                                                               const HaplotypeList& called_haplotypes)
{
    VariantToHaplotypeSetMap result(calls.get_allocator());

    for (pVariant call : calls) {
        if (!is_biallelic_with_one_site_specific_alternate_allele(call)) {
            result.insert({call, {}});
            continue;
        }

        pAllele alt = get_site_specific_alternate_allele(call);
        HaplotypeSet haps_with_allele(calls.get_allocator());
        for (pHaplotype h : called_haplotypes) {
            const Int64ToVariantMap& events = h->event_map()->get_events();
            for (const auto& tup : events) {
                pVariant vc = tup.second;
                const AlleleVector& vc_alt_alleles = vc->alt_alleles();
                bool has_alt = std::find_if(vc_alt_alleles.begin(), vc_alt_alleles.end(), [&](pAllele a) { return a->equals(*alt); }) !=
                               std::end(vc_alt_alleles);
                if (has_alt && vc->get_start() == call->get_start()) {
                    haps_with_allele.insert(h);
                }
            }
        }

        result.insert({call, std::move(haps_with_allele)});
    }
    return result;
}

VariantToInt32PhaseGroupPairMap AssemblyBasedCallerUtils::construct_phase_set_mapping(const VariantVector& calls,
                                                                                      const VariantToHaplotypeSetMap& haplotype_map)
{
    auto alloc = calls.get_allocator();

    HaplotypeSet haplotypes_with_called_variants(alloc);
    for (const auto& item : haplotype_map) {
        haplotypes_with_called_variants.insert(item.second.begin(), item.second.end());
    }

    int32_t unique_counter = 0;
    size_t num_calls = calls.size();
    size_t total_available_haplotypes = haplotypes_with_called_variants.size();

    VariantToInt32PhaseGroupPairMap phase_set_mapping(alloc);

    pVariant call, comp;
    for (size_t i = 0; i < num_calls - 1; ++i) {
        call = calls.at(i);
        const HaplotypeSet& haplotypes_with_call = haplotype_map.at(call);
        if (haplotypes_with_call.empty()) {
            continue;
        }

        // 这个变异出现在所有单倍体上，即所有单倍体都含有这个变异
        bool call_is_on_all_alt_haps = haplotypes_with_call.size() == total_available_haplotypes;

        /*!
         * call_is_on_all_alt_haps为true并不一定意味着该变异位点是纯合子基因型，因为该方法不考虑参考单倍型。
         * 此外，该代码段还维护了一个集合call_haplotypes_available_for_phasing，用于跟踪哪些单倍型可用于下游相位化。
         * 如果变异位点在所有备用等位基因的单倍型上，但我们只使用其中一些单倍型来相位化它，那么我们需要跟踪哪些单倍型仍然可以用于下游相位化
         */
        HaplotypeSet call_haplotypes_available_for_phasing{haplotypes_with_call, alloc};

        for (size_t j = i + 1; j < num_calls; ++j) {
            comp = calls.at(j);
            const HaplotypeSet& haplotypes_with_comp = haplotype_map.at(comp);
            if (haplotypes_with_comp.empty()) {
                continue;
            }

            bool comp_is_on_all_alt_haps = haplotypes_with_comp.size() == total_available_haplotypes;

            if ((haplotypes_with_call.size() == haplotypes_with_comp.size() && contains_all(haplotypes_with_call, haplotypes_with_comp)) ||
                (call_is_on_all_alt_haps && contains_all(call_haplotypes_available_for_phasing, haplotypes_with_comp)) ||
                comp_is_on_all_alt_haps) {
                /*!
                 * 主要用于将变异位点和比较位点分配到一个相位组中
                 * 如果变异位点和比较位点还没有被分配到任何相位组中，则创建一个新的相位组，并将它们分配到该相位组中。
                 * 如果比较位点已经被分配到了某个相位组中，则说明存在另一个变异位点与比较位点相位一致，但与当前变异位点不相位一致，这种情况下，我们需要中止相位化过程。
                 * 如果变异位点和比较位点都还没有被分配到任何相位组中，则将它们分配到一个新的相位组中，并将该相位组标记为"0|1"相位。
                 * 需要注意的是，即使变异位点是纯合子基因型，我们也将其标记为"0|1"相位，因为我们无法确定下游步骤是否会改变基因型。
                 * 因此，我们在这里假设变异位点是杂合子基因型，并在下游步骤中根据需要进行更正。
                 * 最后，我们更新一个集合call_haplotypes_available_for_phasing，用于跟踪哪些单倍型可用于下游相位化。
                 * 如果变异位点是所有备用等位基因的单倍型，则我们需要缩小可用单倍型的范围，以便更准确地相位化其他变异位点。
                 */
                if (!phase_set_mapping.count(call)) {
                    if (phase_set_mapping.count(comp)) {
                        phase_set_mapping.clear();
                        return phase_set_mapping;
                    }
                    phase_set_mapping.insert({call, {unique_counter, PhaseGroup::PHASE_01}});
                    phase_set_mapping.insert({comp, {unique_counter, PhaseGroup::PHASE_01}});
                    retain_all(call_haplotypes_available_for_phasing, haplotypes_with_comp);
                    unique_counter++;
                }
                else if (!phase_set_mapping.count(comp)) {
                    const auto& call_phase = phase_set_mapping.at(call);
                    phase_set_mapping.insert({comp, {call_phase.first, call_phase.second}});
                }
            }
            else if (haplotypes_with_call.size() + haplotypes_with_comp.size() == total_available_haplotypes) {
                /*!
                 * 将一组相同的位点分配到同一个相位组中，并为每个相位组分配一个唯一的id。
                 * 首先判断当前位点是否与所有备选单倍体上的位点都不在同一相位组中。如果是，则为该位点创建一个新的相位组，并为该相位组分配一个唯一的id。
                 * 如果该位点与某些备选单倍体上的位点在同一相位组中，则将该位点分配到已有的相位组中，并使用该相位组的id。
                 */
                HaplotypeSet intersection{haplotypes_with_call, alloc};
                retain_all(intersection, haplotypes_with_comp);
                if (intersection.empty()) {
                    if (!phase_set_mapping.count(call)) {
                        if (phase_set_mapping.count(comp)) {
                            phase_set_mapping.clear();
                            return phase_set_mapping;
                        }
                        phase_set_mapping.insert({call, {unique_counter, PhaseGroup::PHASE_01}});
                        phase_set_mapping.insert({comp, {unique_counter, PhaseGroup::PHASE_10}});
                        unique_counter++;
                    }
                    else if (!phase_set_mapping.count(comp)) {
                        const auto& call_phase = phase_set_mapping.at(call);
                        Int32PhaseGroupPair p{};
                        p.first = call_phase.first;
                        p.second = call_phase.second == PhaseGroup::PHASE_01 ? PhaseGroup::PHASE_10 : PhaseGroup::PHASE_01;
                        phase_set_mapping.insert({comp, p});
                    }
                }
            }
        }
    }

    return phase_set_mapping;
}

VariantVector AssemblyBasedCallerUtils::construct_phase_groups(const VariantVector& calls,
                                                               const VariantToInt32PhaseGroupPairMap& phase_set_mapping, int32_t index_to)
{
    auto alloc = calls.get_allocator();
    VariantVector phased_calls{calls, alloc};
    Int32Vector indexes{calls.size(), alloc};

    for (int32_t count = 0; count < index_to; ++count) {
        indexes.clear();
        for (int32_t index = 0, len = (int32_t)calls.size(); index < len; ++index) {
            pVariant call = calls.at(index);
            if (phase_set_mapping.count(call) && phase_set_mapping.at(call).first == count) {
                indexes.push_back(index);
            }
        }

        CHECK_CONDITION_EXIT(indexes.size() < 2, "somehow we have a group of phased variants that has fewer than 2 members");

        pBases unique_id = create_unique_id(calls.at(indexes.at(0)));

        int32_t phase_set_id = (int32_t)calls.at(indexes.at(0))->get_start();

        for (int32_t index : indexes) {
            pVariant original_call = calls.at(index);
            phased_calls[index] = phase_vc(original_call, unique_id, phase_set_mapping.at(original_call).second, phase_set_id);
        }
    }

    return phased_calls;
}

pBases AssemblyBasedCallerUtils::create_unique_id(pVariant vc)
{
    pMemoryPool pool = vc->alleles().get_allocator().resource();
    int64_t start = vc->get_start();
    pAllele ref = vc->ref_allele();
    pAllele first_alt = vc->alternate_allele_at(0);
    uint32_t digits = uint32_t(std::to_string(start).size()) + 2 + ref->length() + first_alt->length();
    pBases id = new ALLOC_FLEXIBLE_IN_POOL(pool, Bases, digits, uint8_t) Bases{digits};
    // note: 此处可能会有内存越界问题
    sprintf((char*)id->data, "%ld_%s_%s", start, ref->get_display_string()->data, first_alt->get_display_string()->data);
    return id;
}

pVariant AssemblyBasedCallerUtils::phase_vc(pVariant vc, pBases id, PhaseGroup phase_gt, int32_t phase_set_id)
{
    pGenotype g;
    pGenotypesContext genotypes = vc->genotype();
    pMemoryPool pool = vc->alleles().get_allocator().resource();
    for (size_t i = 0, len = genotypes->size(); i < len; ++i) {
        g = genotypes->at(i);
        AlleleVector& alleles_in_genotype = const_cast<AlleleVector&>(g->alleles());
        if (g->is_het() && !is_site_specific_alt_allele(alleles_in_genotype.at(static_cast<int32_t>(phase_gt)))) {
            std::reverse(alleles_in_genotype.begin(), alleles_in_genotype.end());
        }

        pPhasedData data = new ALLOC_TYPE_IN_POOL(pool, PhasedData) PhasedData{};
        data->pid = id;
        data->pgt = GET_DESCRIPTION(phase_gt);
        data->ps = phase_set_id;
        g->set_phased(data);
    }

    return vc;
}

bool AssemblyBasedCallerUtils::is_site_specific_alt_allele(pAllele a)
{
    return !(a->is_reference() || a->is_non_ref_allele() || a->equals(*StaticAllele::get_instance()->_span_del));
}

bool AssemblyBasedCallerUtils::is_biallelic_with_one_site_specific_alternate_allele(pVariant vc)
{
    const AlleleVector& alleles = vc->alleles();
    int32_t count = 0;
    std::for_each(alleles.begin(), alleles.end(), [&](pAllele a) { count += is_site_specific_alt_allele(a) ? 1 : 0; });
    return 1 == count;
}

pAllele AssemblyBasedCallerUtils::get_site_specific_alternate_allele(pVariant vc)
{
    const AlleleVector& alt_alleles = vc->alt_alleles();
    for (const auto& item : alt_alleles) {
        if (is_site_specific_alt_allele(item)) {
            return item;
        }
    }
    return nullptr;
}

bool AssemblyBasedCallerUtils::contains_all(const HaplotypeSet& big, const HaplotypeSet& small)
{
    if (small.size() > big.size()) {
        return false;
    }
    for (const auto& h : small) {
        if (!big.count(h)) {
            return false;
        }
    }
    return true;
}

void AssemblyBasedCallerUtils::retain_all(HaplotypeSet& src, const HaplotypeSet& target)
{
    for (pHaplotype h : src) {
        if (!target.count(h)) {
            src.erase(h);
        }
    }
}

pHaplotype AssemblyBasedCallerUtils::create_reference_haplotype(pRefFragment ref, pSimpleInterval ref_loc, pSimpleInterval padded,
                                                                pMemoryPool pool)
{
    CHECK_CONDITION_EXIT(!ref_loc->contains(*padded), "ref_loc must contains flank_padded");

    int64_t ref_start = padded->get_start() - ref_loc->get_start();
    uint32_t bases_len = uint32_t(padded->get_length_on_reference());

    pBases bases = new ALLOC_FLEXIBLE_IN_POOL(pool, Bases, bases_len, uint8_t) Bases{bases_len};
    memcpy(bases->data, ref->data + ref_start, sizeof(uint8_t) * bases_len);

    pCigar cigar = new ALLOC_FLEXIBLE_IN_POOL(pool, Cigar, 1, uint32_t) Cigar{1};
    cigar->data[0] = bases_len << BAM_CIGAR_SHIFT | BAM_CMATCH;

    pHaplotype h = Haplotype::create(pool);
    h->set_kmer_size(0);
    h->set_cigar(cigar);
    h->init_haplotype(bases, 1);
    h->set_alignment_start_hap_wrt_ref(0);
    h->set_interval(padded);
    return h;
}

}  // namespace rovaca