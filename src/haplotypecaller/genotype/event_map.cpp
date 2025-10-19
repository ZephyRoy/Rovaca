#include "event_map.h"

#include <algorithm>

#include "haplotype.h"
#include "htslib/sam.h"
#include "rovaca_logger.h"
#include "simple_interval.h"
#include "utils/base_utils.h"
#include "variant.h"

namespace rovaca
{

static constexpr size_t s_default_ploidy = 2;
static constexpr size_t s_default_size = 16;

Int64Set EventMap::build_event_maps_for_haplotypes(HaplotypeVector& haplotypes, pRefFragment ref, pSimpleInterval ref_loc,
                                                   uint32_t max_mnp_distance, pMemoryPool pool)
{
    pEventMap event;
    int32_t source_idx = 0;
    Int64Set start_pos_key_set(pool);

    for (pHaplotype h : haplotypes) {
        event = new ALLOC_TYPE_IN_POOL(pool, EventMap) EventMap(source_idx++, pool);
        event->process_cigar_for_initial_events(max_mnp_distance, ref, ref_loc, h);

        for (const auto& tup : event->_events) {
            start_pos_key_set.insert(tup.first);

            // Assert that all of the events discovered have 2 alleles
            CHECK_CONDITION_EXIT(tup.second->allele_num() != s_default_ploidy, "event map variant has too many alleles, h index {}",
                                 source_idx);
        }

        h->set_event_map(event);
    }

    return start_pos_key_set;
}

VariantVector EventMap::get_overlapping_events(int64_t loc) const
{
    VariantVector overlapping_events(_pool);
    VariantVector deletion_events_ending_at_loc(_pool);

    for (auto itr = _events.begin(), upper_bound_itr = _events.upper_bound(loc); itr != upper_bound_itr; ++itr) {
        pVariant vc = itr->second;
        if (vc->get_stop() >= loc) {
            overlapping_events.emplace_back(vc);
            if (vc->is_simple_deletion() && vc->get_stop() == loc) {
                deletion_events_ending_at_loc.emplace_back(vc);
            }
        }
    }

    bool contains_deletion_ending_at_loc = !deletion_events_ending_at_loc.empty();
    bool contains_insertion_at_loc =
        std::any_of(overlapping_events.begin(), overlapping_events.end(), [](pVariant v) { return v->is_simple_insertion(); });

    if (contains_deletion_ending_at_loc && contains_insertion_at_loc) {
        // we are at the end of a deletion and the start of an insertion;
        // only the insertion should be kept in this case.
        auto itr = overlapping_events.begin();
        for (; itr != overlapping_events.end(); ++itr) {
            if (itr.operator*() == deletion_events_ending_at_loc.at(0)) {
                overlapping_events.erase(itr);
                break;
            }
        }
    }

    return overlapping_events;
}

pVariant EventMap::make_block(pVariant vc1, pVariant vc2)
{
    CHECK_CONDITION_EXIT(vc1->get_start() != vc2->get_start(), "vc1 and 2 must have the same start but got {} and {}", vc1->get_start(),
                         vc2->get_start());
    CHECK_CONDITION_EXIT(!vc1->is_biallelic(), "vc1 must be biallelic");
    if (!vc1->is_snp()) {
        CHECK_CONDITION_EXIT(
            !((vc1->is_simple_deletion() && vc2->is_simple_insertion()) || (vc1->is_simple_insertion() && vc2->is_simple_deletion())),
            "can only merge single insertion with deletion");
    }
    else {
        CHECK_CONDITION_EXIT(vc2->is_snp(), "vc1 is snp, vc2 must be a insertion or deletion");
    }

    pAllele ref, alt;
    int64_t new_stop = INVALID_INT;
    if (vc1->is_snp()) {
        // we have to repair the first base, so SNP case is special cased
        if (vc1->ref_allele()->equals(*vc2->ref_allele())) {
            // we've got an insertion, so we just update the alt to have the prev alt
            ref = vc1->ref_allele();
            pBases vc1_alt_bases = vc1->alternate_allele_at(0)->get_display_string();
            pBases vc2_alt_bases = vc2->alternate_allele_at(0)->get_display_string();
            uint32_t alt_bases_len = vc2->alternate_allele_at(0)->length();
            auto new_alt_bases = new ALLOC_FLEXIBLE_IN_POOL(_pool, Bases, alt_bases_len, uint8_t) Bases{alt_bases_len};
            new_alt_bases->data[0] = vc1_alt_bases->data[0];
            memcpy(new_alt_bases->data + 1, vc2_alt_bases->data + 1, (alt_bases_len - 1) * sizeof(uint8_t));
            alt = Allele::create_allele(new_alt_bases, false, _pool);
        }
        else {
            // we're dealing with a deletion, so we patch the ref
            ref = vc2->ref_allele();
            alt = vc1->alternate_allele_at(0);
            new_stop = vc2->get_stop();
        }
    }
    else {
        pVariant insertion = vc1->is_simple_insertion() ? vc1 : vc2;
        pVariant deletion = vc1->is_simple_insertion() ? vc2 : vc1;
        ref = deletion->ref_allele();
        alt = insertion->alternate_allele_at(0);
        new_stop = deletion->get_stop();
    }

    std::pmr::vector<pAllele> new_alleles;
    new_alleles.emplace_back(ref);
    new_alleles.emplace_back(alt);

    pVariant new_vc = Variant::create(_pool);
    new_vc->set_source_id(vc1->source_id());
    new_vc->set_tid(vc1->get_tid());
    new_vc->set_start(vc1->get_start());
    new_vc->set_stop(new_stop == INVALID_INT ? vc1->get_stop() : new_stop);
    new_vc->set_alleles(new_alleles);
    return new_vc;
}

void EventMap::process_cigar_for_initial_events(uint32_t max_mnp_distance, pRefFragment ref, pSimpleInterval ref_loc, pHaplotype h)
{
    int64_t ref_pos = h->alignment_start_hap_wrt_ref();
    if (ref_pos < 0) {
        return;  // protection against SW failures
    }

    std::pmr::vector<pAllele> alleles(_pool);
    std::pmr::vector<pVariant> proposed_events(_pool);
    std::pmr::vector<uint32_t> mismatch_offsets(_pool);
    alleles.reserve(s_default_size);
    proposed_events.reserve(s_default_size);
    mismatch_offsets.reserve(s_default_size);

    pCigar cigar = h->cigar();
    int32_t tid = ref_loc->get_tid();
    pBases alignment = h->get_bases();
    uint32_t alignment_pos = 0, cigar_idx, ce, op, op_len;

    for (cigar_idx = 0; cigar_idx < cigar->num; ++cigar_idx) {
        ce = cigar->data[cigar_idx];
        op = bam_cigar_op(ce);
        op_len = bam_cigar_oplen(ce);
        alleles.clear();
        switch (op) {
            case BAM_CINS: {
                if (ref_pos > 0) {
                    int64_t insertion_start = ref_loc->get_start() + ref_pos - 1;
                    uint8_t ref_byte = ref->data[ref_pos - 1];
                    if (BaseUtils::is_regular_base(ref_byte)) {
                        pAllele ref_allele = Allele::create_allele(ref_byte, true);
                        alleles.emplace_back(ref_allele);
                    }
                    if (0 == cigar_idx || cigar->num - 1 == cigar_idx) {
                        // if the insertion isn't completely resolved in the haplotype, skip it note this used to emit
                        // SYMBOLIC_UNASSEMBLED_EVENT_ALLELE but that seems dangerous
                    }
                    else {
                        auto insertion_bases = new ALLOC_FLEXIBLE_IN_POOL(_pool, Bases, op_len + 1, uint8_t) Bases{op_len + 1};
                        insertion_bases->data[0] = ref_byte;
                        memcpy(insertion_bases->data + 1, alignment->data + alignment_pos, sizeof(uint8_t) * op_len);
                        if (BaseUtils::is_all_regular_bases(insertion_bases)) {
                            pAllele alt_allele = Allele::create_allele(insertion_bases, false, _pool);
                            alleles.emplace_back(alt_allele);
                        }
                    }
                    if (2 == alleles.size()) {
                        pVariant ins_variant = Variant::create(_pool);
                        ins_variant->set_source_id(_source_index);
                        ins_variant->set_tid(tid);
                        ins_variant->set_start(insertion_start);
                        ins_variant->set_stop(insertion_start);
                        ins_variant->set_alleles(alleles);
                        proposed_events.emplace_back(ins_variant);
                    }
                }
                alignment_pos += op_len;
                break;
            }
            case BAM_CSOFT_CLIP: {
                alignment_pos += op_len;
                break;
            }
            case BAM_CDEL: {
                if (ref_pos > 0) {
                    pBases deletion_bases = BaseUtils::bases_create((char*)ref->data + ref_pos - 1, op_len + 1, _pool);
                    int64_t deletion_start = ref_loc->get_start() + ref_pos - 1;
                    uint8_t ref_byte = ref->data[ref_pos - 1];
                    if (BaseUtils::is_regular_base(ref_byte) && BaseUtils::is_all_regular_bases(deletion_bases)) {
                        pAllele ref_allele = Allele::create_allele(deletion_bases, true, _pool);
                        pAllele alt_allele = Allele::create_allele(ref_byte, false);
                        alleles.emplace_back(ref_allele);
                        alleles.emplace_back(alt_allele);

                        pVariant del_variant = Variant::create(_pool);
                        del_variant->set_source_id(_source_index);
                        del_variant->set_tid(tid);
                        del_variant->set_start(deletion_start);
                        del_variant->set_stop(deletion_start + op_len);
                        del_variant->set_alleles(alleles);
                        proposed_events.emplace_back(del_variant);
                    }
                }
                ref_pos += op_len;
                break;
            }
            case BAM_CMATCH:
            case BAM_CEQUAL:
            case BAM_CDIFF: {
                mismatch_offsets.clear();
                uint8_t ref_byte, alt_byte;
                uint32_t offset, start, end;
                pAllele ref_allele, alt_allele;
                pVariant mismatch_vc;
                for (offset = 0; offset < op_len; ++offset) {
                    ref_byte = ref->data[ref_pos + offset];
                    alt_byte = alignment->data[alignment_pos + offset];
                    if (ref_byte != alt_byte && BaseUtils::is_regular_base(ref_byte) && BaseUtils::is_regular_base(alt_byte)) {
                        mismatch_offsets.push_back(offset);
                    }
                }

                std::reverse(mismatch_offsets.begin(), mismatch_offsets.end());

                while (!mismatch_offsets.empty()) {
                    alleles.clear();
                    end = start = mismatch_offsets.back();
                    mismatch_offsets.pop_back();
                    while (!mismatch_offsets.empty() && mismatch_offsets.back() - end <= max_mnp_distance) {
                        end = mismatch_offsets.back();
                        mismatch_offsets.pop_back();
                    }

                    // max_mnp_distance为0时永远成立
                    if (ROVACA_LIKELY(start == end)) {
                        ref_allele = Allele::create_allele(ref->data[ref_pos + start], true);
                        alt_allele = Allele::create_allele(alignment->data[alignment_pos + start], false);
                    }
                    else {
                        auto ref_bases = new ALLOC_FLEXIBLE_IN_POOL(_pool, Bases, end - start + 1, uint8_t) Bases{end - start + 1};
                        auto alt_bases = new ALLOC_FLEXIBLE_IN_POOL(_pool, Bases, end - start + 1, uint8_t) Bases{end - start + 1};
                        memcpy(ref_bases->data, ref->data + ref_pos + start, (end - start + 1) * sizeof(uint8_t));
                        memcpy(alt_bases->data, alignment->data + alignment_pos + start, (end - start + 1) * sizeof(uint8_t));

                        ref_allele = Allele::create_allele(ref_bases, true, _pool);
                        alt_allele = Allele::create_allele(alt_bases, false, _pool);
                    }

                    alleles.emplace_back(ref_allele);
                    alleles.emplace_back(alt_allele);

                    mismatch_vc = Variant::create(_pool);
                    mismatch_vc->set_source_id(_source_index);
                    mismatch_vc->set_tid(tid);
                    mismatch_vc->set_start(ref_loc->get_start() + ref_pos + start);
                    mismatch_vc->set_stop(ref_loc->get_start() + ref_pos + end);
                    mismatch_vc->set_alleles(alleles);
                    proposed_events.emplace_back(mismatch_vc);
                }

                ref_pos += op_len;
                alignment_pos += op_len;
                break;
            }
            case BAM_CREF_SKIP:
            case BAM_CHARD_CLIP:
            case BAM_CPAD:
            default: {
                RovacaLogger::error("unsupported cigar operator created during SW alignment: {}", op);
                exit(EXIT_FAILURE);
            }
        }
    }

    for (pVariant vc : proposed_events) {
        add_variant(vc);
    }
}

void EventMap::add_variant(pVariant vc)
{
    int64_t start_key = vc->get_start();
    if (_events.count(vc->get_start())) {
        pVariant pre = _events.at(start_key);
        _events[start_key] = make_block(pre, vc);
    }
    else {
        _events.insert({vc->get_start(), vc});
    }
}

}  // namespace rovaca