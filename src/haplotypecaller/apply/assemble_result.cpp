#include "assemble_result.h"

#include <algorithm>
#include <iostream>

#include "cigar_builder.h"
#include "haplotype.h"
#include "rovaca_logger.h"
#include "read_record.h"
#include "simple_interval.h"
#include "utils/base_utils.h"
#include "utils/cigar_utils.h"

AssembleRegion::AssembleRegion(p_hc_region_active_storage region, pMemoryPool mem)
    : isActive(region->active != 0)
    , activeSpan(SimpleInterval::create(region->tid, region->activeSpan.start, region->activeSpan.end, mem))
    , paddedSpan(SimpleInterval::create(region->tid, region->paddedSpan.start, region->paddedSpan.end, mem))
    , resource(mem)
{}

AssembleRegion::AssembleRegion(bool active, rovaca::pSimpleInterval original, rovaca::pSimpleInterval padded)
    : isActive(active)
    , activeSpan(original)
    , paddedSpan(padded)
{}

AssembleRegion::AssembleRegion(const AssembleRegion& one, pMemoryPool mem)
    : isActive(one.isActive)
    , activeSpan(SimpleInterval::create(one.activeSpan->get_tid(), one.activeSpan->get_start(), one.activeSpan->get_stop(), mem))
    , paddedSpan(SimpleInterval::create(one.paddedSpan->get_tid(), one.paddedSpan->get_start(), one.paddedSpan->get_stop(), mem))
    , resource(mem)
{}

AssembleResult::AssembleResult(pMemoryPool resource)
    : resource_(resource)
    , padded_ref_()
    , refHaplotype_()
    , padded_ref_loc_()
    , genotyping_region_()
    , reads_(resource)
    , hapltypes_(resource)
    , variant_events_()
    , kemers_()
{}

AssembleResult::AssembleResult(const AssembleResult& one)
    : resource_(std::pmr::new_delete_resource())
    , padded_ref_()
    , refHaplotype_()
    , padded_ref_loc_()
    , genotyping_region_()
    , reads_(resource_)
    , hapltypes_(resource_)
    , variant_events_()
    , kemers_()
{
    // Copy only necessary fields.
    AssembleResult& one_ = const_cast<AssembleResult&>(one);
    genotyping_region_ = new (resource_->allocate(sizeof(AssembleRegion))) AssembleRegion(*one_.get_region(), resource_);
    const auto& reads = one_.get_reads();
    const auto& haplotypes = one_.get_haplotypes();
    size_t reads_count = reads.size();
    size_t hap_count = haplotypes.size();
    reads_.reserve(reads_count);
    hapltypes_.reserve(hap_count);
    for (size_t i = 0; i < reads_count; ++i) {
        // Do deep copy.
        p_hc_apply_one_read assemble_r = new (resource_->allocate(sizeof(hc_apply_one_read))) hc_apply_one_read;

        assemble_r->pos_start = reads[i]->assemble_read()->pos_start;
        assemble_r->pos_end = reads[i]->assemble_read()->pos_end;
        assemble_r->cigar_len = reads[i]->assemble_read()->cigar_len;
        assemble_r->read_len = reads[i]->assemble_read()->read_len;
        assemble_r->ref_len = reads[i]->assemble_read()->ref_len;
        assemble_r->insert_size = reads[i]->assemble_read()->insert_size;

        assemble_r->read_data = new (resource_->allocate(sizeof(uint8_t) * reads[i]->assemble_read()->read.l_data)) uint8_t();
        // assemble_r->cigar_str = new (resource_->allocate(sizeof(char) * strlen(reads[i]->assemble_read()->cigar_str))) char();
        assemble_r->seq = new (resource_->allocate(sizeof(uint8_t) * assemble_r->read_len)) uint8_t();
        assemble_r->qual = new (resource_->allocate(sizeof(uint8_t) * assemble_r->read_len)) uint8_t();
        assemble_r->cigar = reads[i]->assemble_read()->cigar;
        // new (resource_->allocate(sizeof(uint32_t))) uint32_t();
        memcpy(assemble_r->seq, reads[i]->assemble_read()->seq, sizeof(uint8_t) * reads[i]->assemble_read()->read_len);
        memcpy(assemble_r->qual, reads[i]->assemble_read()->qual, sizeof(uint8_t) * reads[i]->assemble_read()->read_len);
        // memcpy(assemble_r->cigar, reads[i]->assemble_read()->cigar, sizeof(uint32_t));
        // memcpy(assemble_r->cigar_str, reads[i]->assemble_read()->cigar_str, sizeof(char) * strlen(reads[i]->assemble_read()->cigar_str));
        memcpy(assemble_r->read_data, reads[i]->assemble_read()->read_data, sizeof(uint8_t) * reads[i]->assemble_read()->read.l_data);
        pReadRecord r = ReadRecord::create(resource_, reads[i]->header(), assemble_r);
        reads_.emplace_back(r);
    }
    for (size_t i = 0; i < hap_count; ++i) {
        pHaplotype h = Haplotype::create(resource_);
        auto base = BaseUtils::bases_create((char*)(haplotypes[i]->get_bases()->data), resource_);
        h->init_haplotype(base, haplotypes[i]->is_reference());
        hapltypes_.emplace_back(h);
    }
}
AssembleResult::~AssembleResult() { resource_->~memory_resource(); }
AssembleResult* AssembleResult::create(pMemoryPool resource)
{
    return new (resource->allocate(sizeof(AssembleResult))) AssembleResult{resource};
}

AssembleResult::AssembleResult(pAssembleRegion assemble_region, uint8_t* ref_with_padded, pSimpleInterval refloc_with_padded,
                               const std::pmr::vector<pReadRecord>& assemble_reads, const std::pmr::vector<pHaplotype>& assemble_hapltypes)
    : padded_ref_(ref_with_padded)
    , padded_ref_loc_(refloc_with_padded)
    , genotyping_region_(assemble_region)
    , reads_(assemble_reads)
    , hapltypes_(assemble_hapltypes)
    , variant_events_({})
{}

AssembleResult* AssembleResult::trim_to(pAssembleRegion trimmedAssemblyRegion, const AssembleArgument& arguments, pMemoryPool resource)
{
    pSimpleInterval padded_span = trimmedAssemblyRegion->paddedSpan;
    HaplotypeHashMap originalByTrimmedHaplotypes = calculateOriginalByTrimmedHaplotypes(padded_span, arguments, resource);
    if (!refHaplotype_) {
        throw std::runtime_error("refHaplotype is null");
    }

    AssembleResult* result = nullptr;
    for (const auto& trimmedPair : originalByTrimmedHaplotypes) {
        pHaplotype trimmed = trimmedPair.first;
        pHaplotype original = trimmedPair.second;
        if (!original) {
            throw std::runtime_error("all trimmed haplotypes must have an original one");
        }
        result->add_result(trimmed);
    }
    result->set_region(trimmedAssemblyRegion);
    result->set_padded_ref(padded_ref_);
    result->set_padded_refloc(padded_ref_loc_);
    if (!result->refHaplotype()) {
        throw std::runtime_error("missing reference haplotype in the trimmed set");
    }
    return result;
}

HaplotypeHashMap AssembleResult::calculateOriginalByTrimmedHaplotypes(pSimpleInterval span, const AssembleArgument& arguments,
                                                                      pMemoryPool resource)
{
    if (arguments.debugAssembly) {}
    std::pmr::vector<pHaplotype> haplotypelist = get_haplotypes();
    HaplotypeHashMap originalByTrimmedHaplotypes = trimDownHaplotypes(span, haplotypelist, resource);
    std::pmr::vector<pHaplotype> trimmedHaplotypes;
    trimmedHaplotypes.reserve(originalByTrimmedHaplotypes.size());
    for (const auto haps : originalByTrimmedHaplotypes) {
        trimmedHaplotypes.emplace_back(haps.first);
    }
    // 按照长度以及字典序排序.
    if (trimmedHaplotypes.size() > 1) {
        std::sort(trimmedHaplotypes.begin(), trimmedHaplotypes.end(), [](const pHaplotype& hap_a, const pHaplotype& hap_b) {
            return hap_a->get_bases()->num < hap_b->get_bases()->num ||
                   strncmp((char*)hap_a->get_bases()->data, (char*)hap_b->get_bases()->data, hap_a->get_bases()->num) < 0;
        });
    }

    HaplotypeHashMap sortedOriginalByTrimmedHaplotypes;
    mapOriginalToTrimmed(originalByTrimmedHaplotypes, sortedOriginalByTrimmedHaplotypes, trimmedHaplotypes);
    if (arguments.debugAssembly) {}
    return sortedOriginalByTrimmedHaplotypes;
}

void AssembleResult::mapOriginalToTrimmed(const HaplotypeHashMap& originalByTrimmedHaplotypes,
                                          HaplotypeHashMap& sortedOriginalByTrimmedHaplotypes,
                                          const std::pmr::vector<pHaplotype>& trimmedHaplotypes)
{
    for (pHaplotype trimmed : trimmedHaplotypes) {
        sortedOriginalByTrimmedHaplotypes.insert({trimmed, originalByTrimmedHaplotypes.at(trimmed)});
    }
}

HaplotypeHashMap AssembleResult::trimDownHaplotypes(pSimpleInterval span, const std::pmr::vector<pHaplotype>& haplotypeList,
                                                    pMemoryPool resource)
{
    HaplotypeHashMap originalByTrimmedHaplotypes;
    for (const pHaplotype hap : haplotypeList) {
        pHaplotype trimmed = trim_haplotype(hap, span, true, resource);
        if (trimmed != nullptr) {
            HaplotypeHashMapIterator trimmed_itr = originalByTrimmedHaplotypes.find(trimmed);
            if (trimmed_itr != originalByTrimmedHaplotypes.end()) {
                if (hap->is_reference()) {
                    originalByTrimmedHaplotypes.erase(trimmed_itr);
                    originalByTrimmedHaplotypes.insert({trimmed, hap});
                }
            }
            else {
                originalByTrimmedHaplotypes.insert({trimmed, hap});
            }
        }
        else if (hap->is_reference()) {
            throw std::runtime_error("trimming eliminates the reference haplotype");
        }
    }
    // Now set reference status originalByTrimmedHaplotypes.
    HaplotypeHashMap fixedOriginalByTrimmedHaplotypes;
    for (const auto& hap_pair : originalByTrimmedHaplotypes) {
        pHaplotype h = hap_pair.first;
        pHaplotype hPair = hap_pair.second;
        if (hPair->is_reference()) {
            pHaplotype fixedHap = Haplotype::create(resource);
            fixedHap->set_cigar(h->cigar());
            fixedHap->set_interval(h->interval());
            fixedHap->set_score(h->score());
            fixedHap->set_alignment_start_hap_wrt_ref(h->alignment_start_hap_wrt_ref());
            fixedOriginalByTrimmedHaplotypes.insert({fixedHap, hPair});
        }
        fixedOriginalByTrimmedHaplotypes.insert(hap_pair);
    }
    return fixedOriginalByTrimmedHaplotypes;
}

bool AssembleResult::add_result(pHaplotype haplotype)
{
    if (haplotype == nullptr) {
        std::cerr << "input haplotype cannot be null" << std::endl;
    }
    if (haplotype->interval()) {
        std::cerr << "haplotype genomeLocation cannot be null" << std::endl;
    }
    for (const auto hap : hapltypes_) {
        if (hap->equals(*haplotype)) {
            return false;
        }
    }
    hapltypes_.emplace_back(haplotype);
    updateReferenceHaplotype(haplotype);
    return true;
}

void AssembleResult::updateReferenceHaplotype(pHaplotype newHaplotype)
{
    if (!newHaplotype->is_reference()) {
        return;
    }
    if (refHaplotype() == nullptr) {
        refHaplotype_ = newHaplotype;
    }
    else {  // assumes that we have checked wether the haplotype is already in the collection and so is no need to check equality.
        throw std::runtime_error("the assembly-result-set already have a reference haplotype that is different");
    }
}

pHaplotype AssembleResult::trim_haplotype(pHaplotype original, pSimpleInterval loc, bool ignoreRefState, pMemoryPool pool)
{
    if (loc == nullptr) {
        return nullptr;
    }
    if (original->interval() == nullptr) {
        return nullptr;
    }
    if (original->cigar() == nullptr) {
        return nullptr;
    }

    int newStart = loc->get_start() - original->interval()->get_start();
    int newStop = newStart + loc->get_stop() - loc->get_start();
    pBases newBases = get_bases_covering_ref_interval(newStart, newStop, original->get_bases(), 0, original->cigar(), pool);
    if (newBases == nullptr || newBases->num == 0) {
        return nullptr;
    }
    pCigar newCigar = trim_cigar_by_reference(original->cigar(), newStart, newStop, pool);
    bool leadingInsertion = !(bam_cigar_type(newCigar->data[0]) & 0x2);
    bool trailingInsertion = !(bam_cigar_type(newCigar->data[newCigar->num - 1]) & 0x2);
    int firstIndexToKeepInclusive = leadingInsertion ? 1 : 0;
    int lastIndexToKeepExclusive = newCigar->num - (trailingInsertion ? 1 : 0);

    if (lastIndexToKeepExclusive <= firstIndexToKeepInclusive) {  // edge case of entire cigar is insertion
        return nullptr;
    }
    pCigar leadingIndelTrimmedNewCigar =
        !(leadingInsertion || trailingInsertion)
            ? newCigar
            : CigarBuilder::create(pool)
                  ->add_all(newCigar->data + firstIndexToKeepInclusive, lastIndexToKeepExclusive - firstIndexToKeepInclusive)
                  ->make();

    pHaplotype ret = Haplotype::create(pool);
    ret->init_haplotype((const char*)(newBases->data), ignoreRefState ? false : original->is_reference(), pool);
    ret->set_cigar(leadingIndelTrimmedNewCigar);
    ret->set_interval(loc);
    ret->set_score(original->score());
    ret->set_kmer_size(original->kmer_size());
    ret->set_alignment_start_hap_wrt_ref(newStart + original->alignment_start_hap_wrt_ref());
    return ret;
}

pBases AssembleResult::get_bases_covering_ref_interval(int64_t ref_start, int64_t ref_end, pBases bases, int64_t bases_start_on_ref,
                                                       pCigar bases_to_ref_cigar, pMemoryPool pool)
{
    CHECK_CONDITION_EXIT(ref_start < 0 || ref_end < ref_start, "bad start {} and/or stop {}", ref_start, ref_end);
    CHECK_CONDITION_EXIT(bases_start_on_ref < 0, "bases_start_on_ref must be >= 0");
    CHECK_CONDITION_EXIT(nullptr == bases || nullptr == bases_to_ref_cigar, "nullptr == bases || nullptr == bases_to_ref_cigar");
    CHECK_CONDITION_EXIT(bases->num != bam_cigar2rlen((int32_t)bases_to_ref_cigar->num, bases_to_ref_cigar->data),
                         "mismatch in length between reference and cigar");

    int64_t ref_pos = bases_start_on_ref;
    int64_t bases_pos = 0;
    int64_t bases_start = INVALID_INT;
    int64_t bases_stop = INVALID_INT;
    bool done = false;

    uint32_t ce, op, op_len;
    for (uint32_t iii = 0; !done && iii < bases_to_ref_cigar->num; iii++) {
        ce = bases_to_ref_cigar->data[iii];
        op = bam_cigar_op(ce);
        op_len = bam_cigar_oplen(ce);
        switch (op) {
            case BAM_CINS: {
                bases_pos += op_len;
                break;
            }
            case BAM_CDEL: {
                for (uint32_t i = 0; i < op_len; i++) {
                    if (ref_pos == ref_end || ref_pos == ref_start) {
                        return nullptr;
                    }
                    ref_pos++;
                }
                break;
            }
            case BAM_CMATCH:
            case BAM_CEQUAL:
            case BAM_CDIFF: {
                for (uint32_t i = 0; i < op_len; i++) {
                    if (ref_pos == ref_start) bases_start = bases_pos;
                    if (ref_pos == ref_end) {
                        bases_stop = bases_pos;
                        done = true;
                        break;
                    }
                    ref_pos++;
                    bases_pos++;
                }
                break;
            }
            case BAM_CREF_SKIP:
            case BAM_CSOFT_CLIP:
            case BAM_CHARD_CLIP:
            case BAM_CPAD:
            case BAM_CBACK:
            default: {
                throw std::runtime_error("unsupported operator");
            }
        }
    }

    CHECK_CONDITION_EXIT(bases_start == INVALID_INT || bases_stop == INVALID_INT, "never found start or stop");
    char* start_ptr = (char*)bases->data + bases_start;
    uint32_t len = bases_stop - bases_start + 1;
    return BaseUtils::bases_create(start_ptr, len, pool);
}

pCigar AssembleResult::trim_cigar_by_reference(pCigar cigar, int64_t start, int64_t end, pMemoryPool pool)
{
    CHECK_CONDITION_EXIT(start < 0, "start position can't be negative");
    CHECK_CONDITION_EXIT(end < start, "end is smaller than start");

    pCigarBuilder builder = CigarBuilder::create(pool);

    // these variables track the inclusive start and exclusive end of the current cigar element in reference (if byReference) or read
    // (otherwise) coordinates
    int64_t element_start;    // inclusive
    int64_t element_end = 0;  // exclusive -- start of next element
    int64_t overlap_length;
    uint32_t element;
    for (uint32_t i = 0; i < cigar->num; i++) {
        element = cigar->data[i];
        element_start = element_end;
        element_end = element_start + consumes_ref_bases(bam_cigar_op(element)) ? bam_cigar_oplen(element) : 0;

        // we are careful to include zero-length elements at both ends, that is, elements with elementStart == elementEnd == start and
        // elementStart == elementEnd == end + 1
        if (element_end < start || (element_end == start && element_start < start)) {
            continue;
        }
        else if (element_start > end && element_end > end + 1) {
            break;
        }

        overlap_length =
            element_end == element_start ? bam_cigar_oplen(element) : std::min(end + 1, element_end) - std::max(start, element_start);

        builder->add(overlap_length << BAM_CIGAR_SHIFT | bam_cigar_op(element));
    }
    CHECK_CONDITION_EXIT(element_end <= end, "cigar elements don't reach end position (inclusive) {}", end);
    return builder->make_and_record_deletions_removed_result().cigar;
}
