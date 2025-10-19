#include "adapter_utils.h"

#include <algorithm>

#include "constants_str.hpp"
#include "rovaca_logger.h"
#include "event_map.h"
#include "genotype_argument.h"
#include "genotypes_context.hpp"
#include "haplotype.h"
#include "info_data.hpp"
#include "read_clipper.h"
#include "read_record.h"
#include "simple_interval.h"
#include "utils/debug_utils.h"
#include "utils/haplotype_utils.h"
#include "utils/rovaca_variant_context_utils.h"
#include "variant.h"

namespace rovaca
{

#define READ_LENGTH_FILTER_THRESHOLD (10)

static const int32_t k_hom_gt[HOM_GT_CUNT]{2, 2};

static pConstantsStr k_constan = ConstantsStr::get_instance();

static const std::string k_a_non_ref_str{"A,<NON_REF>"};
static const std::string k_t_non_ref_str{"T,<NON_REF>"};
static const std::string k_g_non_ref_str{"G,<NON_REF>"};
static const std::string k_c_non_ref_str{"C,<NON_REF>"};
static const std::string k_n_non_ref_str{"N,<NON_REF>"};

const char* base2allele_str(uint8_t b)
{
    switch (b) {
        case 'A': return k_a_non_ref_str.c_str();
        case 'T': return k_t_non_ref_str.c_str();
        case 'G': return k_g_non_ref_str.c_str();
        case 'C': return k_c_non_ref_str.c_str();
        case 'N': return k_n_non_ref_str.c_str();
        default: {
            RovacaLogger::error("invalid bases: {}", b);
            exit(EXIT_FAILURE);
        }
    }
}

struct HaplotypeVariantContextComparator
{
    bool operator()(pVariant v1, pVariant v2) const
    {
        if (v1->get_start() != v2->get_start()) {
            return v1->get_start() < v2->get_start();
        }
        if (v1->ref_allele()->length() != v2->ref_allele()->length()) {
            return v1->ref_allele()->length() < v2->ref_allele()->length();
        }
        return v1->alternate_allele_at(0) < v2->alternate_allele_at(0);
    }
};

IntervalPair AdapterUtils::non_variant_left_flank_region(pSimpleInterval original, pSimpleInterval original_padded, pSimpleInterval variant,
                                                         int64_t padding, int64_t chr_len, pMemoryPool pool)
{
    pSimpleInterval new_span{nullptr}, new_span_padded{nullptr};
    if (original->get_start() < variant->get_start()) {
        pSimpleInterval flank = SimpleInterval::create(original->get_tid(), original->get_start(), variant->get_start() - 1, pool);
        pSimpleInterval flank_padded = flank->expand_within_contig(padding, chr_len, pool);
        CHECK_CONDITION_EXIT(!flank_padded->contains(*flank), "the left_flank_padded must fully contain the left_flank");
        new_span = original->intersect(*flank, pool);
        new_span_padded = original_padded->intersect(*flank_padded, pool);
    }
    return {new_span, new_span_padded};
}

IntervalPair AdapterUtils::non_variant_right_flank_region(pSimpleInterval original, pSimpleInterval original_padded,
                                                          pSimpleInterval variant, int64_t padding, int64_t chr_len, pMemoryPool pool)
{
    pSimpleInterval new_span{nullptr}, new_span_padded{nullptr};
    if (variant->get_stop() < original->get_stop()) {
        pSimpleInterval flank = SimpleInterval::create(original->get_tid(), variant->get_stop() + 1, original->get_stop(), pool);
        pSimpleInterval flank_padded = flank->expand_within_contig(padding, chr_len, pool);
        CHECK_CONDITION_EXIT(!flank_padded->contains(*flank), "the left_flank_padded must fully contain the left_flank");
        new_span = original->intersect(*flank, pool);
        new_span_padded = original_padded->intersect(*flank_padded, pool);
    }
    return {new_span, new_span_padded};
}

IntervalPair AdapterUtils::trim_region(HaplotypeVector& haplotypes, pRefFragment ref, pSimpleInterval ref_loc, pSimpleInterval original,
                                       pSimpleInterval original_padded, pHCArgs args, pMemoryPool pool)
{
    EventMap::build_event_maps_for_haplotypes(haplotypes, ref, ref_loc, args->max_mnp_distance, pool);

    std::pmr::set<pVariant, HaplotypeVariantContextComparator> sorted_variant{pool};
    for (pHaplotype h : haplotypes) {
        const auto& events = h->event_map()->get_events();
        for (const auto& tup : events) {
            sorted_variant.insert(tup.second);
        }
    }

    std::pmr::list<pVariant> variants_in_region{pool};
    std::for_each(sorted_variant.begin(), sorted_variant.end(), [&](pVariant vc) {
        if (original->overlaps(*vc)) {
            variants_in_region.push_back(vc);
        }
    });
    if (variants_in_region.empty()) {
        return {nullptr, nullptr};
    }

    int64_t min_start = INT32_MAX, max_end = INT32_MIN;
    std::for_each(variants_in_region.begin(), variants_in_region.end(), [&](pVariant vc) {
        min_start = std::min(min_start, vc->get_start());
        max_end = std::max(max_end, vc->get_stop());
    });

    pSimpleInterval variant_span = SimpleInterval::create(original->get_tid(), min_start, max_end, pool)->intersect(*original, pool);

    bool ret;
    int64_t padding, start_index;
    RepeatUnitsResult num_repeats_and_unit{pool};
    RefFragment ref_bases_starting_at_vc_without_pad;
    for (pVariant vc : variants_in_region) {
        padding = args->snp_padding_for_genotyping;
        if (vc->is_indel1()) {
            num_repeats_and_unit.lengths.clear();
            padding = args->indel_padding_for_genotyping;
            start_index = vc->get_start() + 1 - ref_loc->get_start();

            ref_bases_starting_at_vc_without_pad.len = ref_loc->get_stop() - vc->get_start();
            ref_bases_starting_at_vc_without_pad.data = ref->data + start_index;

            ret = ROVACAVariantContextUtils::get_num_tandem_repeat_units(vc, &ref_bases_starting_at_vc_without_pad, &num_repeats_and_unit,
                                                                       pool);

            if (ret && num_repeats_and_unit.bases.data && num_repeats_and_unit.bases.len) {
                auto repeat_length = int32_t(num_repeats_and_unit.bases.len);
                int32_t most_repeats = *std::max_element(num_repeats_and_unit.lengths.begin(), num_repeats_and_unit.lengths.end());
                int32_t longest_str = most_repeats * repeat_length;
                padding = args->str_padding_for_genotyping + longest_str;
            }
        }

        min_start = std::min(min_start, std::max(vc->get_start() - padding, (int64_t)1));
        max_end = std::max(max_end, vc->get_stop() + padding);
    }

    pSimpleInterval padded_span = SimpleInterval::create(original->get_tid(), min_start, max_end, pool)->intersect(*original_padded, pool);

    return {variant_span, padded_span};
}

ReadHashSet AdapterUtils::trim_reads_by_region(const ReadHashSet& reads, pSimpleInterval padded, pMemoryPool pool, pBamDataPool bpool,
                                               std::pmr::list<bam1_t*>& extra_memory_reads)
{
    std::pmr::set<pReadRecord> sort_reads{pool};

    pReadRecord new_r{};
    int64_t start = padded->get_start(), stop = padded->get_stop();
    for (pReadRecord r : reads) {
        new_r = ReadClipper{r, pool, bpool}.hard_clip_to_region(start, stop);
        if (new_r && !new_r->is_empty() && new_r->overlaps(*padded)) {
            sort_reads.insert(new_r);
            if (new_r != r && new_r->raw_data() != r->raw_data() && (bam_get_mempolicy(new_r->raw_data()) & BAM_USER_OWNS_DATA) == 0) {
                extra_memory_reads.push_back(new_r->raw_data());
            }
        }
    }
    return {{sort_reads.begin(), sort_reads.end()}, pool};
}

static inline bool mate_on_same_contig_or_no_mapped_mate(pReadRecord read)
{
    return !read->is_paired() || read->mate_is_unmapped() || (!read->is_unmapped() && read->get_tid() == read->mate_tid());
}

struct HaplotypeLess
{
    bool operator()(const pHaplotype& l, const pHaplotype& r) const
    {
        rovaca::pBases l_bases = l->get_display_string();
        rovaca::pBases r_bases = r->get_display_string();
        return rovaca::less_bases(l_bases, r_bases);
    }
};

HaplotypeVector AdapterUtils::trim_haplotype_by_region(const HaplotypeVector& haplotypes, pSimpleInterval padded, pMemoryPool pool)
{
    if (haplotypes.empty()) {
        return {};
    }

    size_t ref_index = 0;
    while (!haplotypes.at(ref_index)->is_reference()) {
        ++ref_index;
    }

    pHaplotype ref_hap = haplotypes.at(ref_index);
    pHaplotype trimed_ref_hap = HaplotypeUtils::trim(ref_hap, padded, pool);

    std::pmr::set<pHaplotype, HaplotypeLess> uniqie_haplotypes{pool};
    uniqie_haplotypes.insert(trimed_ref_hap);

    for (size_t i = 0, hap_size = haplotypes.size(); i < hap_size; ++i) {
        if (i == ref_index) continue;

        pHaplotype trimed_h = HaplotypeUtils::trim(haplotypes.at(i), padded, true, pool);
        if (!trimed_h) continue;

        if (!uniqie_haplotypes.count(trimed_h)) {
            uniqie_haplotypes.insert(trimed_h);
        }
    }

    HaplotypeVector result{pool};
    result.reserve(uniqie_haplotypes.size());
    result.push_back(trimed_ref_hap);
    std::copy_if(uniqie_haplotypes.begin(), uniqie_haplotypes.end(), std::back_inserter(result),
                 [](pHaplotype h) { return !h->is_reference(); });

    return result;
}

ReadList AdapterUtils::filter_non_passing_reads1(ReadHashSet& reads, int32_t minimum_read_length_after_trimming, pMemoryPool pool)
{
    ReadList filter_reads{pool};
    for (pReadRecord read : reads) {
        if (read->unclipped_read_length() < minimum_read_length_after_trimming) {
            filter_reads.push_back(read);
            reads.erase(read);
        }
    }
    return filter_reads;
}

ReadList AdapterUtils::filter_non_passing_reads2(ReadHashSet& reads, uint8_t mapping_quality_threshold, pMemoryPool pool)
{
    ReadList filter_reads{pool};
    for (pReadRecord read : reads) {
        if (read->unclipped_read_length() < READ_LENGTH_FILTER_THRESHOLD || read->mapping_quality() < mapping_quality_threshold ||
            !mate_on_same_contig_or_no_mapped_mate(read)) {
            filter_reads.push_back(read);
            reads.erase(read);
        }
    }
    return filter_reads;
}

bcf1_t* AdapterUtils::variant2bcf(bam_hdr_t* bam_header, bcf_hdr_t* vcf_header, pVariant vc, bcf1_t* rec)
{
    int ret;
    bcf_clear1(rec);

    // .. CHROM
    rec->rid = bcf_hdr_name2id(vcf_header, bam_header->target_name[vc->get_tid()]);
    // .. POS
    rec->pos = vc->get_start() - 1;  // 0 bases
    // .. ID
    ret = bcf_update_id(vcf_header, rec, vc->db_id() ? (char*)vc->db_id()->data : nullptr);
    CHECK_CONDITION_EXIT(ret != 0, "bcf_update_id");
    // .. REF and ALT
    const AlleleVector& alleles = vc->alleles();
    std::pmr::vector<uint8_t*> alleles_str(alleles.get_allocator());
    alleles_str.reserve(alleles.size());
    std::for_each(alleles.begin(), alleles.end(), [&](pAllele a) { alleles_str.push_back(a->get_display_string()->data); });
    ret = bcf_update_alleles(vcf_header, rec, (const char**)alleles_str.data(), int32_t(alleles.size()));
    CHECK_CONDITION_EXIT(ret != 0, "bcf_update_id");
    // .. QUAL
    if (ROVACA_LIKELY(vc->has_log10_error())) {
        rec->qual = static_cast<float>(vc->get_phred_scaled_qual());
    }
    else {
        bcf_float_set_missing(rec->qual);
    }
    // .. INFO
    vc->info()->info2bcf(vcf_header, rec);
    // .. FORMAT
    pGenotypesContext gc = vc->genotype();
    if (gc) {
        pGenotype g = !gc->empty() ? gc->at(0) : nullptr;
        if (g) {
            std::pmr::map<pAllele, int32_t> mm = DebugUtils::bcf_build_allele_map(vc->alleles());
            g->genotype2bcf(vcf_header, rec, mm);
        }
    }

    return rec;
}

bcf1_t* AdapterUtils::variant2bcf(bam_hdr_t* bam_header, bcf_hdr_t* vcf_header, int32_t tid, int64_t pos, pAllele ref, int32_t block_end,
                                  int32_t dp, int32_t gq, int32_t min_dp, const std::vector<int32_t>& pls, bcf1_t* rec)
{
    int ret;
    bcf_clear1(rec);

    // CHROM
    rec->rid = bcf_hdr_name2id(vcf_header, bam_header->target_name[tid]);
    // .. POS
    rec->pos = pos - 1;  // 0 bases
    // .. ID
    bcf_update_id(vcf_header, rec, nullptr);
    // .. REF and ALT
    const char* allele_str = base2allele_str(ref->get_bases()->data[0]);
    ret = bcf_update_alleles_str(vcf_header, rec, allele_str);
    CHECK_CONDITION_EXIT(ret != 0, "bcf_update_alleles_str");
    // .. QUAL
    bcf_float_set_missing(rec->qual);
    // .. FILTER
    bcf_update_filter(vcf_header, rec, nullptr, 0);
    // .. INFO
    ret = bcf_update_info_int32(vcf_header, rec, k_constan->k_end_key.c_str(), &block_end, 1);
    CHECK_CONDITION_EXIT(ret != 0, "bcf_update_info_int32: {}", k_constan->k_end_key.c_str());
    // .. FORMAT
    ret = bcf_update_genotypes(vcf_header, rec, k_hom_gt, bcf_hdr_nsamples(vcf_header) * 2);
    CHECK_CONDITION_EXIT(ret != 0, "bcf_update_genotypes");

    ret = bcf_update_format_int32(vcf_header, rec, k_constan->k_depth_key.c_str(), &dp, 1);
    CHECK_CONDITION_EXIT(ret != 0, "bcf_update_info_int32: {}", k_constan->k_depth_key.c_str());

    int32_t filter_gq = std::min(gq, MAX_GENOTYPE_QUAL);
    ret = bcf_update_format_int32(vcf_header, rec, k_constan->k_genotype_quality_key.c_str(), &filter_gq, 1);
    CHECK_CONDITION_EXIT(ret != 0, "bcf_update_info_int32: {}", k_constan->k_genotype_quality_key.c_str());

    ret = bcf_update_format_int32(vcf_header, rec, k_constan->k_min_dp_key.c_str(), &min_dp, 1);
    CHECK_CONDITION_EXIT(ret != 0, "bcf_update_info_int32: {}", k_constan->k_min_dp_key.c_str());

    ret = bcf_update_format_int32(vcf_header, rec, k_constan->k_genotype_pl_key.c_str(), pls.data(), pls.size());
    CHECK_CONDITION_EXIT(ret != 0, "bcf_update_info_int32: {}", k_constan->k_genotype_pl_key.c_str());

    return rec;
}

}  // namespace rovaca
