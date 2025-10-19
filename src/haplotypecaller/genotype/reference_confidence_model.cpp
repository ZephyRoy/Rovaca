#include "reference_confidence_model.h"

#include <algorithm>

#include "allele_likelihoods.hpp"
#include "rovaca_logger.h"
#include "genotype_likelihoods.h"
#include "genotype_likelihoods_cache.h"
#include "genotypes_context.hpp"
#include "haplotype.h"
#include "math_utils.h"
#include "pileup_element.h"
#include "quality_utils.h"
#include "read_pileup.h"
#include "read_record.h"
#include "simple_interval.h"
#include "utils/alignment_utils.h"
#include "utils/assembly_based_caller_utils.h"
#include "utils/cigar_utils.h"
#include "utils/rovaca_variant_context_utils.h"
#include "variant.h"

namespace rovaca
{

#define IDX_HOM_REF                   (0)
#define MAX_N_INDEL_INFORMATIVE_READS (40)
#define MAX_INDEL_SIZE                (10)

static constexpr uint8_t m_mask[256]{
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 2, 0, 0, 0, 4, 0, 0, 0, 0, 0, 0, 15, 0, 0, 0, 0, 0, 8, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  0, 0, 0, 0, 0};

/*!
 * @brief
 * @param q: read_byte
 * @param r: ref_byte
 */
#define IS_MISMATCH_AND_NOT_AN_ALIGNMENT_GAP(q, r) ((0 == ((m_mask[q]) & (m_mask[r]))) && ((q) != AlignmentUtils::s_gap_base_character))

pGenotypeLikelihoodsCache ReferenceConfidenceModel::s_cache = GenotypeLikelihoodsCache::get_instance();

VariantVector ReferenceConfidenceModel::calculate_ref_confidence(const VariantVector& calls, pHaplotype ref_h, pSimpleInterval flank,
                                                                 pSimpleInterval flank_padded, int32_t ploidy, const ReadHashSet& reads)
{
    CHECK_CONDITION_EXIT(!ref_h || !flank || !flank_padded, "nullptr");
    CHECK_CONDITION_EXIT(int64_t(ref_h->length()) != flank_padded->size(), "ref_haplotype->length() != pad_active_loc->size()");

    VariantVector results{_pool};
    results.reserve(flank_padded->size());

    int32_t sample_idx = 0;
    pBases ref_bases = ref_h->get_bases();
    int64_t global_ref_offset = flank->get_start() - flank_padded->get_start();

    ReadPileupVector ref_pileups = AssemblyBasedCallerUtils::get_pileups_over_reference(reads, flank, ref_bases, global_ref_offset, _pool);

    int64_t offset;
    pVariant overlapping_site;
    for (pReadPileup pileup : ref_pileups) {
        offset = pileup->get_start() - flank->get_start();
        overlapping_site = ROVACAVariantContextUtils::get_overlapping_variant_context(pileup, calls);
        if (nullptr != overlapping_site && overlapping_site->get_start() == pileup->get_start()) {
            results.push_back(overlapping_site);
        }
        else {
            results.push_back(make_reference_confidence_variant_context(ploidy, ref_bases, sample_idx, global_ref_offset, pileup, offset));
        }
    }

    std::for_each(reads.begin(), reads.end(), [](pReadRecord r) { r->clear_informative_set(); });
    return results;
}

pVariant ReferenceConfidenceModel::make_reference_confidence_variant_context(int32_t ploidy, pBases ref, int32_t sample_idx,
                                                                             int64_t global_ref_offset, pReadPileup pileup, int64_t offset)
{
    int64_t ref_offset = offset + global_ref_offset;
    uint8_t ref_base = ref->data[ref_offset];

    pRefVsAnyResult hom_ref_calc = calc_genotype_likelihoods_of_ref_vs_any(ploidy, pileup);

    do_indel_ref_conf_calc(ploidy, ref, pileup, ref_offset, hom_ref_calc);

    pAllele ref_allele = Allele::create_allele(ref_base, 1);

    AlleleVector genotype_alleles = ROVACAVariantContextUtils::homozygous_allele_list(ref_allele, ploidy, _pool);
    pGenotype g = Genotype::create(_pool);
    g->set_id(sample_idx);
    g->set_ad({{hom_ref_calc->ref_depth, hom_ref_calc->non_ref_depth}, _pool});
    g->set_dp(hom_ref_calc->ref_depth + hom_ref_calc->non_ref_depth);
    g->set_alleles(std::move(genotype_alleles));
    g->set_gq(ROVACAVariantContextUtils::calculate_gq_from_pls(hom_ref_calc->final_phred_scaled_genotype_likelihoods));
    g->set_pl(std::move(hom_ref_calc->final_phred_scaled_genotype_likelihoods));

    pGenotypesContext gc = GenotypesContext::create(_pool);
    gc->add(g);

    AlleleVector vc_alleles{{ref_allele, StaticAllele::get_instance()->_non_ref_allele.get()}, _pool};
    pVariant vc = Variant::create(_pool);
    vc->set_tid(pileup->get_tid());
    vc->set_start(pileup->get_start());
    vc->set_stop(pileup->get_stop());
    vc->set_alleles(vc_alleles);
    vc->set_genotype(gc);

    return vc;
}

pRefVsAnyResult ReferenceConfidenceModel::calc_genotype_likelihoods_of_ref_vs_any(int32_t ploidy, pReadPileup pileup)
{
    int32_t likelihood_count = ploidy + 1;
    double log10ploidy = MathUtils::log10(ploidy);

    pRefVsAnyResult result = new ALLOC_TYPE_IN_POOL(_pool, RefVsAnyResult) RefVsAnyResult{_pool};

    uint16_t read_count = 0;
    for (int32_t i = 0; i < FLAG_STATUS; i++) {
        uint16_t max = pileup->qual_max[i] + 1;
        for (uint16_t j = BASE_QUAL_THRESHOLD + 1; j < max; j++) {
            uint32_t count = pileup->qual_his[i][j];
            if (0 == count) {
                continue;
            }
            apply_pileup_element_ref_vs_non_ref_likelihood_and_count(i == FS_NON_REF, likelihood_count, log10ploidy, result, uint8_t(j),
                                                                     count);
            read_count += count;
        }
    }

    double denominator = read_count * log10ploidy;
    for (int32_t i = 0; i < likelihood_count; ++i) {
        result->genotype_likelihoods[i] -= denominator;
    }

    return result;
}

void ReferenceConfidenceModel::do_indel_ref_conf_calc(int32_t ploidy, pBases ref, pReadPileup pileup, int64_t ref_offset,
                                                      pRefVsAnyResult result)
{
    DoubleVector lk = get_genotype_likelihoods_capped_by_hom_ref_likelihood(result);
    pGenotypeLikelihoods snp_gls = GenotypeLikelihoods::create(std::move(lk), _pool);
    int32_t n_informative = calc_nreads_with_no_plausible_indels_reads(pileup, ref_offset, ref);
    pGenotypeLikelihoods indel_gls = get_indel_pls(ploidy, n_informative);

    // now that we have the snp and indel gls, we take the one with the least confidence, as this is the most conservative estimate of our
    // certainty that we are hom-ref. for example, if the snp pls are 0,10,100 and the indel pls are 0,100,1000 we are very certain that
    // there's no indel here, but the snp confidence imply that we are far less confident that the ref base is actually the only thing here.
    // so we take 0,10,100 as our gls for the site.
    pGenotypeLikelihoods least_confidence_gls = get_gl_with_worst_gq(indel_gls, snp_gls);

    result->final_phred_scaled_genotype_likelihoods.reserve(least_confidence_gls->_pls.size());
    std::copy(least_confidence_gls->_pls.begin(), least_confidence_gls->_pls.end(),
              std::back_inserter(result->final_phred_scaled_genotype_likelihoods));
}

void ReferenceConfidenceModel::apply_pileup_element_ref_vs_non_ref_likelihood_and_count(bool is_alt, int32_t likelihood_count,
                                                                                        double log10_ploidy, pRefVsAnyResult result,
                                                                                        uint8_t qual, uint32_t count)
{
    double reference_likelihood, non_ref_likelihood;
    if (is_alt) {
        non_ref_likelihood = QualityUtils::qual_to_prob_log10(qual);
        reference_likelihood = QualityUtils::qual_to_error_prob_log10(qual) + math_utils::s_log10_one_third;
        result->non_ref_depth += int32_t(count);
    }
    else {
        reference_likelihood = QualityUtils::qual_to_prob_log10(qual);
        non_ref_likelihood = QualityUtils::qual_to_error_prob_log10(qual) + math_utils::s_log10_one_third;
        result->ref_depth += int32_t(count);
    }

    result->genotype_likelihoods[0] += double(count * (reference_likelihood + log10_ploidy));
    result->genotype_likelihoods[likelihood_count - 1] += double(count * (non_ref_likelihood + log10_ploidy));

    for (int32_t i = 1, j = likelihood_count - 2; i < likelihood_count - 1; ++i, --j) {
        result->genotype_likelihoods[i] += count * (MathUtils::approximate_log10sum_log10(reference_likelihood + MathUtils::log10(j),
                                                                                          non_ref_likelihood + MathUtils::log10(i)));
    }
}

int32_t ReferenceConfidenceModel::calc_nreads_with_no_plausible_indels_reads(pReadPileup pileup, int64_t ref_offset, pBases ref)
{
    int32_t n_informative = 0;
    int64_t offset;
    const std::pmr::list<pPileupElement>& pes = pileup->get_pileup_element();
    for (pPileupElement pe : pes) {
        if (pe->is_before_deletion_start() || pe->is_before_insertion() || pe->is_deletion()) {
            continue;
        }
        offset = get_cigar_modified_offset(pe);
        if (read_has_no_plausible_ideals_of_size(pe->get_read(), offset, ref, ref_offset)) {
            n_informative++;
            if (n_informative > MAX_N_INDEL_INFORMATIVE_READS) {
                return MAX_N_INDEL_INFORMATIVE_READS;
            }
        }
    }
    return n_informative;
}

int64_t ReferenceConfidenceModel::get_cigar_modified_offset(pPileupElement pe)
{
    uint32_t op = bam_cigar_op(pe->_current_cigar_element);
    auto offset = int64_t((consumes_ref_bases(op) || cigar_op_is_soft_clip(op)) ? pe->_offset_in_current_cigar : 0);
    uint32_t element;
    for (uint32_t i = 0; i < pe->_current_cigar_offset; ++i) {
        element = pe->get_read()->cigar_i(i);
        op = bam_cigar_op(element);
        if (consumes_ref_bases(op) || cigar_op_is_soft_clip(op)) {
            offset += int64_t(bam_cigar_oplen(element));
        }
    }
    return offset;
}

pGenotypeLikelihoods ReferenceConfidenceModel::get_indel_pls(int32_t ploidy, int32_t n)
{
    return s_cache->get_indel_pls(ploidy, n > MAX_N_INDEL_INFORMATIVE_READS ? MAX_N_INDEL_INFORMATIVE_READS : n);
}

pGenotypeLikelihoods ReferenceConfidenceModel::get_snp_pls(DoubleVector&& d, pMemoryPool pool)
{
    return GenotypeLikelihoods::create(std::forward<DoubleVector>(d), pool);
}

DoubleVector ReferenceConfidenceModel::get_genotype_likelihoods_capped_by_hom_ref_likelihood(pRefVsAnyResult result)
{
    DoubleVector ret{_pool};
    ret.reserve(TWO_PLOIDY_LIKELIHOOD_CAPACITY);
    ret.push_back(result->genotype_likelihoods[0]);
    for (int32_t i = 1; i < TWO_PLOIDY_LIKELIHOOD_CAPACITY; ++i) {
        ret.push_back(std::min(ret.front(), result->genotype_likelihoods[i]));
    }
    return ret;
}

bool ReferenceConfidenceModel::read_has_no_plausible_ideals_of_size(pReadRecord read, int64_t read_start, pBases ref_bases,
                                                                    int64_t ref_start)
{
    pInformativeSet informative = read->get_informative_set();
    if (nullptr == informative) {
        CHECK_CONDITION_EXIT(read_start < 0, "read_start must >= 0");
        CHECK_CONDITION_EXIT(ref_start < 0, "ref_start must >= 0");
        unsigned long read_len = (unsigned long)read->seq_length();
        informative = new ALLOC_TYPE_IN_POOL(_pool, InformativeSet) InformativeSet{read_len, 0};
        if (int64_t(read->seq_length()) - read_start >= MAX_INDEL_SIZE && int64_t(ref_bases->num) - ref_start >= MAX_INDEL_SIZE) {
            int64_t secondary_read_break_position = read->seq_length() - MAX_INDEL_SIZE;

            // We are safe to use the faster no-copy versions of getBases and getBaseQualities here,
            // since we're not modifying the returned arrays in any way. This makes a small difference
            // in the HaplotypeCaller profile, since this method is a major hotspot.
            std::pair<pBases, pBases> bq = AlignmentUtils::get_bases_and_base_qualities_aligned_one_to_one(read, _pool);
            pBases read_bases = bq.first;
            pBases read_qualities = bq.second;

            // Need to check for closeness to the end of the read again as the array size may be different than read.Len() due to deletions
            // in the cigar
            if (int64_t(read_bases->num) - read_start > MAX_INDEL_SIZE) {
                // Compute where the end of marking would have been given the above two break conditions so we can stop marking there for
                // our cached results
                bool reference_was_shorter;
                int64_t last_read_base_to_mark_as_indel_relevant;
                if (int64_t(read_bases->num) < int64_t(ref_bases->num) - ref_start + read_start + 1) {
                    // if the read ends first, then we don't mark the last max_indel_size bases from it as relevant
                    last_read_base_to_mark_as_indel_relevant = int64_t(read_bases->num) - MAX_INDEL_SIZE;
                    reference_was_shorter = false;
                }
                else {
                    // if the reference ends first, then we don't mark the last max_indel_size bases from it as relevant
                    last_read_base_to_mark_as_indel_relevant = int64_t(ref_bases->num) - ref_start + read_start - MAX_INDEL_SIZE + 1;
                    reference_was_shorter = true;
                }

                // Compute the absolute baseline sum against which to test
                Int32Vector mis_match_sums = calculate_baseline_mmqualities(read_bases, read_qualities, read_start, ref_bases, ref_start);

                // consider each indel size up to max in term, checking if an indel that deletes either the ref bases (deletion)
                // or read bases (insertion) would fit as well as the origin baseline sum of mismatching quality scores. These scores
                // are computed starting from the last base in the read/reference that would be offset by the indel and compared against
                // the mismatch cost for the same base of the reference. Once the sum of mismatch qualities counting from the back for
                // one indel size exceeds the global indel mismatch cost, the code stops as it will never find a better mismatch value.
                for (int64_t indel_size = 1; indel_size <= MAX_INDEL_SIZE; indel_size++) {
                    // Computing mismatches corresponding to a deletion
                    traverse_end_of_read_for_indel_mismatches(informative, read_start, read_bases, read_qualities,
                                                              last_read_base_to_mark_as_indel_relevant, secondary_read_break_position,
                                                              ref_start, ref_bases, mis_match_sums, indel_size, false);

                    // Computing mismatches corresponding to an insertion
                    traverse_end_of_read_for_indel_mismatches(informative, read_start, read_bases, read_qualities,
                                                              last_read_base_to_mark_as_indel_relevant, secondary_read_break_position,
                                                              ref_start, ref_bases, mis_match_sums, indel_size, true);
                }

                // flip the bases at the front of the read (the ones not within max_indel_size of the end as those are never
                // informative) these must be flipped because thus far we have marked reads for which there were plausible indels with a
                // true value in the bitset. this method returns false for cases where we have discovered plausible indels so we must
                // flip them. this is done in part to preserve a sensible default behavior for bases not considered by this approach.
                if (last_read_base_to_mark_as_indel_relevant <= secondary_read_break_position) {
                    for (int64_t k = 0; k < last_read_base_to_mark_as_indel_relevant; ++k) {
                        informative->flip(k);
                    }
                    // resolve the fact that the old approach would always mark the last base examined as being indel uninformative when
                    // the reference ends first despite it corresponding to a comparison of zero bases against the read
                    if (reference_was_shorter) {
                        informative->set(last_read_base_to_mark_as_indel_relevant - 1, false);
                    }
                }
                else {
                    for (int64_t k = 0; k < secondary_read_break_position + 1; ++k) {
                        informative->flip(k);
                    }
                }
            }
        }
        read->set_informative_set(informative);
    }
    return (size_t)read_start >= informative->size() ? false : informative->test(read_start);
}

Int32Vector ReferenceConfidenceModel::calculate_baseline_mmqualities(pBases read_bases, pBases read_quals, int64_t read_start, pBases ref,
                                                                     int64_t ref_start)
{
    int64_t n = std::min(int64_t(read_bases->num) - read_start, int64_t(ref->num) - ref_start);
    Int32Vector results{size_t(n), _pool};
    int32_t sum = 0;

    uint8_t ref_byte, read_byte;
    for (int64_t i = n - 1; i >= 0; --i) {
        read_byte = read_bases->data[read_start + i];
        ref_byte = ref->data[ref_start + i];
        if (IS_MISMATCH_AND_NOT_AN_ALIGNMENT_GAP(read_byte, ref_byte)) {
            sum += read_quals->data[read_start + i];
        }
        results[i] = sum;
    }
    return results;
}

void ReferenceConfidenceModel::traverse_end_of_read_for_indel_mismatches(pInformativeSet informative, int64_t read_start, pBases read_bases,
                                                                         pBases read_quals,
                                                                         int64_t last_read_base_to_mark_as_indel_relevant,
                                                                         int64_t secondary_read_break_position, int64_t ref_start,
                                                                         pBases ref_bases, const Int32Vector& baseline_mmsums,
                                                                         int64_t indel_size, bool insertion)
{
    int32_t base_quality_sum = 0;
    int32_t global_mismatch_cost_for_read_aligned_to_reference = baseline_mmsums.at(0);

    // compute how many bases forward we should compare taking into account reference/read overhang
    int64_t insertion_length = !insertion ? 0 : indel_size;
    int64_t deletion_length = insertion ? 0 : indel_size;

    // based on the offsets and the indel_size we are considering, how many bases until we fall off the end of the read/reference arrays?
    int64_t first = int64_t(read_bases->num) - read_start - insertion_length;
    int64_t second = int64_t(ref_bases->num) - ref_start - deletion_length;
    int64_t number_of_bases_to_directly_compare = std::min(first, second);

    for (int64_t read_offset = number_of_bases_to_directly_compare + insertion_length - 1,
                 ref_offset = number_of_bases_to_directly_compare + deletion_length - 1;
         read_offset >= 0 && ref_offset >= 0; read_offset--, ref_offset--) {
        // calculate the real base offset for the read:
        uint8_t read_base = read_bases->data[read_start + read_offset];
        uint8_t ref_base = ref_bases->data[ref_start + ref_offset];
        if (IS_MISMATCH_AND_NOT_AN_ALIGNMENT_GAP(read_base, ref_base)) {
            base_quality_sum += read_quals->data[read_start + read_offset];
            if (base_quality_sum > global_mismatch_cost_for_read_aligned_to_reference) {
                // abort early if we are over our global mismatch cost
                break;
            }
        }

        // the hypothetical "read_offset" that corresponds to the comparison we are currently making
        int64_t site_of_real_comparison_point = std::min(read_offset, ref_offset);

        // if it's a real character and the cost isn't greater than the non-indel cost, label it as uninformative
        if (read_bases->data[read_start + site_of_real_comparison_point] != AlignmentUtils::s_gap_base_character &&
            // use less than here because last_read_base_to_mark_as_indel_relevant is the exclusive site where we flip bases later on.
            read_start + site_of_real_comparison_point < last_read_base_to_mark_as_indel_relevant &&
            // resolving the edge case involving read.get_length() disagreeing with the realigned indel length
            read_start + site_of_real_comparison_point <= secondary_read_break_position &&
            baseline_mmsums[site_of_real_comparison_point] >= base_quality_sum) {
            // label with true here because we flip these results later
            informative->set(read_start + site_of_real_comparison_point, true);
        }
    }
}

pGenotypeLikelihoods ReferenceConfidenceModel::get_gl_with_worst_gq(pGenotypeLikelihoods gl1, pGenotypeLikelihoods gl2)
{
    return get_gq_for_hom_ref(gl1) > get_gq_for_hom_ref(gl2) ? gl1 : gl2;
}

double ReferenceConfidenceModel::get_gq_for_hom_ref(pGenotypeLikelihoods gl)
{
    return GenotypeLikelihoods::get_gq_log10from_likelihoods(IDX_HOM_REF, gl->_log10likelihoods, _pool);
}

}  // namespace rovaca