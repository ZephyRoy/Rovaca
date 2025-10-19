#ifndef ROVACA_HC_REFERENCE_CONFIDENCE_MODEL_H_
#define ROVACA_HC_REFERENCE_CONFIDENCE_MODEL_H_
#include "forward.h"
#include "genotype_macors.h"

namespace rovaca
{

#define TWO_PLOIDY_LIKELIHOOD_CAPACITY (3)

typedef struct RefVsAnyResult
{
    int32_t ref_depth;
    int32_t non_ref_depth = 0;
    DoubleVector genotype_likelihoods;
    Int32Vector final_phred_scaled_genotype_likelihoods;

    explicit RefVsAnyResult(pMemoryPool pool)
        : ref_depth(0)
        , non_ref_depth(0)
        , genotype_likelihoods(TWO_PLOIDY_LIKELIHOOD_CAPACITY, pool)
        , final_phred_scaled_genotype_likelihoods(pool)
    {}
} RefVsAnyResult, *pRefVsAnyResult;

class ReferenceConfidenceModel
{
    /* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
private:
    static pGenotypeLikelihoodsCache s_cache;

    /* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
private:
    pMemoryPool _pool;

    /* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
public:
    DISALLOW_COPY_AND_ASSIGN(ReferenceConfidenceModel);
    explicit ReferenceConfidenceModel(pMemoryPool pool)
        : _pool(pool)
    {}

    VariantVector calculate_ref_confidence(const VariantVector& calls, pHaplotype ref_h, pSimpleInterval flank,
                                           pSimpleInterval flank_padded, int32_t ploidy, const ReadHashSet& reads);

    pVariant make_reference_confidence_variant_context(int32_t ploidy, pBases ref, int32_t sample_idx, int64_t global_ref_offset,
                                                       pReadPileup pileup, int64_t offset);

    pRefVsAnyResult calc_genotype_likelihoods_of_ref_vs_any(int32_t ploidy, pReadPileup pileup);

    void do_indel_ref_conf_calc(int32_t ploidy, pBases ref, pReadPileup pileup, int64_t ref_offset, pRefVsAnyResult result);

    /* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
private:
    static void apply_pileup_element_ref_vs_non_ref_likelihood_and_count(bool is_alt, int32_t likelihood_count, double log10_ploidy,
                                                                         pRefVsAnyResult result, uint8_t qual, uint32_t count);

    int32_t calc_nreads_with_no_plausible_indels_reads(pReadPileup pileup, int64_t ref_offset, pBases ref);

    static int64_t get_cigar_modified_offset(pPileupElement pe);
    static pGenotypeLikelihoods get_indel_pls(int32_t ploidy, int32_t n);
    static pGenotypeLikelihoods get_snp_pls(DoubleVector&& d, pMemoryPool pool);
    DoubleVector get_genotype_likelihoods_capped_by_hom_ref_likelihood(pRefVsAnyResult result);

    bool read_has_no_plausible_ideals_of_size(pReadRecord read, int64_t read_start, pBases ref, int64_t ref_start);

    Int32Vector calculate_baseline_mmqualities(pBases read_bases, pBases read_quals, int64_t read_start, pBases ref, int64_t ref_start);

    static void traverse_end_of_read_for_indel_mismatches(pInformativeSet informative, int64_t read_start, pBases read_bases,
                                                          pBases read_quals, int64_t last_read_base_to_mark_as_indel_relevant,
                                                          int64_t secondary_read_break_position, int64_t ref_start, pBases ref_bases,
                                                          const Int32Vector& baseline_mmsums, int64_t indel_size, bool insertion);

    /**
     * get the GenotypeLikelihoods with the least strong corresponding GQ value
     * @param gl1 first to consider (cannot be nullptr)
     * @param gl2 second to consider (cannot be nullptr)
     * @return gl1 or gl2, whichever has the worst gq
     */
    pGenotypeLikelihoods get_gl_with_worst_gq(pGenotypeLikelihoods gl1, pGenotypeLikelihoods gl2);
    double get_gq_for_hom_ref(pGenotypeLikelihoods gl);
};

}  // namespace rovaca

#endif  // ROVACA_HC_REFERENCE_CONFIDENCE_MODEL_H_
