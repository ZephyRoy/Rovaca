#ifndef ROVACA_HC_GERMLINE_GENOTYING_ENGINE_H_
#define ROVACA_HC_GERMLINE_GENOTYING_ENGINE_H_
#include <algorithm>
#include <functional>
#include <list>
#include <queue>

#include "common/base/include/smithwaterman_common.h"
#include "forward.h"
#include "genotype_macors.h"
#include "htslib/vcf.h"
#include "interface/interface_locatable.hpp"

namespace rovaca
{

typedef struct OutputAlleleSubset OutputAlleleSubset, *pOutputAlleleSubset;

class GermlineGenotyingEngine
{
    /* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
private:
    // 以下属性由外部管理, 传入Engine
    pHCArgs _args{nullptr};
    bam_hdr_t* _header{nullptr};
    pMemoryPool _pool{nullptr};
    pBamDataPool _bam_pool{nullptr};
    pInterfaceSampleList _sample_list{nullptr};
    pInterfacePloidyModel _ploidy_model{nullptr};

    // 以下属性需要由 Engine 管理
    p_lib_sw_avx sw_;
    pBlockCombiner _block_combiner{nullptr};
    pVariantAnnotatorEngine _annotator_engine{nullptr};
    pGenotypePriorCalculator _prior_calculator{nullptr};
    pIndependentSampleGenotypesModel _genotyping_model{nullptr};
    pAlleleFrequencyCalculator _allele_frequency_calculator{nullptr};

    // 跨 region deletion 存储
    std::set<pInterfaceLocatable> _upstream_deletions_loc;
    std::vector<pInterfaceLocatable> _deletions_locs;

    // dbsnp 相关
    int32_t db_offset_{};
    std::vector<bcf1_t*>* db_data_{nullptr};

    /* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
public:
    DISALLOW_COPY_AND_ASSIGN(GermlineGenotyingEngine);
    GermlineGenotyingEngine();
    ~GermlineGenotyingEngine();

    void init_engine_per_loop(pHCArgs args, pMemoryPool pool, pBamDataPool bpool, bam_hdr_t* header, bcf_hdr_t* bcf_hdr,
                              pInterfaceSampleList samples, pInterfacePloidyModel pm);

    CalledHaplotypes assign_genotype_likelihoods(pRHLikelihoods rh_likelihoods, pRefFragment ref, pSimpleInterval ref_loc,
                                                 pSimpleInterval active_region, const Int32ToReadVectorMap& per_sample_filtered_read);

    /*!
     * @brief gvcf 模式计算非 active 位点信息
     * @param ref_h 用于 assign_genotype_likelihoods 计算的单倍体中的 ref
     * @param ref 用于 assign_genotype_likelihoods 计算的 ref
     * @param ref_loc ref 对应的 loc
     * @param original 最原始的
     * @param original_padded original 左右扩100
     * @param variant 包含所有变异点的最小区间
     * @param variant_padded variant 左右扩 20
     * @param calls assign_genotype_likelihoods 的返回结果
     * @param reads
     * @return
     */
    VariantVector call_non_active_site(pHaplotype ref_h, pRefFragment ref, pSimpleInterval ref_loc, pSimpleInterval original,
                                       pSimpleInterval original_padded, pSimpleInterval variant, pSimpleInterval variant_padded,
                                       const CalledHaplotypes& calls, const ReadHashSet& original_reads, const ReadHashSet& genotype_reads,
                                       std::pmr::list<bam1_t*>& extra_memory_reads);

    /*!
     * @brief
     * @param ref
     * @param ref_loc original_padded基础上左右扩 500
     * @param active
     * @param padded
     * @param ploidy
     * @param reads
     * @return
     */
    VariantVector reference_model_for_no_variation(pRefFragment ref, pSimpleInterval ref_loc, pSimpleInterval active,
                                                   pSimpleInterval padded, int32_t ploidy, ReadHashSet& reads);

    /*! @brief 将start不等于loc的Variant替换为start等于loc且allele为{ref_allele, span_del} */
    VariantVector replace_span_dels(const VariantVector& events_at_this_loc, pAllele ref_allele, int64_t loc);

    /*!
     * @brief main entry function to calculate genotypes of a given vc with corresponding gls that is shared across genotypers (namely
     * ggvcfs and hc). completes a variant context with genotype calls and associated annotations given the genotype likelihoods and the
     * model that need to be applied.  hom-ref likelihoods can be approximated from gqs, but if no genotype has likelihoods then that
     * variant is either all-ref or contains variants with no likelihoods, and in both cases we want to exit.
     * @note 对于HC, 此函数仅用于计算genotype中的Allele,GQ 以及info列中的MLEAC、MLEAF，很多冗余代码
     * @return
     */
    pVariant calculate_genotypes(pVariant merged_vc, const AlleleVector& given_alleles);

    static bool contains_calls(const VariantVector& calls);

    void clear_upstream_deletions_loc();
    void record_deletions(pVariant vc, const AlleleVector& alleles);
    bool is_vc_covered_by_deletion(pVariant vc);

    bool emit_reference_confidence() const;

    p_lib_sw_avx sw() { return sw_; }
    pBlockCombiner block_combiner() { return _block_combiner; }

    /// @brief 设置dbsnp相关数据，每个pack设置一次即可
    /// @param db_offset
    /// @param db_data
    void set_dbsnp(int32_t db_offset, std::vector<bcf1_t*>* db_data)
    {
        db_offset_ = db_offset;
        db_data_ = db_data;
    }

    /* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
private:
    static pGenotypePriorCalculator resolve_genotype_prior_calculator(pHCArgs args);
    static bool cannot_be_genotyped(pVariant vc);

    static bool is_spanning_deletion(pAllele allele);
    static bool no_alleles_or_first_allele_is_not_non_ref(const AlleleVector& alleles);
    bool force_keep_allele(pAllele allele);
    bool emit_all_active_sites() const;
    bool passes_call_threshold(double conf) const;
    bool passes_emit_threshold(double conf, bool best_guess_is_ref) const;

    pVariant make_annotated_call(pRefFragment ref, pVariant call, pRALikelihoods ra_likelihoods, size_t num_before_remove);

    /*!
     * @brief For a particular event described in inputVC, form PL vector for each sample by looking into allele read map and filling
     * likelihood matrix for each allele
     */
    pGenotypesContext calculate_gls_for_this_event(pRALikelihoods likelihoods, pVariant merged_vc, pRefFragment ref,
                                                   int64_t offset_for_ref_into_event);

    /**
     * provided the exact mode computations it returns the appropriate subset of alleles that progress to genotyping.
     * @param af_calculation_result the allele fraction calculation result.
     * @param vc the variant context
     * @param forced_alleles alleles from the vc input that are consistent with forced alleles in the assembly region {@link
     * assembly_based_caller_utils#get_alleles_consistent_with_given_alleles}
     * @return information about the alternative allele subsetting {@code null}.
     */
    pOutputAlleleSubset calculate_output_allele_subset(pAFCalculationResult af_calculation_result, pVariant merged_vc,
                                                       const AlleleSet& forced_alleles);

    void compose_call_attributes(pVariant vc, Int32Vector& mleac, pAFCalculationResult af_result, const AlleleVector& all_alleles_to_use,
                                 pGenotypesContext genotypes, pInfoData info);

    // Builds the read-likelihoods collection to use for annotation considering user arguments and the collection used for genotyping
    pRALikelihoods prepare_read_allele_likelihoods_for_annotation(pRHLikelihoods rh_likelihoods, const Int32ToReadVectorMap& filtered_read,
                                                                  bool erc, const AlleleToHaplotypeVectorMap& allele_mapper,
                                                                  pRALikelihoods ra_likelihoods, pVariant call, pSimpleInterval loc);

    Int32ToReadVectorMap overlapping_filtered_reads(const Int32ToReadVectorMap& filtered_read, pSimpleInterval loc);

    static pVariant remove_alt_alleles_if_too_many_genotypes(int32_t ploidy, const AlleleToHaplotypeVectorMap& allele_mapper,
                                                             pVariant merged_vc);

    void init_engine();
};

}  // namespace rovaca

#endif  // ROVACA_HC_GERMLINE_GENOTYING_ENGINE_H_
