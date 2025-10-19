#ifndef ROVACA_HC_GENOTYPE_ARGUMENT_H_
#define ROVACA_HC_GENOTYPE_ARGUMENT_H_
#include <cmath>
#include <cstdint>
#include <string>
#include <vector>

#include "genotype/genotype_enum.h"

namespace rovaca
{

typedef struct GenotypeArgument HCArgs, *pHCArgs;

struct GenotypeArgument
{
    /* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
public:  // 当前一直使用默认值的参数
    bool applyBQD{false};
    bool applyFRD{false};
    int32_t sample_ploidy{2};
    uint32_t max_mnp_distance{0};
    int32_t max_alternate_alleles{6};
    double contamination_fraction{0.0};
    bool dont_use_dragstr_priors{false};
    const char* dragstr_params{nullptr};
    int64_t assembly_region_padding{100};
    uint8_t mapping_quality_threshold{20};
    bool map_has_contamination_set{false};
    double snp_heterozygosity{0.0010000000};
    double indel_heterozygosity{0.0001250000};
    int64_t informative_read_overlap_margin{2};
    int64_t snp_padding_for_genotyping{20};
    int64_t indel_padding_for_genotyping{75};
    int64_t str_padding_for_genotyping{75};
    OutputMode output_mode{EMIT_VARIANTS_ONLY};
    bool disable_spanning_event_genotyping{false};
    int32_t minimum_read_length_after_trimming{10};
    bool use_filtered_read_map_for_annotations{false};
    double heterozygosity_standard_deviation{0.0100000000};
    bool use_posterior_probabilities_to_calculate_qual{false};
    GenotypeAssignmentMethod genotype_assignment_method{USE_PLS_TO_ASSIGN};

    /* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
public:  // 需要初始化的参数
    const char* output;
    std::string idx_name{};
    const char* tool_name;
    std::string command_line;
    std::vector<int32_t> gvcf_gq_bands;
    ReferenceConfidenceMode reference_confidence_mode{ReferenceConfidenceMode::NONE};

    /* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
public:                                      // 依赖参数，依赖其他参数
    double standard_confidence_for_calling;  // gvcf 0, vcf 30
    bool annotate_all_sites_with_pls;        // gvcf true, vcf false
    bool do_not_run_physical_phasing;        // gvcf false, vcf true

    /* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
public:
    /*! @return true if there is some sample contamination present, false otherwise */
    bool is_sample_contamination_present() const
    {
        return contamination_fraction_is_set(contamination_fraction) || map_has_contamination_set;
    }

    bool emit_reference_confidence() const { return reference_confidence_mode != ReferenceConfidenceMode::NONE; }

    void init_reference_confidence_mode(ReferenceConfidenceMode mode)
    {
        reference_confidence_mode = mode;
        standard_confidence_for_calling = reference_confidence_mode == ReferenceConfidenceMode::GVCF ? 0 : 30;
        annotate_all_sites_with_pls = reference_confidence_mode == ReferenceConfidenceMode::GVCF;
        do_not_run_physical_phasing = reference_confidence_mode != ReferenceConfidenceMode::GVCF;
    }

    /* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
private:
    static bool contamination_fraction_is_set(double fraction) { return !std::isnan(fraction) && fraction > 0.0; }
};

}  // namespace rovaca

#endif  // ROVACA_HC_GENOTYPE_ARGUMENT_H_
