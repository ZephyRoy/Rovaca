#ifndef HAPLOTYPECALLER_ARGS_COLLECTOR_H
#define HAPLOTYPECALLER_ARGS_COLLECTOR_H

#include <assemble_argument.h>

#include <iostream>
#include <string>
#include <vector>

#include "common/enum.h"
#include "genotype/genotype_argument.h"

struct RegionArgument
{
    std::string ASSEMBLY_REGION_OUT_LONG_NAME = "assembly-region-out";
    std::string FORCE_ACTIVE_REGIONS_LONG_NAME = "force-active";
    long serialVersionUID = 1L;
    std::string MIN_ASSEMBLY_LONG_NAME = "min-assembly-region-size";
    std::string MAX_ASSEMBLY_LONG_NAME = "max-assembly-region-size";
    std::string ASSEMBLY_PADDING_LONG_NAME = "assembly-region-padding";
    std::string MAX_STARTS_LONG_NAME = "max-reads-per-alignment-start";
    std::string THRESHOLD_LONG_NAME = "active-probability-threshold";
    std::string PROPAGATION_LONG_NAME = "max-prob-propagation-distance";
    int DEFAULT_MIN_ASSEMBLY_REGION_SIZE = 50;
    int DEFAULT_MAX_ASSEMBLY_REGION_SIZE = 300;
    int DEFAULT_ASSEMBLY_REGION_PADDING = 100;
    int DEFAULT_MAX_READS_PER_ALIGNMENT = 50;
    double DEFAULT_ACTIVE_PROB_THRESHOLD = 0.002;
    int DEFAULT_MAX_PROB_PROPAGATION_DISTANCE = 50;
    std::string INDEL_PADDING_LONG_NAME = "padding-around-indels";
    std::string SNP_PADDING_LONG_NAME = "padding-around-snps";
    std::string STR_PADDING_LONG_NAME = "padding-around-strs";
};

struct PairhmmArgument
{
    int32_t base_quality_score_threshold{18};  // default 18
    PcrIndelModel pcr_option{PcrIndelModel::CONSERVATIVE};
};

class HaplotypeCallerArgs
{
public:
    RegionArgument region_arguments_;
    AssembleArgument assemble_arguments_;
    PairhmmArgument pairhmm_arguments_;
    rovaca::GenotypeArgument genotype_args_;

    enum FlowMode { NONE, STANDARD, ADVANCED };
    struct DbsnpArgument
    {};

    DbsnpArgument dbsnp;
    std::vector<int32_t> GQBands_{};
    FlowMode flowMode = NONE;
    bool applyBQD;
    bool applyFRD;
    bool floorBlocks;
    bool disableOptimizations;
    bool dragenMode;
    bool disableSpanningEventGenotyping;
    bool transformDRAGENMapQ;
    bool justDetermineActiveRegions;
    bool dontGenotype;
    bool doNotRunPhysicalPhasing;
    bool doNotCorrectOverlappingBaseQualities;
    bool useFilteredReadMapForAnnotations;
    bool stepwiseFiltering;
    int mappingQualityThreshold;
    int maxEffectiveDepthAdjustment;
    int indelSizeToEliminateInRefModel;

    PcrIndelModel pcrIndelModel;
    ReferenceConfidenceMode referenceConfidenceMode;

    std::string keepRG;
    std::string assemblyStateOutput;
    std::string genotyperDebugOutStream;
    std::string sampleNameToUse;
    std::string getDragenNameValuePairs() { return ""; };

public:
    HaplotypeCallerArgs();
    ~HaplotypeCallerArgs();
    void initGQBands(const std::vector<int32_t>& specified_gq);
};
#endif  // HAPLOTYPECALLER_ARGS_COLLECTOR_H