#include "haplotypecaller_args.h"

#include <algorithm>

HaplotypeCallerArgs::HaplotypeCallerArgs()
    : region_arguments_()
    , assemble_arguments_(ASSEMBLE_DEFAULT_ARGUMENTS)
    , genotype_args_()
    , applyBQD(false)
    , applyFRD(false)
    , floorBlocks(false)
    , disableOptimizations(false)
    , dragenMode(false)
    , disableSpanningEventGenotyping(false)
    , transformDRAGENMapQ(false)
    , justDetermineActiveRegions(false)
    , dontGenotype(false)
    , doNotRunPhysicalPhasing(false)
    , doNotCorrectOverlappingBaseQualities(false)
    , useFilteredReadMapForAnnotations(false)
    , stepwiseFiltering(false)
    , mappingQualityThreshold(20)
    , maxEffectiveDepthAdjustment(0)
    , indelSizeToEliminateInRefModel(10)
    , keepRG()
    , assemblyStateOutput()
    , genotyperDebugOutStream()
    , sampleNameToUse()
{}

HaplotypeCallerArgs::~HaplotypeCallerArgs() {}

void HaplotypeCallerArgs::initGQBands(const std::vector<int32_t>& specified_gq)
{
    if (specified_gq.empty()) {
        for (int32_t i = 1; i <= 60; ++i) {
            GQBands_.push_back(i);
        }
        GQBands_.push_back(70);
        GQBands_.push_back(80);
        GQBands_.push_back(90);
        GQBands_.push_back(99);
    }
    else {
        std::copy(specified_gq.begin(), specified_gq.end(), std::back_inserter(GQBands_));
    }
}
