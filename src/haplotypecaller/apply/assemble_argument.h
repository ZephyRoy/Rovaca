#ifndef ASSEMBLE_ARGUMENT_H
#define ASSEMBLE_ARGUMENT_H

#include <cmath>
#include <string>
#include <vector>
#include <cstdint>

static const double DEFAULT_PRUNING_LOG_ODDS_THRESHOLD = 0.0;          // log10ToLog(1.0)
static const double DEFAULT_PRUNING_SEEDING_LOG_ODDS_THRESHOLD = 0.0;  // log10ToLog(4.0)

static const float PREFILTER_QUAL_THRESHOLD = 30;
static const float PREFILTER_SOR_THRESHOLD = 3;

struct ReadThreadingArgument
{
    std::vector<uint32_t> kmer;
    bool dontIncreaseKmerSizesForCycles;        // Disable iterating over kmer sizes when graph cycles are detected..
    bool allowNonUniqueKmersInRef;              // Allow graphs that have non-unique kmers in the reference..
    bool recoverAllDanglingBranches;            // Recover all dangling branches..
    bool useLinkedDeBruijnGraph;                // Never construct a Linked De Bruijn graph to recover better haplotypes..
    bool disableArtificialHaplotypeRecovery;    // disable recovery of haplotypes based on graph edges that are not in junction trees..
    bool debugGraphTransformations;             // Write DOT formatted graph files out of the assembler for only this graph size.
    bool captureAssemblyFailureBAM;             // Write a BAM called assemblyFailure.bam for all failure reads..
    bool errorCorrectReads;                     // Use an exploratory algorithm to error correct the kmers used during assembly.
    uint32_t numPruningSamples;                 // Number of samples that must pass the minPruning threshold.
    uint32_t minDanglingBranchLength;           // Minimum length of a dangling branch to attempt recovery.
    uint32_t maxNumHaplotypesInPopulation;      // Maximum number of haplotypes to consider for your population.
    uint32_t minPruneFactor;                    // Minimum support to not prune paths in the graph.
    uint32_t maxUnprunedVariants;               // Maximum number of variants in graph the adaptive pruner will allow.
    uint32_t kmerLengthForReadErrorCorrection;  // Use an exploratory algorithm to error correct the kmers used during assembly.
    uint32_t minObservationsForKmerToBeSolid;   // A k-mer must be seen at least these times for it considered to be solid.
    int minMatchingBasesToDanglingEndRecovery;  //.
    double initialErrorRateForPruning;          // Initial base error rate estimate for adaptive pruning.
    double pruningLogOddsThreshold;             // Ln likelihood ratio threshold for adaptive pruning algorithm.
    double pruningSeedingLogOddsThreshold;      // threshold for seeding subgraph of good variation in adaptive pruning algorithm.
    double pileupErrorCorrectionLogOdds;        // Log odds threshold for pileup error correction.  Off by default.
    std::string debugAssemblyVariantsOut;       // Write a VCF from the event map for each active region.
    std::string graphOutput;                    // Write debug assembly graph information to this file.
    std::string haplotypeHistogramOutput;       // Write debug assembly graph information to this file.
    ReadThreadingArgument()
        : kmer({10, 25})
        , dontIncreaseKmerSizesForCycles(false)
        , allowNonUniqueKmersInRef(false)
        , recoverAllDanglingBranches(false)
        , useLinkedDeBruijnGraph(false)
        , disableArtificialHaplotypeRecovery(false)
        , debugGraphTransformations(false)
        , captureAssemblyFailureBAM(false)
        , errorCorrectReads(false)
        , numPruningSamples(1)
        , minDanglingBranchLength(4)
        , maxNumHaplotypesInPopulation(128)
        , minPruneFactor(2)
        , maxUnprunedVariants(100)
        , kmerLengthForReadErrorCorrection(25)
        , minObservationsForKmerToBeSolid(20)
        , minMatchingBasesToDanglingEndRecovery(-1)
        , initialErrorRateForPruning(0.001)
        , pruningLogOddsThreshold(DEFAULT_PRUNING_LOG_ODDS_THRESHOLD)
        , pruningSeedingLogOddsThreshold(DEFAULT_PRUNING_SEEDING_LOG_ODDS_THRESHOLD)
        , pileupErrorCorrectionLogOdds(-INFINITY)
        , debugAssemblyVariantsOut("")
        , graphOutput("")
        , haplotypeHistogramOutput("")
    {}
};
struct SmithWatermanArgument
{
    int smithWatermanDanglingEndMatchValue;
    int smithWatermanDanglingEndMismatchPenalty;
    int smithWatermanDanglingEndGapOpenPenalty;
    int smithWatermanDanglingEndGapExtendPenalty;
    int smithWatermanHaplotypeToReferenceMatchValue;
    int smithWatermanHaplotypeToReferenceMismatchPenalty;
    int smithWatermanHaplotypeToReferenceGapOpenPenalty;
    int smithWatermanHaplotypeToReferenceGapExtendPenalty;
    int smithWatermanReadToHaplotypeMatchValue;
    int smithWatermanReadToHaplotypeMismatchPenalty;
    int smithWatermanReadToHaplotypeGapOpenPenalty;
    int smithWatermanReadToHaplotypeGapExtendPenalty;
    SmithWatermanArgument()
        : smithWatermanDanglingEndMatchValue()
        , smithWatermanDanglingEndMismatchPenalty()
        , smithWatermanDanglingEndGapOpenPenalty()
        , smithWatermanDanglingEndGapExtendPenalty()
        , smithWatermanHaplotypeToReferenceMatchValue()
        , smithWatermanHaplotypeToReferenceMismatchPenalty()
        , smithWatermanHaplotypeToReferenceGapOpenPenalty()
        , smithWatermanHaplotypeToReferenceGapExtendPenalty()
        , smithWatermanReadToHaplotypeMatchValue()
        , smithWatermanReadToHaplotypeMismatchPenalty()
        , smithWatermanReadToHaplotypeGapOpenPenalty()
        , smithWatermanReadToHaplotypeGapExtendPenalty()
    {}
};

struct AssembleArgument
{
    ReadThreadingArgument read_threading_argument;
    SmithWatermanArgument sw_argument;
    bool debugAssembly;
    bool dontUseSoftClippedBases;
    bool overrideSoftclipFragmentCheck;
    bool forceCallFiltered;
    bool softClipLowQualityEnds;
    bool flowAssemblyCollapsePartialMode;
    bool filterAlleles;
    bool filterLoneAlleles;
    bool writeFilteringGraphs;
    uint8_t minBaseQualityScore;
    uint8_t refModelDelQual;
    int maxMnpDistance;
    int informativeReadOverlapMargin;
    int flowAssemblyCollapseHKerSize;
    float prefilterQualThreshold;
    float prefilterSorThreshold;

    AssembleArgument()
        : read_threading_argument()
        , sw_argument()
        , debugAssembly(false)
        , dontUseSoftClippedBases(false)
        , overrideSoftclipFragmentCheck(false)
        , forceCallFiltered(false)
        , softClipLowQualityEnds(false)
        , flowAssemblyCollapsePartialMode(false)
        , filterAlleles(false)
        , filterLoneAlleles(false)
        , writeFilteringGraphs(false)
        , minBaseQualityScore(10)
        , refModelDelQual(30)
        , maxMnpDistance()  // from super
        , informativeReadOverlapMargin(2)
        , flowAssemblyCollapseHKerSize(0)
        , prefilterQualThreshold(30)
        , prefilterSorThreshold(3.0)
    {}
};

static AssembleArgument ASSEMBLE_DEFAULT_ARGUMENTS{};
#endif  // ASSEMBLE_ARGUMENT_H