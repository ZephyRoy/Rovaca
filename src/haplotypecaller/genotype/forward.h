#ifndef ROVACA_HC_FORWARD_H_
#define ROVACA_HC_FORWARD_H_
#include <bitset>
#include <boost/dynamic_bitset.hpp>
#include <list>
#include <map>
#include <memory_resource>
#include <set>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <vector>

struct bam1_t;
struct sam_hdr_t;
typedef struct sam_hdr_t bam_hdr_t;
typedef struct lib_sw_avx_st lib_sw_avx, *p_lib_sw_avx;
namespace rovaca
{

typedef std::pmr::memory_resource MemoryPool, *pMemoryPool;
typedef class BamDataPool BamDataPool, *pBamDataPool;

typedef struct Bases Bases, *pBases;
typedef struct Cigar Cigar, *pCigar;
typedef struct GenotypeArgument HCArgs, *pHCArgs;
typedef struct RefFragment RefFragment, *pRefFragment;
typedef struct GenotypeCounts GenotypeCounts, *pGenotypeCounts;

typedef class Allele Allele, *pAllele;
typedef class Variant Variant, *pVariant;
typedef class EventMap EventMap, *pEventMap;
typedef class InfoData InfoData, *pInfoData;
typedef class Genotype Genotype, *pGenotype;
typedef class Haplotype Haplotype, *pHaplotype;
typedef class ReadRecord ReadRecord, *pReadRecord;
typedef class ReadPileup ReadPileup, *pReadPileup;
typedef class HomRefBlock HomRefBlock, *pHomRefBlock;
typedef class ConstantsStr ConstantsStr, *pConstantsStr;
typedef class CigarBuilder CigarBuilder, *pCigarBuilder;
typedef class PileupElement PileupElement, *pPileupElement;
typedef class BlockCombiner BlockCombiner, *pBlockCombiner;
typedef class SimpleInterval SimpleInterval, *pSimpleInterval;
typedef class GenotypesContext GenotypesContext, *pGenotypesContext;
typedef class InterfaceLocatable InterfaceLocatable, *pInterfaceLocatable;
typedef class GenotypeLikelihoods GenotypeLikelihoods, *pGenotypeLikelihoods;
typedef class AFCalculationResult AFCalculationResult, *pAFCalculationResult;
typedef class InterfaceSampleList InterfaceSampleList, *pInterfaceSampleList;
typedef class InterfacePloidyModel InterfacePloidyModel, *pInterfacePloidyModel;
typedef class GenotypeAlleleCounts GenotypeAlleleCounts, *pGenotypeAlleleCounts;
typedef class VariantAnnotatorEngine VariantAnnotatorEngine, *pVariantAnnotatorEngine;
typedef class GenotypePriorCalculator GenotypePriorCalculator, *pGenotypePriorCalculator;
typedef class GermlineGenotyingEngine GermlineGenotyingEngine, *pGermlineGenotyingEngine;
typedef class GenotypeLikelihoodsCache GenotypeLikelihoodsCache, *pGenotypeLikelihoodsCache;
typedef class AlleleFrequencyCalculator AlleleFrequencyCalculator, *pAlleleFrequencyCalculator;
typedef class GenotypeAlleleCountsManger GenotypeAlleleCountsManger, *pGenotypeAlleleCountsManger;
typedef class InterfaceGenotypeAnnotation InterfaceGenotypeAnnotation, *pInterfaceGenotypeAnnotation;
typedef class InterfaceInfoFieldAnnotation InterfaceInfoFieldAnnotation, *pInterfaceInfoFieldAnnotation;
typedef class GenotypeLikelihoodCalculator GenotypeLikelihoodCalculator, *pGenotypeLikelihoodCalculator;
typedef class IndependentSampleGenotypesModel IndependentSampleGenotypesModel, *pIndependentSampleGenotypesModel;

template <typename EVIDENCE, typename A>
class AlleleLikelihoods;
typedef AlleleLikelihoods<pReadRecord, pAllele> RALikelihoods, *pRALikelihoods;
typedef AlleleLikelihoods<pReadRecord, pHaplotype> RHLikelihoods, *pRHLikelihoods;

typedef boost::dynamic_bitset<unsigned long> InformativeSet, *pInformativeSet;
typedef std::pmr::set<int64_t> Int64Set;
typedef std::pmr::set<pAllele> AlleleSet;
typedef std::pmr::vector<bool> BoolVector;
typedef std::pmr::vector<long> LongVector;
typedef std::pmr::list<pReadRecord> ReadList;
typedef std::pmr::vector<int32_t> Int32Vector;
typedef std::pmr::vector<double> DoubleVector;
typedef std::pmr::vector<pAllele> AlleleVector;
typedef std::pmr::set<pHaplotype> HaplotypeSet;
typedef std::pmr::vector<uint32_t> Uint32Vector;
typedef std::pmr::vector<pReadRecord> ReadVector;
typedef std::pmr::vector<ReadVector> ReadVector2D;
typedef std::pmr::vector<pVariant> VariantVector;
typedef std::pmr::list<pHaplotype> HaplotypeList;
typedef std::pmr::set<std::string> StringSet;
typedef std::pmr::vector<Int32Vector> Int32Vector2D;
typedef std::pair<int64_t, uint32_t> Int64Uint32Pair;
typedef std::pmr::vector<pHaplotype> HaplotypeVector;
typedef std::pmr::vector<DoubleVector> DoubleVector2D;
typedef std::pmr::vector<pReadPileup> ReadPileupVector;
typedef std::pmr::vector<DoubleVector2D> DoubleVector3D;
typedef std::pmr::unordered_set<pReadRecord> ReadHashSet;
typedef std::pmr::map<pAllele, double> AlleleToDoubleMap;
typedef std::pmr::map<int64_t, pVariant> Int64ToVariantMap;
typedef std::pmr::map<int32_t, ReadVector> Int32ToReadVectorMap;
typedef std::pair<VariantVector, HaplotypeList> CalledHaplotypes;
typedef std::pair<pSimpleInterval, pSimpleInterval> IntervalPair;
typedef std::pmr::map<pVariant, HaplotypeSet> VariantToHaplotypeSetMap;
typedef std::pmr::map<pAllele, HaplotypeVector> AlleleToHaplotypeVectorMap;

typedef std::pair<bool, double> OptionalDouble;
typedef std::pair<bool, uint8_t> OptionalUint8;

}  // namespace rovaca

#endif  // ROVACA_HC_FORWARD_H_
