#ifndef ROVACA_HC_UNIT_TEST_HEADER_H_
#define ROVACA_HC_UNIT_TEST_HEADER_H_
#include <algorithm>
#include <string>
#include <unordered_map>

#include "forward.h"
#include "genotype_struct.h"

using namespace rovaca;

struct FindNumberOfRepetitionsData
{
    std::string repeat_unit_;
    std::string test_string_;
    bool leading_repeats_;
    int32_t expected_;

    FindNumberOfRepetitionsData(const char* repeat_unit, const char* test_string, bool leading_repeats, int32_t expected)
        : repeat_unit_(repeat_unit)
        , test_string_(test_string)
        , leading_repeats_(leading_repeats)
        , expected_(expected)
    {}
};

struct MergeAllelesTest
{
    std::pmr::vector<std::pmr::vector<pAllele>> inputs;
    std::pmr::vector<pAllele> expected;

    MergeAllelesTest(const std::vector<std::vector<pAllele>>& i, const std::vector<pAllele>& e, pMemoryPool pool)
        : inputs(pool)
        , expected(pool)
    {
        inputs.resize(i.size());
        for (size_t j = 0, len = i.size(); j < len; ++j) {
            const std::vector<pAllele>& ii = i.at(j);
            std::for_each(ii.begin(), ii.end(), [&](pAllele a) { inputs[j].push_back(a); });
        }

        expected.reserve(e.size());
        for (auto j : e) {
            expected.push_back(j);
        }
    }
};

struct MergeAllelesData
{
    int32_t ploidy_;
    GenotypeAssignmentMethod mode_;
    DoubleVector likelihoods_;
    AlleleVector originalGT_;
    AlleleVector allelesToUse_;
    AlleleVector expectedAlleles_;

    MergeAllelesData(pMemoryPool pool, int32_t ploidy, GenotypeAssignmentMethod mode, const std::vector<double>& likelihoods,
                     const std::vector<pAllele>& originalGT, const std::vector<pAllele>& allelesToUse,
                     const std::vector<pAllele>& expectedAlleles)
        : ploidy_(ploidy)
        , mode_(mode)
        , likelihoods_(pool)
        , originalGT_(pool)
        , allelesToUse_(pool)
        , expectedAlleles_(pool)
    {
        std::copy(likelihoods.begin(), likelihoods.end(), std::back_inserter(likelihoods_));
        std::copy(originalGT.begin(), originalGT.end(), std::back_inserter(originalGT_));
        std::copy(allelesToUse.begin(), allelesToUse.end(), std::back_inserter(allelesToUse_));
        std::copy(expectedAlleles.begin(), expectedAlleles.end(), std::back_inserter(expectedAlleles_));
    }
};

struct ClipAllelesData
{
    std::vector<std::string> alleleStrings;
    std::vector<std::string> expected;
    int32_t numLeftClipped;
};

enum LineType {
    k_num = 0,
    k_original,
    k_original_padded,
    k_variant,
    k_variant_padded,
    k_ref_loc,
    k_ref_bases,
    k_original_reads,
    k_trimed_reads,
    k_filtered_reads1,
    k_filtered_reads2,
    k_pairhmm_reads,
    k_genotype_reads,
    k_original_haplotype,
    k_trimed_haplotype,
    k_likelihoods
};

extern const std::unordered_map<std::string, LineType> s_str2type{{"num", k_num},
                                                                  {"original", k_original},
                                                                  {"original_padded", k_original_padded},
                                                                  {"variant", k_variant},
                                                                  {"variant_padded", k_variant_padded},
                                                                  {"ref_loc", k_ref_loc},
                                                                  {"ref_bases", k_ref_bases},
                                                                  {"original_reads", k_original_reads},
                                                                  {"trimed_reads", k_trimed_reads},
                                                                  {"filtered_reads1", k_filtered_reads1},
                                                                  {"filtered_reads2", k_filtered_reads2},
                                                                  {"pairhmm_reads", k_pairhmm_reads},
                                                                  {"genotype_reads", k_genotype_reads},
                                                                  {"original_haplotype", k_original_haplotype},
                                                                  {"trimed_haplotype", k_trimed_haplotype},
                                                                  {"likelihoods", k_likelihoods}};

struct TestData
{
    pSimpleInterval original{nullptr};
    pSimpleInterval original_padded{nullptr};
    pSimpleInterval variant{nullptr};
    pSimpleInterval variant_padded{nullptr};
    pSimpleInterval ref_loc{nullptr};
    RefFragment ref_bases{0, nullptr};

    ReadVector original_reads;
    ReadVector trimed_reads;
    ReadVector filtered_reads1;
    ReadVector filtered_reads2;
    ReadVector pairhmm_reads;
    ReadVector genotype_reads;

    HaplotypeVector original_haplotype;
    HaplotypeVector trimed_haplotype;

    DoubleVector2D likelihoods;

    explicit TestData(int32_t original_reads_num, int32_t trimed_reads_num, int32_t filtered_reads1_num, int32_t pairhmm_reads_num,
                      int32_t filtered_reads2_num, int32_t genotype_reads_num, int32_t original_haplotype_num, int32_t trimed_haplotype_num,
                      pMemoryPool pool)
        : original_reads(pool)
        , trimed_reads(pool)
        , filtered_reads1(pool)
        , filtered_reads2(pool)
        , pairhmm_reads(pool)
        , genotype_reads(pool)
        , original_haplotype(pool)
        , trimed_haplotype(pool)
        , likelihoods(pool)
    {
        original_reads.reserve(original_reads_num);
        trimed_reads.reserve(trimed_reads_num);
        filtered_reads1.reserve(filtered_reads1_num);
        filtered_reads2.reserve(pairhmm_reads_num);
        pairhmm_reads.reserve(filtered_reads2_num);
        genotype_reads.reserve(genotype_reads_num);

        original_haplotype.reserve(original_haplotype_num);
        trimed_haplotype.reserve(trimed_haplotype_num);

        likelihoods.resize(trimed_haplotype_num);
        std::for_each(likelihoods.begin(), likelihoods.end(), [&](DoubleVector& dd) { dd.resize(genotype_reads_num); });
    }
};

#endif  // ROVACA_HC_UNIT_TEST_HEADER_H_
