#include <gtest/gtest.h>

#include <random>

#include "allele.h"
#include "forward.h"
#include "genotype.h"
#include "genotype_struct.h"
#include "genotypes_context.hpp"
#include "math_utils.h"
#include "simple_interval.h"
#include "unit_test_header.h"
#include "unit_test_utils.hpp"
#include "utils/debug_utils.h"
#include "utils/rovaca_variant_context_utils.h"
#include "variant.h"

using namespace std;
using namespace rovaca;

static constexpr size_t s_buffer_size = 1024 * 1024 * 10;

class RovacaVariantContextUtilsUnitTest : public ::testing::Test
{
public:
    pVariant makeVC(const AlleleVector& alleles, int64_t start);
    pVariant makeVC(int32_t source, const AlleleVector& alleles, pGenotypesContext genotypes);
    pVariant makeVC(pAllele ref, int64_t start) { return makeVC({{ref}, _pool}, start); }
    pVariant makeVC(pAllele ref, pAllele alt, int64_t start) { return makeVC({{ref, alt}, _pool}, start); }
    pVariant makeVC(int32_t source, const AlleleVector& alleles) { return makeVC(source, alleles, nullptr); }
    pVariant makeFromAlleles(int32_t source, int32_t tid, int64_t start, const std::vector<std::string>& alleleStrings);

    std::vector<std::tuple<VariantVector, pSimpleInterval, pAllele>> testDetermineReferenceAlleleData();
    std::pmr::vector<MergeAllelesTest> testMergeAllelesData();
    std::vector<MergeAllelesData> testMakeGenotypeCallData();
    static std::vector<ClipAllelesData> testClipAllelesData();
    std::vector<AlleleMap> testAlleleRemappingData();

    static void assertGenotypesAreEqual(pGenotype actual, pGenotype expected);

    static int32_t findNumberOfRepetitionsAdapter(const char* repeat_unit, const char* test_string, bool leading_repeats);

protected:
    void SetUp() override;
    void TearDown() override;

    uint8_t* _buffer{};
    pMemoryPool _pool{};

    pAllele Aref{}, Cref{}, Gref{}, Tref{}, A{}, T{}, C{}, G{}, ATC{}, ATCref{}, ATCATC{}, ATCATCT{}, ATref{}, Anoref{}, GT{};
    pSimpleInterval START_AT_1{}, START_AT_2{};
};

TEST_F(RovacaVariantContextUtilsUnitTest, testDetermineReferenceAllele)
{
    std::vector<std::tuple<VariantVector, pSimpleInterval, pAllele>> data = testDetermineReferenceAlleleData();
    for (size_t i = 0, len = data.size(); i < len; ++i) {
        const auto& tup = data.at(i);
        const VariantVector& vcs = get<0>(tup);
        pSimpleInterval loc = get<1>(tup);
        pAllele expectedRef = get<2>(tup);

        pAllele ref = ROVACAVariantContextUtils::determine_reference_allele(vcs, loc);

        if (expectedRef == nullptr) {
            ASSERT_TRUE(ref == nullptr) << i;
        }
        else {
            ASSERT_TRUE(ref->is_reference()) << i;
            ASSERT_EQ(DebugUtils::bases2str(ref->get_bases()), DebugUtils::bases2str(expectedRef->get_bases())) << i;
        }
    }
}

TEST_F(RovacaVariantContextUtilsUnitTest, testMergeAlleles)
{
    std::pmr::vector<MergeAllelesTest> data = testMergeAllelesData();

    for (size_t i = 0, len = data.size(); i < len; ++i) {
        const MergeAllelesTest& td = data.at(i);
        VariantVector inputs{_pool};
        int32_t source = 0;
        for (const AlleleVector& alleles : td.inputs) {
            inputs.push_back(makeVC(++source, alleles));
        }

        pVariant merged = ROVACAVariantContextUtils::simple_merge(inputs, FilteredRecordMergeType::KEEP_IF_ANY_UNFILTERED,
                                                                GenotypeMergeType::PRIORITIZE, false);

        ASSERT_EQ(merged->allele_num(), td.expected.size()) << i;
        for (size_t j = 0, allele_num = merged->allele_num(); j < allele_num; ++j) {
            std::string l = DebugUtils::bases2str(merged->alleles().at(j)->get_bases());
            std::string r = DebugUtils::bases2str(td.expected.at(j)->get_bases());
            ASSERT_EQ(l, r) << i << " " << j;
        }
    }
}

TEST_F(RovacaVariantContextUtilsUnitTest, testSubsetToRef)
{
    int32_t dp = 10, gq = 30;
    pGenotype l, r;
    AlleleVector as{{Aref, C}, _pool};
    std::pmr::vector<std::pmr::vector<pAllele>> allAlleles{{{Aref}, {C}, {Aref, C}, {Aref, C, C}}, _pool};
    for (const std::pmr::vector<pAllele>& alleles : allAlleles) {
        for (int32_t i : {1, 2}) {
            l = Genotype::create(_pool);
            l->set_id(i);
            l->set_alleles({alleles, _pool});
            l->set_dp(dp);
            l->set_gq(gq);
            l->set_ad(alleles.size() == 1 ? Int32Vector{{1}, _pool}
                                          : (alleles.size() == 2 ? Int32Vector{{1, 2}, _pool} : Int32Vector{{1, 2, 3}, _pool}));
            l->set_pl(alleles.size() == 1 ? Int32Vector{{1}, _pool}
                                          : (alleles.size() == 2 ? Int32Vector{{1, 2}, _pool} : Int32Vector{{1, 2, 3}, _pool}));

            pVariant vc = Variant::create(_pool);
            vc->set_source_id(0);
            vc->set_tid(20);
            vc->set_start(1);
            vc->set_stop(1);
            vc->set_alleles(as);
            pGenotypesContext gg = GenotypesContext::create(_pool);
            gg->add(l);
            vc->set_genotype(gg);

            // 实际就是根据ploidy将Genotype中的Allele设置为纯合子
            pGenotypesContext gc = ROVACAVariantContextUtils::subset_to_ref_only(vc, 2, _pool);

            AlleleVector refs{alleles.size(), Aref, _pool};
            r = Genotype::create(_pool);
            r->set_id(i);
            r->set_alleles(std::move(refs));
            r->set_dp(dp);
            r->set_gq(gq);
            r->set_ad(alleles.size() == 1 ? Int32Vector{{1}, _pool}
                                          : (alleles.size() == 2 ? Int32Vector{{1, 2}, _pool} : Int32Vector{{1, 2, 3}, _pool}));
            r->set_pl(alleles.size() == 1 ? Int32Vector{{1}, _pool}
                                          : (alleles.size() == 2 ? Int32Vector{{1, 2}, _pool} : Int32Vector{{1, 2, 3}, _pool}));

            ASSERT_EQ(gc->size(), 1) << i;
            assertGenotypesAreEqual(gc->at(0), r);
        }
    }
}

TEST_F(RovacaVariantContextUtilsUnitTest, testMakeGenotypeCall)
{
    std::vector<MergeAllelesData> data = testMakeGenotypeCallData();
    for (size_t i = 0, len = data.size(); i < len; ++i) {
        MergeAllelesData& d = data.at(i);

        ASSERT_EQ(d.ploidy_, d.originalGT_.size()) << i;
        pGenotype g = Genotype::create(_pool);

        AlleleVector originalGT{_pool};
        DoubleVector likelihoods = math_utils::normalize_log10(d.likelihoods_);
        ROVACAVariantContextUtils::make_genotype_call(d.ploidy_, g, d.mode_, likelihoods, d.allelesToUse_, originalGT, nullptr);

        ASSERT_EQ(g->alleles().size(), d.expectedAlleles_.size()) << i;
        for (size_t j = 0, allele_num = g->alleles().size(); j < allele_num; ++j) {
            std::string l = DebugUtils::bases2str(g->alleles().at(j)->get_bases());
            std::string r = DebugUtils::bases2str(d.expectedAlleles_.at(j)->get_bases());

            ASSERT_EQ(l, r) << i;
        }
    }
}

TEST_F(RovacaVariantContextUtilsUnitTest, testClipAlleles)
{
    std::vector<ClipAllelesData> data = testClipAllelesData();

    int64_t start = 10;
    for (size_t i = 0, len = data.size(); i < len; ++i) {
        const ClipAllelesData& d = data.at(i);
        int32_t allele_num = int32_t(d.alleleStrings.size());
        pVariant unclipped = makeFromAlleles(1, 20, start, d.alleleStrings);
        pVariant clipped = ROVACAVariantContextUtils::trim_alleles(unclipped, true, true);

        ASSERT_EQ(clipped->get_start(), start + d.numLeftClipped) << i;
        for (int32_t j = 0; j < allele_num; ++j) {
            std::string l = DebugUtils::bases2str(clipped->alleles().at(j)->get_display_string());
            ASSERT_EQ(l, d.expected.at(j)) << i;
        }
    }
}

TEST_F(RovacaVariantContextUtilsUnitTest, testAlleleRemapping)
{
    std::mt19937 gen{GATK_RANDOM_SEED};
    std::uniform_int_distribution<> dis(0, 1);

    std::vector<AlleleMap> data = testAlleleRemappingData();
    for (size_t i = 0, len = data.size(); i < len; ++i) {
        const AlleleMap& mm = data.at(i);
        AlleleVector alleleKeys{_pool};
        for (const auto& item : mm.old2new_map) {
            alleleKeys.push_back(item.first);
        }

        AlleleVector alleles{_pool};
        alleles.push_back(alleleKeys.at(dis(gen)));
        alleles.push_back(alleleKeys.at(dis(gen)));

        AlleleVector alleles_copy{alleles, _pool};

        pGenotype g = Genotype::create(_pool);
        g->set_alleles(std::move(alleles));

        pGenotypesContext originalGC = GenotypesContext::create(_pool);
        originalGC->add(g);

        pGenotypesContext remappedGC = ROVACAVariantContextUtils::update_genotypes_with_mapped_alleles(originalGC, mm);

        pGenotype originalG = originalGC->at(0);
        pGenotype remappedG = remappedGC->at(0);

        ASSERT_EQ(originalG->alleles().size(), remappedG->alleles().size()) << i;
        for (size_t j = 0, jlen = originalG->alleles().size(); j < jlen; ++j) {
            std::string l = DebugUtils::bases2str(remappedG->alleles().at(j)->get_display_string());
            std::string r = DebugUtils::bases2str(mm.old2new_map.at(alleles_copy.at(j))->get_display_string());
            ASSERT_EQ(l, r) << i << " " << j;
        }
    }
}

TEST_F(RovacaVariantContextUtilsUnitTest, testRepeatAllele)
{
    pAllele nullR = Aref;
    pAllele nullA = Allele::create_allele("A", false, _pool);
    pAllele atc = Allele::create_allele("AATC", false, _pool);
    pAllele atcatc = Allele::create_allele("AATCATC", false, _pool);
    pAllele ccccR = Allele::create_allele("ACCCC", true, _pool);
    pAllele cc = Allele::create_allele("ACC", false, _pool);
    pAllele cccccc = Allele::create_allele("ACCCCCC", false, _pool);
    pAllele gagaR = Allele::create_allele("AGAGA", true, _pool);
    pAllele gagagaga = Allele::create_allele("AGAGAGAGA", false, _pool);

    int64_t insLocStart = 20;
    int64_t insLocStop = 20;

    pVariant vc;
    auto result = new ALLOC_TYPE_IN_POOL(_pool, RepeatUnitsResult) RepeatUnitsResult{_pool};

    RefFragment ref1{(uint32_t)strlen("ATG"), (uint8_t*)"ATG"};
    RefFragment ref2{(uint32_t)strlen("AAA"), (uint8_t*)"AAA"};
    RefFragment ref3{(uint32_t)strlen("CACACAC"), (uint8_t*)"CACACAC"};
    RefFragment ref4{(uint32_t)strlen("CACACA"), (uint8_t*)"CACACA"};
    RefFragment ref5{(uint32_t)strlen("CATGCATG"), (uint8_t*)"CATGCATG"};
    RefFragment ref6{(uint32_t)strlen("AATAATA"), (uint8_t*)"AATAATA"};
    ASSERT_EQ(ROVACAVariantContextUtils::find_repeated_substring(&ref1), 3);
    ASSERT_EQ(ROVACAVariantContextUtils::find_repeated_substring(&ref2), 1);
    ASSERT_EQ(ROVACAVariantContextUtils::find_repeated_substring(&ref3), 7);
    ASSERT_EQ(ROVACAVariantContextUtils::find_repeated_substring(&ref4), 2);
    ASSERT_EQ(ROVACAVariantContextUtils::find_repeated_substring(&ref5), 4);
    ASSERT_EQ(ROVACAVariantContextUtils::find_repeated_substring(&ref6), 7);

    RefFragment refBytes{(uint32_t)strlen("ATCATCATCGGA"), (uint8_t*)"ATCATCATCGGA"};

    // A*,ATC, context = ATC ATC ATC : (ATC)3 -> (ATC)4
    result->lengths.clear();
    result->bases.len = 0;
    result->bases.data = nullptr;
    vc = Variant::create(_pool);
    vc->set_source_id(0);
    vc->set_tid(0);
    vc->set_start(insLocStart);
    vc->set_stop(insLocStop);
    vc->set_alleles({{nullR, atc}, _pool});
    ROVACAVariantContextUtils::get_num_tandem_repeat_units(vc, &refBytes, result, _pool);
    ASSERT_EQ(result->lengths[0], 3);
    ASSERT_EQ(result->lengths[1], 4);
    ASSERT_EQ(result->bases.len, 3);

    // ATC*,A,ATCATC
    result->lengths.clear();
    result->bases.len = 0;
    result->bases.data = nullptr;
    vc = Variant::create(_pool);
    vc->set_source_id(0);
    vc->set_tid(0);
    vc->set_start(insLocStart);
    vc->set_stop(insLocStart + 3);
    vc->set_alleles({{Allele ::create_allele("AATC", 1, _pool), nullA, atcatc}, _pool});
    ROVACAVariantContextUtils::get_num_tandem_repeat_units(vc, &refBytes, result, _pool);
    ASSERT_EQ(result->lengths[0], 3);
    ASSERT_EQ(result->lengths[1], 2);
    ASSERT_EQ(result->lengths[2], 4);
    ASSERT_EQ(result->bases.len, 3);

    // simple non-tandem deletion: CCCC*, -
    result->lengths.clear();
    result->bases.len = 0;
    result->bases.data = nullptr;
    refBytes.data = (uint8_t*)"CCCCCCCCATG";
    refBytes.len = strlen("CCCCCCCCATG");
    vc = Variant::create(_pool);
    vc->set_source_id(0);
    vc->set_tid(0);
    vc->set_start(10);
    vc->set_stop(14);
    vc->set_alleles({{ccccR, nullA}, _pool});
    ROVACAVariantContextUtils::get_num_tandem_repeat_units(vc, &refBytes, result, _pool);
    ASSERT_EQ(result->lengths[0], 8);
    ASSERT_EQ(result->lengths[1], 4);
    ASSERT_EQ(result->bases.len, 1);

    // CCCC*,CC,-,CCCCCC, context = CCC: (C)7 -> (C)5,(C)3,(C)9
    result->lengths.clear();
    result->bases.len = 0;
    result->bases.data = nullptr;
    refBytes.data = (uint8_t*)"CCCCCCCAGAGAGAG";
    refBytes.len = strlen("CCCCCCCAGAGAGAG");
    vc = Variant::create(_pool);
    vc->set_source_id(0);
    vc->set_tid(0);
    vc->set_start(insLocStart);
    vc->set_stop(insLocStart + 4);
    vc->set_alleles({{ccccR, cc, nullA, cccccc}, _pool});
    ROVACAVariantContextUtils::get_num_tandem_repeat_units(vc, &refBytes, result, _pool);
    ASSERT_EQ(result->lengths[0], 7);
    ASSERT_EQ(result->lengths[1], 5);
    ASSERT_EQ(result->lengths[2], 3);
    ASSERT_EQ(result->lengths[3], 9);
    ASSERT_EQ(result->bases.len, 1);

    // GAGA*,-,GAGAGAGA
    result->lengths.clear();
    result->bases.len = 0;
    result->bases.data = nullptr;
    refBytes.data = (uint8_t*)"GAGAGAGAGATTT";
    refBytes.len = strlen("GAGAGAGAGATTT");
    vc = Variant::create(_pool);
    vc->set_source_id(0);
    vc->set_tid(0);
    vc->set_start(insLocStart);
    vc->set_stop(insLocStart + 4);
    vc->set_alleles({{gagaR, nullA, gagagaga}, _pool});
    ROVACAVariantContextUtils::get_num_tandem_repeat_units(vc, &refBytes, result, _pool);
    ASSERT_EQ(result->lengths[0], 5);
    ASSERT_EQ(result->lengths[1], 3);
    ASSERT_EQ(result->lengths[2], 7);
    ASSERT_EQ(result->bases.len, 2);
}

TEST_F(RovacaVariantContextUtilsUnitTest, testFindNumberOfRepetitions)
{
    ASSERT_EQ(findNumberOfRepetitionsAdapter("AT", "GATAT", false), 2);
    ASSERT_EQ(findNumberOfRepetitionsAdapter("AT", "GATAT", true), 0);
    ASSERT_EQ(findNumberOfRepetitionsAdapter("A", "ATATG", true), 1);
    ASSERT_EQ(findNumberOfRepetitionsAdapter("AT", "ATATG", true), 2);
    ASSERT_EQ(findNumberOfRepetitionsAdapter("CCC", "CCCCCCCC", true), 2);
    ASSERT_EQ(findNumberOfRepetitionsAdapter("CCC", "CCCCCCCC", false), 2);

    ASSERT_EQ(findNumberOfRepetitionsAdapter("ATG", "ATGATGATGATG", true), 4);
    ASSERT_EQ(findNumberOfRepetitionsAdapter("G", "ATGATGATGATG", true), 0);
    ASSERT_EQ(findNumberOfRepetitionsAdapter("T", "T", true), 1);
    ASSERT_EQ(findNumberOfRepetitionsAdapter("AT", "ATGATGATCATG", true), 1);
    ASSERT_EQ(findNumberOfRepetitionsAdapter("CCCCCCCC", "CCC", true), 0);
    ASSERT_EQ(findNumberOfRepetitionsAdapter("AT", "AT", true), 1);
    ASSERT_EQ(findNumberOfRepetitionsAdapter("AT", "", true), 0);

    ASSERT_EQ(findNumberOfRepetitionsAdapter("ATG", "ATGATGATGATG", false), 4);
    ASSERT_EQ(findNumberOfRepetitionsAdapter("G", "ATGATGATGATG", false), 1);
    ASSERT_EQ(findNumberOfRepetitionsAdapter("T", "T", false), 1);
    ASSERT_EQ(findNumberOfRepetitionsAdapter("AT", "ATGATGATCATG", false), 0);
    ASSERT_EQ(findNumberOfRepetitionsAdapter("CCCCCCCC", "CCC", false), 0);
    ASSERT_EQ(findNumberOfRepetitionsAdapter("AT", "AT", false), 1);
    ASSERT_EQ(findNumberOfRepetitionsAdapter("AT", "", false), 0);  // empty test string
}

void RovacaVariantContextUtilsUnitTest::SetUp()
{
    _buffer = new uint8_t[s_buffer_size]{};
    _pool = new std::pmr::monotonic_buffer_resource(_buffer, s_buffer_size, std::pmr::null_memory_resource());

    Aref = Allele::create_allele('A', 1);
    Cref = Allele::create_allele('C', 1);
    Gref = Allele::create_allele('G', 1);
    Tref = Allele::create_allele('T', 1);
    A = Allele::create_allele('A', 0);
    T = Allele::create_allele('T', 0);
    C = Allele::create_allele('C', 0);
    G = Allele::create_allele('G', 0);
    ATC = Allele::create_allele("ATC", 0, _pool);
    ATCref = Allele::create_allele("ATC", 1, _pool);
    ATCATC = Allele::create_allele("ATCATC", 0, _pool);
    ATCATCT = Allele::create_allele("ATCATCT", 0, _pool);
    ATref = Allele::create_allele("AT", 1, _pool);
    Anoref = Allele::create_allele('A', 0);
    GT = Allele::create_allele("GT", 0, _pool);

    START_AT_1 = SimpleInterval::create(0, 1, 1, _pool);
    START_AT_2 = SimpleInterval::create(0, 2, 2, _pool);
}

void RovacaVariantContextUtilsUnitTest::TearDown()
{
    delete _pool;
    delete[] _buffer;
}

pVariant RovacaVariantContextUtilsUnitTest::makeVC(const AlleleVector& alleles, int64_t start)
{
    int64_t stop = start + alleles.at(0)->length() - 1;
    pVariant ret = Variant::create(_pool);
    ret->set_source_id(0);
    ret->set_tid(1);
    ret->set_start(start);
    ret->set_stop(stop);
    ret->set_alleles(alleles);
    return ret;
}

pVariant RovacaVariantContextUtilsUnitTest::makeVC(int32_t source, const AlleleVector& alleles, pGenotypesContext genotypes)
{
    int64_t start = 10;
    int64_t stop = start + alleles.at(0)->length() - 1;
    pVariant ret = Variant::create(_pool);
    ret->set_source_id(source);
    ret->set_tid(1);
    ret->set_start(start);
    ret->set_stop(stop);
    ret->set_alleles(alleles);
    if (nullptr != genotypes) {
        ret->set_genotype(genotypes);
    }
    return ret;
}

std::vector<std::tuple<VariantVector, pSimpleInterval, pAllele>> RovacaVariantContextUtilsUnitTest::testDetermineReferenceAlleleData()
{
    return {{{}, nullptr, nullptr},
            {{}, START_AT_1, nullptr},
            {{{makeVC(Aref, 1)}, _pool}, START_AT_1, Aref},
            {{{makeVC(Aref, 1)}, _pool}, START_AT_2, nullptr},
            {{{makeVC(Aref, 1), makeVC(ATref, 1)}, _pool}, START_AT_1, ATref},
            {{{makeVC(ATref, 1), makeVC(Aref, 1)}, _pool}, START_AT_1, ATref},
            {{{makeVC(Aref, 1), makeVC(ATref, 1)}, _pool}, nullptr, ATref},
            {{{makeVC(ATref, 1), makeVC(Aref, 1)}, _pool}, START_AT_2, nullptr},
            {{{makeVC(Aref, 1), makeVC(ATref, 2)}, _pool}, START_AT_1, Aref},
            {{{makeVC(Aref, C, 1), makeVC(ATref, ATCATC, 1)}, _pool}, START_AT_1, ATref},
            {{{makeVC(Aref, 1), makeVC(ATCref, 1), makeVC(ATref, 1)}, _pool}, START_AT_1, ATCref}};
}

std::pmr::vector<MergeAllelesTest> RovacaVariantContextUtilsUnitTest::testMergeAllelesData()
{
    std::pmr::vector<MergeAllelesTest> ret{_pool};
    ret.push_back({{{Aref}}, {Aref}, _pool});
    ret.push_back({{{Aref}, {Aref}}, {Aref}, _pool});
    ret.push_back({{{Aref}, {Aref, T}}, {Aref, T}, _pool});
    ret.push_back({{{Aref, C}, {Aref, T}}, {Aref, C, T}, _pool});
    ret.push_back({{{Aref, T}, {Aref, C}}, {Aref, T, C}, _pool});
    ret.push_back({{{Aref, C, T}, {Aref, C}}, {Aref, C, T}, _pool});
    ret.push_back({{{Aref, C, T}}, {Aref, C, T}, _pool});
    ret.push_back({{{Aref, T, C}}, {Aref, T, C}, _pool});
    ret.push_back({{{Aref, T, C}, {Aref, C}}, {Aref, T, C}, _pool});
    ret.push_back({{{Aref}, {Aref, ATC}}, {Aref, ATC}, _pool});
    ret.push_back({{{Aref}, {Aref, ATC, ATCATC}}, {Aref, ATC, ATCATC}, _pool});
    ret.push_back({{{Aref, ATCATC}, {Aref, ATC, ATCATC}}, {Aref, ATCATC, ATC}, _pool});
    ret.push_back({{{Aref, ATC}, {Aref, ATCATC}}, {Aref, ATC, ATCATC}, _pool});
    ret.push_back({{{ATref, ATC, Anoref, G}, {Aref, ATCATC, G}}, {ATref, ATC, Anoref, G, ATCATCT, GT}, _pool});
    return ret;
}

void RovacaVariantContextUtilsUnitTest::assertGenotypesAreEqual(pGenotype actual, pGenotype expected)
{
    ASSERT_EQ(actual->sample_id(), expected->sample_id());
    ASSERT_EQ(actual->alleles().size(), expected->alleles().size());
    for (size_t i = 0, len = actual->alleles().size(); i < len; ++i) {
        std::string l = DebugUtils::bases2str(actual->alleles().at(i)->get_bases());
        std::string r = DebugUtils::bases2str(expected->alleles().at(i)->get_bases());
        ASSERT_EQ(l, r);
    }
    ASSERT_EQ(actual->get_type(), expected->get_type());

    ASSERT_EQ(actual->has_ad(), expected->has_ad());
    ASSERT_EQ(actual->get_dp(), expected->get_dp());
    ASSERT_EQ(actual->has_ad(), expected->has_ad());
    ASSERT_EQ(actual->ad(), expected->ad());
    ASSERT_EQ(actual->has_gq(), expected->has_gq());
    ASSERT_EQ(actual->get_gq(), expected->get_gq());
    ASSERT_EQ(actual->has_likelihoods(), expected->has_likelihoods());
    ASSERT_EQ(actual->pl(), expected->pl());

    ASSERT_EQ(actual->is_phased(), expected->is_phased());
    ASSERT_EQ(actual->get_ploidy(), expected->get_ploidy());
}

std::vector<MergeAllelesData> RovacaVariantContextUtilsUnitTest::testMakeGenotypeCallData()
{
    std::vector<MergeAllelesData> ret;

    std::vector<pAllele> AA{Aref, Aref};
    std::vector<pAllele> AC{Aref, C};
    std::vector<pAllele> CC{C, C};
    std::vector<pAllele> AG{Aref, G};
    std::vector<pAllele> CG{C, G};
    std::vector<pAllele> GG{G, G};
    std::vector<pAllele> ACG{Aref, C, G};
    std::vector<std::vector<pAllele>> allDiploidSubsetAlleles{AC, AG, ACG};

    // for P=2 and N=1, the ordering is 00,01,11
    std::vector<double> homRefPL{0.9, 0.09, 0.01};
    std::vector<double> hetPL{0.09, 0.9, 0.01};
    std::vector<double> homVarPL{0.01, 0.09, 0.9};
    std::vector<double> uninformative{0.33, 0.33, 0.33};
    std::vector<std::vector<double>> allDiploidPLs{homRefPL, hetPL, homVarPL, uninformative};

    std::vector<pAllele> DiploidNocalls{StaticAllele::get_instance()->_no_call.get(), StaticAllele::get_instance()->_no_call.get()};

    for (const std::vector<pAllele>& alleles : allDiploidSubsetAlleles) {
        ret.emplace_back(_pool, 2, GenotypeAssignmentMethod::SET_TO_NO_CALL, allDiploidPLs.at(0), AA, alleles, DiploidNocalls);
    }

    for (const std::vector<pAllele>& originalGT : {AA, AC, CC, AG, CG, GG}) {
        ret.emplace_back(_pool, 2, GenotypeAssignmentMethod::USE_PLS_TO_ASSIGN, homRefPL, originalGT, AC, AA);
        ret.emplace_back(_pool, 2, GenotypeAssignmentMethod::USE_PLS_TO_ASSIGN, hetPL, originalGT, AC, AC);
        ret.emplace_back(_pool, 2, GenotypeAssignmentMethod::USE_PLS_TO_ASSIGN, homVarPL, originalGT, AC, CC);
    }

    return ret;
}

std::vector<ClipAllelesData> RovacaVariantContextUtilsUnitTest::testClipAllelesData()
{
    return {{{"ACC", "AC", "<NON_REF>"}, {"AC", "A", "<NON_REF>"}, 0},
            {{"ACC", "AC", "*"}, {"AC", "A", "*"}, 0},
            {{"ACC", "AC"}, {"AC", "A"}, 0},
            {{"ACGC", "ACG"}, {"GC", "G"}, 2},
            {{"ACGC", "ACGA"}, {"C", "A"}, 3},
            {{"ACGC", "AGC"}, {"AC", "A"}, 0},
            {{"AT", "AC", "AG"}, {"T", "C", "G"}, 1},
            {{"AT", "AC", "ACG"}, {"T", "C", "CG"}, 1},
            {{"AC", "ACT", "ACG"}, {"C", "CT", "CG"}, 1},
            {{"ACG", "ACGT", "ACGTA"}, {"G", "GT", "GTA"}, 2},
            {{"ACG", "ACGT", "ACGCA"}, {"G", "GT", "GCA"}, 2},
            {{"ACGTT", "ACCTT"}, {"G", "C"}, 2},
            {{"ACGTT", "ACCCTT"}, {"G", "CC"}, 2},
            {{"ACGTT", "ACGCTT"}, {"G", "GC"}, 2}};
}

pVariant RovacaVariantContextUtilsUnitTest::makeFromAlleles(int32_t source, int32_t tid, int64_t start,
                                                          const vector<std::string>& alleleStrings)
{
    uint8_t first = 1;
    AlleleVector alleles{_pool};
    size_t length = alleleStrings.at(0).size();
    for (const std::string& alleleString : alleleStrings) {
        alleles.push_back(Allele::create_allele(alleleString.c_str(), first, _pool));
        first = 0;
    }

    pVariant ret = Variant::create(_pool);
    ret->set_source_id(source);
    ret->set_tid(tid);
    ret->set_start(start);
    ret->set_stop(start + int64_t(length) - 1);
    ret->set_alleles(alleles);

    return ret;
}

std::vector<AlleleMap> RovacaVariantContextUtilsUnitTest::testAlleleRemappingData()
{
    std::vector<AlleleMap> ret;
    pAllele originalBase1 = Allele::create_allele('A', 0);
    pAllele originalBase2 = Allele::create_allele('T', 0);
    pAllele a1, a2;
    for (char bases1 : {'A', 'C', 'G', 'T'}) {
        for (char bases2 : {'A', 'C', 'G', 'T'}) {
            AlleleMap mm{_pool};
            a1 = Allele::create_allele(bases1, 0);
            a2 = Allele::create_allele(bases2, 0);
            mm.new_arr.push_back(a1);
            mm.new_arr.push_back(a2);
            mm.old2new_map.insert({originalBase1, a1});
            mm.old2new_map.insert({originalBase2, a2});
            ret.push_back(std::move(mm));
        }
    }
    return ret;
}

int32_t RovacaVariantContextUtilsUnitTest::findNumberOfRepetitionsAdapter(const char* repeat_unit, const char* test_string,
                                                                        bool leading_repeats)
{
    RefFragment ref{(uint32_t)strlen(repeat_unit), (uint8_t*)repeat_unit};
    RefFragment test{(uint32_t)strlen(test_string), (uint8_t*)test_string};
    return ROVACAVariantContextUtils::find_number_of_repetitions(&ref, &test, leading_repeats);
}
