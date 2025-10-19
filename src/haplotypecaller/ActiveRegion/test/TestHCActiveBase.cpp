

#include <boost/dynamic_bitset.hpp>
#include <memory_resource>
#include <vector>

#include "HcActiveBase.h"
#include "gtest/gtest.h"
#include "htslib/faidx.h"
#include "htslib/sam.h"
#include "math_utils.h"
#include "quality_utils.h"

using namespace rovaca;
using namespace math_utils;
using namespace QualityUtils;
static constexpr size_t s_max_bufer_size = 1024 * 1024 * 20;

typedef std::pmr::memory_resource MemoryPool, *pMemoryPool;

class HCLikelihoodWGS : public ::testing::Test
{
protected:
    static uint8_t* _buffer;
    static pMemoryPool _pool;
    static char* ref;
    static int ref_len;
    static void SetUpTestSuite()
    {
        faidx_t* fai = fai_load("/data/pipelines/WGS_bgionline/db/GRCh37/ref/GRCh37_no_alt.fna");
        ref_len = faidx_seq_len(fai, "chr22");
        int len = 0;
        ref = faidx_fetch_seq(fai, "chr22", 0, ref_len, &len);

        _buffer = new uint8_t[s_max_bufer_size]{};
        _pool = new std::pmr::monotonic_buffer_resource(_buffer, s_max_bufer_size, std::pmr::null_memory_resource());
        fai_destroy(fai);
        QualityUtils::init();
    }

    static void TearDownTestSuite()
    {
        delete _pool;
        delete[] _buffer;
    }
    void SetUp() override { _bit = new boost::dynamic_bitset<>(1024, 0x0); }

    void TearDown() override { delete _bit; }

    boost::dynamic_bitset<>* _bit;
};

uint8_t* HCLikelihoodWGS::_buffer = nullptr;
pMemoryPool HCLikelihoodWGS::_pool = nullptr;
int HCLikelihoodWGS::ref_len = 0;
char* HCLikelihoodWGS::ref = nullptr;

TEST_F(HCLikelihoodWGS, common_func_test0)
{
    HcActiveBase base(_pool, 2, 0, 0, 1024, _bit, 10, NULL, 0);

    EXPECT_EQ(1024, base.get_size());
    EXPECT_EQ(1024, base.get_stop());
    EXPECT_EQ(10, base.get_work_id());
    EXPECT_EQ(base.get_offset(0, 1), -1);
}

TEST_F(HCLikelihoodWGS, common_func_test1)
{
    for (int i = 0; i < 924; i++) _bit->set(i, true);

    HcActiveBase base(_pool, 2, 0, 100, 1024, _bit, 10, NULL, 0);

    EXPECT_EQ(924, base.get_size());
    EXPECT_EQ(1024, base.get_stop());
    EXPECT_EQ(10, base.get_work_id());
    EXPECT_EQ(base.get_offset(0, 100), 0);
    EXPECT_EQ(base.get_offset(0, 99), -1);
    EXPECT_EQ(base.get_offset(1, 100), -1);
    EXPECT_EQ(base.get_offset(0, 1024), -1);
    EXPECT_EQ(base.get_offset(0, 1025), -1);
}

TEST_F(HCLikelihoodWGS, slot_likelihood_case)
{
    for (int i = 0; i < 924; i++) _bit->set(i, true);

    HcActiveBase base(_pool, 2, 0, 100, 1024, _bit, 10, NULL, 0);

    std::vector<std::pair<HCActiveBaseStatus, uint8_t>> inputbase({{HCActiveBaseStatus::NONREF, 32}});
    std::vector<double> expected_values({-2.741077772782652E-4, -0.3012127050374378, -3.677121254719663});
    for (auto input_pair : inputbase) {
        base.increase_hist(0, input_pair.first, input_pair.second);
    }
    base.compute_genotype_PL(0);
    auto likelihoods = base.get_current_likelihood();
    EXPECT_EQ(likelihoods.size(), expected_values.size());

    for (int i = 0; i < 3; i++) {
        EXPECT_NEAR(likelihoods[i], expected_values[i], 1E-8);
    }
}

TEST_F(HCLikelihoodWGS, final_likelihood)
{
    for (int i = 0; i < 924; i++) _bit->set(i, true);

    HcActiveBase base(_pool, 2, 0, 100, 1024, _bit, 10, NULL, 0);

    std::vector<std::pair<HCActiveBaseStatus, uint8_t>> inputbase({{HCActiveBaseStatus::NONREF, 41},
                                                                   {HCActiveBaseStatus::NONREF, 41},
                                                                   {HCActiveBaseStatus::NONREF, 41},
                                                                   {HCActiveBaseStatus::VariantIDX, 22},
                                                                   {HCActiveBaseStatus::NONREF, 41},
                                                                   {HCActiveBaseStatus::NONREF, 37},
                                                                   {HCActiveBaseStatus::NONREF, 41},
                                                                   {HCActiveBaseStatus::NONREF, 41},
                                                                   {HCActiveBaseStatus::NONREF, 32},
                                                                   {HCActiveBaseStatus::NONREF, 41},
                                                                   {HCActiveBaseStatus::VariantIDX, 41},
                                                                   {HCActiveBaseStatus::VariantIDX, 41},
                                                                   {HCActiveBaseStatus::VariantIDX, 12},
                                                                   {HCActiveBaseStatus::VariantIDX, 41}

    });
    std::vector<double> expected_values({-18.08620853341253, -4.235384193484483, -39.92524846633521});
    for (auto input_pair : inputbase) {
        base.increase_hist(0, input_pair.first, input_pair.second);
    }
    base.compute_genotype_PL(0);
    auto likelihoods = base.get_current_likelihood();
    EXPECT_EQ(likelihoods.size(), expected_values.size());
    for (int i = 0; i < 3; i++) {
        EXPECT_NEAR(likelihoods[i], expected_values[i], 1E-8);
    }
    double ative_value = base.compute_biallelic_non_ref_posterior();
    EXPECT_NEAR(0.9999999999930759, ative_value, 1E-8);
}

TEST_F(HCLikelihoodWGS, HQ_CASE1)
{
    for (int i = 0; i < 924; i++) _bit->set(i, true);
    HcActiveBase base(_pool, 2, 0, 100, 1024, _bit, 10, NULL, 0);
    base.insert_high_count(0, 10);
    base.insert_high_count(0, 2);
    base.insert_high_count(0, 6);
    base.compute_extension_length(0, 1.0);

    auto result = base.get_active_result(0);

    EXPECT_EQ(result.acitve_value, 1.0);
    EXPECT_EQ(result.state, BaseState::HC_REGION_BASES_ACTIVE_ONE_SLOT);
    EXPECT_EQ(result.extend_length, 0);
}

TEST_F(HCLikelihoodWGS, HQ_CASE2)
{
    for (int i = 0; i < 924; i++) _bit->set(i, true);
    HcActiveBase base(_pool, 2, 0, 100, 1024, _bit, 10, NULL, 0);
    base.insert_high_count(0, 10);
    base.insert_high_count(0, 4);
    base.insert_high_count(0, 7);
    base.compute_extension_length(0, 1.0);

    auto result = base.get_active_result(0);

    EXPECT_EQ(result.acitve_value, 1.0);
    EXPECT_EQ(result.state, BaseState::HC_REGION_BASES_ACTIVE_HQ_SOFT_CLIPS);
    EXPECT_EQ(result.extend_length, 7);
}

TEST_F(HCLikelihoodWGS, HQ_CASE3)
{
    for (int i = 0; i < 924; i++) _bit->set(i, true);
    HcActiveBase base(_pool, 2, 0, 100, 1024, _bit, 10, NULL, 0);
    base.insert_high_count(0, 10);
    base.insert_high_count(0, 4);
    base.insert_high_count(0, 7);
    base.compute_extension_length(0, 0);

    auto result = base.get_active_result(0);

    EXPECT_EQ(result.acitve_value, 0);
    EXPECT_EQ(result.state, BaseState::HC_REGION_BASES_ACTIVE_NONE);
    EXPECT_EQ(result.extend_length, 0);
}
TEST_F(HCLikelihoodWGS, HQ_CASE4)
{
    for (int i = 0; i < 924; i++) _bit->set(i, true);
    HcActiveBase base(_pool, 2, 0, 100, 1024, _bit, 10, NULL, 0);
    base.insert_high_count(0, 10);
    base.insert_high_count(0, 5);
    base.compute_extension_length(0, 1.0);

    auto result = base.get_active_result(0);

    EXPECT_EQ(result.acitve_value, 1.0);
    EXPECT_EQ(result.state, BaseState::HC_REGION_BASES_ACTIVE_HQ_SOFT_CLIPS);
    EXPECT_EQ(result.extend_length, 7);
}

TEST_F(HCLikelihoodWGS, HQ_CASE5)
{
    for (int i = 0; i < 924; i++) _bit->set(i, true);
    HcActiveBase base(_pool, 2, 0, 100, 1024, _bit, 10, NULL, 0);
    base.insert_high_count(0, 100);
    base.insert_high_count(0, 101);
    base.compute_extension_length(0, 1.0);

    auto result = base.get_active_result(0);

    EXPECT_EQ(result.acitve_value, 1.0);
    EXPECT_EQ(result.state, BaseState::HC_REGION_BASES_ACTIVE_HQ_SOFT_CLIPS);
    EXPECT_EQ(result.extend_length, 100);
}

TEST_F(HCLikelihoodWGS, High_QUAL_CLIPS)
{
    const char sam[] =
        "data:,"
        "@SQ\tSN:chr1\tLN:249250621\n"
        "1\t129\tchr1\t10000\t0\t10S10M\t*\t0\t0\tAAAAAAAAAAAAAAAAAAAA\tAAAAAAAAAAAAAAAAAAAA\n"
        "1\t129\tchr1\t10000\t0\t10S5M5S\t*\t0\t0\tAAAAAAAAAAAAAAAAAAAA\tAAAAAAAAAAAAAAAAAAAA\n"
        "1\t129\tchr1\t10000\t0\t10S5M5S\t*\t0\t0\tAAAAAAAAAAAAAAAAAAAA\t!AAAAAAAAAAAAAAAAAAA\n"
        "1\t129\tchr1\t10000\t0\t10S5M5S\t*\t0\t0\tAAAAAAAAAAAAAAAAAAAA\t!AAAAAAAAAAAAA!AAAAA\n";

    samFile* in = sam_open(sam, "r");
    bam_hdr_t* h = sam_hdr_read(in);
    bam1_t* bam = bam_init1();

    EXPECT_TRUE(sam_read1(in, h, bam) >= 0);
    EXPECT_EQ(HcActiveBase::get_high_quality_soft_clips(bam), 10);
    EXPECT_TRUE(sam_read1(in, h, bam) >= 0);
    EXPECT_EQ(HcActiveBase::get_high_quality_soft_clips(bam), 15);
    EXPECT_TRUE(sam_read1(in, h, bam) >= 0);
    EXPECT_EQ(HcActiveBase::get_high_quality_soft_clips(bam), 14);
    EXPECT_TRUE(sam_read1(in, h, bam) >= 0);
    EXPECT_EQ(HcActiveBase::get_high_quality_soft_clips(bam), 14);
    bam_destroy1(bam);
    bam_hdr_destroy(h);
    sam_close(in);
}

TEST_F(HCLikelihoodWGS, GET_Adaptor)
{
    const char sam[] =
        "data:,"
        "@SQ\tSN:chr1\tLN:249250621\n"
        "1\t35\tchr1\t998\t60\t8M\t=\t1000\t10\tACGTACGT\tACGTACGT\n"
        "1\t35\tchr1\t1002\t60\t8M\t=\t1000\t10\tACGTACGT\tACGTACGT\n"
        "1\t19\tchr1\t1002\t60\t8M\t=\t1000\t10\tACGTACGT\tACGTACGT\n"
        "1\t19\tchr1\t998\t60\t8M\t=\t1000\t10\tACGTACGT\tACGTACGT\n"
        "1\t19\tchr1\t998\t60\t8M\t=\t1000\t0\tACGTACGT\tACGTACGT\n"
        "1\t35\tchr1\t998\t60\t8M\t=\t1000\t0\tACGTACGT\tACGTACGT\n"
        "1\t7\tchr1\t998\t60\t8M\t=\t1000\t0\tACGTACGT\tACGTACGT\n"
        "1\t19\tchr1\t980\t60\t8M\t=\t1000\t20\tACGTACGT\tACGTACGT\n"
        "1\t1\tchr1\t10000\t60\t8M\t=\t1000\t20\tACGTACGT\tACGTACGT\n"
        "1\t35\tchr1\t998\t60\t8M\t=\t1000\t10\tACGTACGT\tACGTACGT\n"
        "1\t19\tchr1\t998\t60\t8M\t=\t1000\t10\tACGTACGT\tACGTACGT\n"
        "1\t3\tchr1\t998\t60\t8M\t=\t1000\t10\tACGTACGT\tACGTACGT\n"
        "1\t51\tchr1\t998\t60\t8M\t=\t1000\t10\tACGTACGT\tACGTACGT\n";

    ;

    samFile* in = sam_open(sam, "r");
    bam_hdr_t* h = sam_hdr_read(in);
    bam1_t* bam = bam_init1();
    EXPECT_TRUE(sam_read1(in, h, bam) >= 0);
    EXPECT_EQ(HcActiveBase::get_adaptor_boundary(bam), 998 + 10);
    EXPECT_TRUE(sam_read1(in, h, bam) >= 0);
    EXPECT_EQ(HcActiveBase::get_adaptor_boundary(bam), 1002 + 10);
    EXPECT_TRUE(sam_read1(in, h, bam) >= 0);
    EXPECT_EQ(HcActiveBase::get_adaptor_boundary(bam), 999);
    EXPECT_TRUE(sam_read1(in, h, bam) >= 0);
    EXPECT_EQ(HcActiveBase::get_adaptor_boundary(bam), 999);
    EXPECT_TRUE(sam_read1(in, h, bam) >= 0);
    EXPECT_EQ(HcActiveBase::get_adaptor_boundary(bam), INT_MIN);
    EXPECT_TRUE(sam_read1(in, h, bam) >= 0);
    EXPECT_EQ(HcActiveBase::get_adaptor_boundary(bam), INT_MIN);
    EXPECT_TRUE(sam_read1(in, h, bam) >= 0);
    EXPECT_EQ(HcActiveBase::get_adaptor_boundary(bam), INT_MIN);
    EXPECT_TRUE(sam_read1(in, h, bam) >= 0);
    EXPECT_EQ(HcActiveBase::get_adaptor_boundary(bam), INT_MIN);
    EXPECT_TRUE(sam_read1(in, h, bam) >= 0);
    EXPECT_EQ(HcActiveBase::get_adaptor_boundary(bam), INT_MIN);
    EXPECT_TRUE(sam_read1(in, h, bam) >= 0);
    EXPECT_NE(HcActiveBase::get_adaptor_boundary(bam), INT_MIN);
    EXPECT_TRUE(sam_read1(in, h, bam) >= 0);
    EXPECT_NE(HcActiveBase::get_adaptor_boundary(bam), INT_MIN);
    EXPECT_TRUE(sam_read1(in, h, bam) >= 0);
    EXPECT_EQ(HcActiveBase::get_adaptor_boundary(bam), INT_MIN);
    EXPECT_TRUE(sam_read1(in, h, bam) >= 0);
    EXPECT_EQ(HcActiveBase::get_adaptor_boundary(bam), INT_MIN);

    bam_destroy1(bam);
    bam_hdr_destroy(h);
    sam_close(in);
}

TEST_F(HCLikelihoodWGS, process_bam_to_slot_case1)
{
    const char sam[] =
        "data:,"
        "@SQ\tSN:chr1\tLN:249250621\n"
        "@SQ\tSN:chr2\tLN:243199373\n"
        "@SQ\tSN:chr3\tLN:198022430\n"
        "@SQ\tSN:chr4\tLN:191154276\n"
        "@SQ\tSN:chr5\tLN:180915260\n"
        "@SQ\tSN:chr6\tLN:171115067\n"
        "@SQ\tSN:chr7\tLN:159138663\n"
        "@SQ\tSN:chr8\tLN:146364022\n"
        "@SQ\tSN:chr9\tLN:141213431\n"
        "@SQ\tSN:chr10\tLN:135534747\n"
        "@SQ\tSN:chr11\tLN:135006516\n"
        "@SQ\tSN:chr12\tLN:133851895\n"
        "@SQ\tSN:chr13\tLN:115169878\n"
        "@SQ\tSN:chr14\tLN:107349540\n"
        "@SQ\tSN:chr15\tLN:102531392\n"
        "@SQ\tSN:chr16\tLN:90354753\n"
        "@SQ\tSN:chr17\tLN:81195210\n"
        "@SQ\tSN:chr18\tLN:78077248\n"
        "@SQ\tSN:chr19\tLN:59128983\n"
        "@SQ\tSN:chr20\tLN:63025520\n"
        "@SQ\tSN:chr21\tLN:48129895\n"
        "@SQ\tSN:chr22\tLN:51304566\n"
        "SRR14724454.133886660\t163\tchr22\t16050006\t27\t150M\t="
        "\t16050296\t440\tGATAAGTCCCAGGACTTCAGAAGAGCTGTGCGACCTTGGCCAAGTCACTTCCTCCTTCAGGAACATTGCAGTGGGCCTAAGTGCCTCCTCTCGGGACTGGTATGGGGACGGTC"
        "ATGCAATCTGGACAACATTCACCTTTAAAAGTTTATT\tAAFFFJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJAJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJ"
        "JJJJJJFJJJJJJJJJFJJJFJJJJJJJJJJJJJJJJJJJJJJJJFJJJJJJF<7<FJJ\tXA:Z:chr14,-19792846,150M,0;\tMC:Z:150M\tMD:Z:30A119\tRG:Z:NA12878."
        "1\tNM:i:1\tAS:i:145\tXS:i:150\n"
        "SRR14724454.285163302\t163\tchr22\t16050016\t27\t150M\t="
        "\t16050210\t344\tAGGACTTCAGAAGAGCTGTGAGACCTTGGCCAAGTCACTTCCTCCTTCAGGAACATTGCAGTGGGCCTAAGTGCCTCCTCTCGGGACTGGTATAGGGACGGTCATGCAATCTG"
        "GACAACATTAACCTTTAAAAGTTTATTGATCTTTTGT\tAAFF-"
        "FJJJJJJJJJJJJJJJJJFJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJFJJJJJJJJFJJAJJ-<AJ7FAJ7FJJJJJJFJJJJFJJFFJJJ-"
        "7FJJJJJFAJJFF<A<FJFJJAJJJJJ\tXA:Z:chr14,-19792836,150M,3;\tMC:Z:150M\tMD:Z:93G28C27\tRG:Z:NA12878.1\tNM:i:2\tAS:i:140\tXS:i:135\n"
        "SRR14724454.215616012\t163\tchr22\t16050042\t27\t150M\t="
        "\t16050404\t512\tTGGCCAAGTCACTTCCTCCTTCAGGAACATTGCAGTGGGCCTAAGTGCCTCCTCTCGGGACTGGTATGGGGACGGTCATGCAATCTGGACAACATTCACCTTTAAAAGTTTAT"
        "TGATCTTTTGTGACATGCACGTGGGTTCCCAGTAGCA\tA<"
        "AFFJAFJJFJJJJJJJJJJJJJFJJJJJFJJJJJAJJJJJJJJFJJFAJF7AFJJJJJJJJJJJJFJJFAJFJJJJFFJJJFJJJFJJJJJJJJJJJJJJJJJJFFAJAJJFFJJJJFJJJJJAAFJJ<"
        "FFJJJJJFAJA<FJJJJJJ\tXA:Z:chr14,-19792810,150M,0;\tMC:Z:150M\tMD:Z:150\tRG:Z:NA12878.1\tNM:i:0\tAS:i:150\tXS:i:150\n"
        "SRR14724454.80367511\t163\tchr22\t16050044\t27\t150M\t="
        "\t16050245\t351\tGCCAAGTCACTTCCTCCTTCAGGAACATTGCAGTGGGCCTAAGTGCCTCCTCTCGGGACTGGTATGGGGACGCTCATGCAATCTGGACAACATTAACCTTTAAAAGTTTATTG"
        "ATCTTTTGTGACATGCACGTGGGTTCCCAGTAGCAAG\tA<-AFF-AAFJFAFJ<7F<FJFFFAJJ-F7F-7-7FJ<<FF--A<J<JJ<AAJJ---A--AAJA7-7-FAJA7AA7-FF7F-<FJF-"
        "7F77FAF--AAAAJFAF-F--A7F<FJFJJ<<JA-7F-AFJ-AJ7<)A<)<)7)7A)<-7AFF\tXA:Z:chr14,-19792808,150M,2;\tMC:Z:150M\tMD:Z:72G21C55\tRG:Z:"
        "NA12878.1\tNM:i:2\tAS:i:140\tXS:i:140\n"
        "SRR14724454.187408578\t163\tchr22\t16050052\t27\t150M\t="
        "\t16050238\t336\tACTTCCTCCTTCAGGAACATTGCAGTGGGCCTAAGTGCCTCCTCTCGGGACTGGTATGGGGACGGTCATGCAATCTGGACAACATTCACCTTTAAAAGTTTATTGATCTTTTG"
        "TGACATGCACGTGGGTTCCCAGTAGCAAGAAACTAAA\tAAFFFJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJ"
        "JJJFJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJFJJJJJJJJJJJJJJJJ\tXA:Z:chr14,-19792800,150M,0;\tMC:Z:150M\tMD:Z:150\tRG:Z:NA12878."
        "1\tNM:i:0\tAS:i:150\tXS:i:150\n"
        "SRR14724454.100241018\t163\tchr22\t16050053\t40\t150M\t="
        "\t16050331\t428\tCTTCCTCCTTCAGGAACATTGCAGTGGGCCTAAGTGCCCCCTCTCGGGACTGGTATGGGGACGGTCATGCAATCTGGACAACATTCACCTTTAAAAGTTTATTGATCTTTTGT"
        "GACATGCACGTGGGTTCCCAGTAGCAAGAAACTAAAG\tAAFFFJFJJJJJJJJJJJJJJJJJJJJJJJFJJJJJJJ-7FJJJJJJJFJJJJFJJFJJJJJJFJFJJJJJJFJJJJJJJJJ<"
        "JJJJJJJJJJJJAAJJJJJJFJJFAJJJJJFFAJJFJFFAAFFAFFFJAAJJAJJJJJJJJJAJAFJ\tXA:Z:chr14,-19792799,150M,1;\tMC:Z:150M\tMD:Z:38T111\tRG:Z:"
        "NA12878.1\tNM:i:1\tAS:i:145\tXS:i:145\n"
        "SRR14724454.45348318\t99\tchr22\t16050054\t27\t150M\t="
        "\t16050233\t329\tTTCCTCCTTCAGGAACATTGCAGTGGGCCTAAGTGCCTCCTCTCGGGACTGGTATGGGGACGGTCATGCAATCTGGACAACATTCACCTTTAAAAGTTTATTGATTTTTTGTG"
        "ACATGCACGTGGGTTCCCAGTAGCAAGAAACTAAAGG\tAAFFFJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJAJJFJJJJJJJJJJJJJJFJJJJJJJJJJJJJJJJJJJJ"
        "AFJAFFFFJJJJ7AFFJJJJJJJJJFFJJFFFJJJJAAFJJJJ-<FFJA----7FA7AF\tXA:Z:chr14,-19792798,150M,1;\tMC:Z:150M\tMD:Z:105C44\tRG:Z:NA12878."
        "1\tNM:i:1\tAS:i:145\tXS:i:145\n"
        "SRR14724454.222980101\t163\tchr22\t16050058\t40\t150M\t="
        "\t16050325\t417\tTCCTTCAGGAACATTGCAGTGGGCCTAAGTGCCTCCTCTCGGGACTGGTATGGGGACGGTCATGCAATCTGGACAACATTCACCTTTAAAAGTTTATTGATTTTTTGTGACAT"
        "GCACGTGGGTTCCCAGTAGCAAGAAACTAAAGGGTCG\tAAFFFJJJJJJJJJJJJFJJJJJJFJJJJJJJJJJFJJJ<"
        "JJJJJJJJJJFJJJJJJJJAJFJJJFJJJJJJ7FJJJJJJAJJJJJFJJJJJFJFJJJJJJJJJJJJJJFJJJJJJFFJF<AAFJ7F<JJJJAJFFJJJJJJJJJJFAA<\tXA:Z:chr14,-"
        "19792794,150M,1;\tMC:Z:150M\tMD:Z:101C48\tRG:Z:NA12878.1\tNM:i:1\tAS:i:145\tXS:i:145\n"
        "SRR14724454.65269742\t99\tchr22\t16050068\t40\t150M\t="
        "\t16050320\t369\tACATTGCAGTGGGCCTAAGTGCCTCCTCTCGGGACTGGTATGGGGACGGTCATGCAATCTGGACAACATTCACCTTTAAAAGTTTATTGATCTTTTGTGACATGCACGTGGGT"
        "TCCCAGTAGCAAGAAACTAAAGGGTCGCAGGCCGGTT\t<AFFFJJJJJJJJJJJJJJJJJJJJJJFF<"
        "FJJJJJJJJJJJJJFAJJAJAJJJJJJFJJJJFJFJFJJJJFJJFFJJJJJJJJJJJJJJFJJJJJJJJAJ7F7FFFJJFAFJJAF7FJFAF--7A7A-A7<AFF-A7AF)77F-<F<JJ\tXA:Z:"
        "chr14,-19792784,150M,0;\tMC:Z:117M33S\tMD:Z:150\tRG:Z:NA12878.1\tNM:i:0\tAS:i:150\tXS:i:150\n"
        "SRR14724454.171029920\t99\tchr22\t16050076\t40\t150M\t="
        "\t16050351\t425\tGTGGGCCTAAGTGCCTCCTCTCGGGACTGGTATGGGGACGGTCATGCAATCTGGACAACATTCACCTTTAAAAGTTTATTGATCTTTTGTGACATGCACGTGGGTTCCCAGTA"
        "GCAAGAAACTAAAGGGTCGCAGGCCGGTTTCTGCTAA\tAAFFFJJJJJJJJJJJJJJJJJJJFJ<"
        "JJFAAJJJJJJJJJJJJJJJJJJFJJFFJJJAFJJJFFJFJJFJJJJFFJFJJJJJJJJJJFJJJJJJJJJJ<FJJJJJJJJ<<JJJJFAAJJJJJ7JJJJFF<JJ<JFAJJAJJJJJAFAFF\tXA:Z:"
        "chr14,-19792776,150M,0;\tMC:Z:150M\tMD:Z:150\tRG:Z:NA12878.1\tNM:i:0\tAS:i:150\tXS:i:150\n"
        "SRR14724454.122030801\t99\tchr22\t16050085\t21\t150M\t="
        "\t16050229\t294\tAGTGCCTCCTCTCGGGACTGGTATGGGGACGCTCATGCAATCTGGACAACATTCACCTTTAAAAGTTTATTGATCTTTTGTGACATGCACGTGGGTTCCCAGTAGCAAGAAAC"
        "TAAAGGGTCGCAGGCCGGTTTCTGCTAATTTCTTTAA\tAAFFFJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJ"
        "JJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJFJFJJJJJJJJJJJJJJJ\tXA:Z:chr14,-19792767,150M,1;\tMC:Z:150M\tMD:Z:31G118\tRG:Z:NA12878."
        "1\tNM:i:1\tAS:i:145\tXS:i:145\n"
        "SRR14724454.11297339\t163\tchr22\t16050098\t27\t150M\t="
        "\t16050240\t292\tGGGACTGGTATGGGGACGCTCATGCAATCTGGACAACATTCACCTTTAAAAGTTTATTGATCTTTTGTGACATGCACGTGGGTTCCCAGTAGCAAGAAACTAAAGGGTCGCAG"
        "GCCGGTTTCTGCTAATTTCTTTAATTCCAAGACAGTC\tAAFFFJJJJJJJJJJFFJJJJJJJJJJJJJJJJJFJJJJJJJJJJJJJJJJFFFJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJFJJJJFA"
        "JJAJJFJJJJJJJJFJFFJ-JJJJJJJJJJJJJJAF<FFFJJJ<JJJJJJJA<-7AA<J\tXA:Z:chr14,-19792754,150M,1;\tMC:Z:150M\tMD:Z:18G131\tRG:Z:NA12878."
        "1\tNM:i:1\tAS:i:145\tXS:i:145\n"
        "SRR14724454.329966935\t163\tchr22\t16050098\t40\t91M59S\t="
        "\t16050365\t417\tGGGACCGGTATGGAGACGATCATGCAATCTGGACAACATTCACCTTCAAAAGCTTATTGATTTTTTGTGACATGCACGTCGGTTCCCAGTACCAAGCGACTCCAGAACCGCAG"
        "CTCGACATAGCCTCGGTGCAGTTATTGACAGACAGTG\tAA-AF-<--A--7-AFJA--<F77FFF<F-<<FJFJJJF-<FJJ<--<FFJJ-<-F-7JJFF-7-FFF-<AJ-7FAF7---7<A7AJA--<"
        "-77-<---77-7<77--7-7-7--7--------7---7777)7-----7--7--77---\tXA:Z:chr14,-19792813,59S91M,7;\tMC:Z:150M\tMD:Z:"
        "5T7G4G27T5T8C17G11\tRG:Z:NA12878.1\tNM:i:7\tAS:i:56\tXS:i:56\n"
        "SRR14724454.94699838\t99\tchr22\t16050110\t27\t150M\t="
        "\t16050244\t284\tGGGACGCTCATGCAATCTGGACAACATTCACCTTTAAAAGTTTATTGATCTTTTGTGACATGCACGTGGGTTCCCAGTAGCAAGAAACTAAAGGGTCGCAGGCCGGTTTCTGC"
        "TAATTTCTTTAATTCCAAGACAGTCTCAAATATTTTC\tAAFFFJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJ"
        "JJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJ\tXA:Z:chr14,-19792742,150M,1;\tMC:Z:150M\tMD:Z:6G143\tRG:Z:NA12878."
        "1\tNM:i:1\tAS:i:145\tXS:i:145\n";
    std::vector<double> expected_values({-18.08620853341253, -4.235384193484483, -39.92524846633521});

    // 16050116
    samFile* in = sam_open(sam, "r");
    bam_hdr_t* h = sam_hdr_read(in);
    bam1_t* bam = bam_init1();
    for (int i = 0; i < 100; i++) _bit->set(i, true);

    HcActiveBase base(_pool, 2, 21, 16050115, 16050215, _bit, 0, ref, ref_len);

    HcActiveBase compute_base(_pool, 2, 21, 16050115, 16050215, _bit, 0, ref, ref_len);

    while (sam_read1(in, h, bam) >= 0) {
        base.process_bam_to_slot(bam, 0, 0);
        compute_base.process_bam_to_slot(bam, 0, 0);
    }
    base.compute_genotype_PL(0);
    auto likelihoods = base.get_current_likelihood();
    EXPECT_EQ(likelihoods.size(), expected_values.size());
    for (int i = 0; i < 3; i++) {
        EXPECT_NEAR(likelihoods[i], expected_values[i], 1E-8);
    }

    compute_base.compulte_all_likelihood();
    auto result = compute_base.get_active_result(0);

    EXPECT_NEAR(result.acitve_value, 0.9999999999930759, 1E-8);
    bam_destroy1(bam);
    bam_hdr_destroy(h);
    sam_close(in);
}

TEST_F(HCLikelihoodWGS, process_bam_to_slot_case2)
{
    const char sam[] =
        "data:,"
        "@SQ\tSN:chr1\tLN:249250621\n"
        "@SQ\tSN:chr2\tLN:243199373\n"
        "@SQ\tSN:chr3\tLN:198022430\n"
        "@SQ\tSN:chr4\tLN:191154276\n"
        "@SQ\tSN:chr5\tLN:180915260\n"
        "@SQ\tSN:chr6\tLN:171115067\n"
        "@SQ\tSN:chr7\tLN:159138663\n"
        "@SQ\tSN:chr8\tLN:146364022\n"
        "@SQ\tSN:chr9\tLN:141213431\n"
        "@SQ\tSN:chr10\tLN:135534747\n"
        "@SQ\tSN:chr11\tLN:135006516\n"
        "@SQ\tSN:chr12\tLN:133851895\n"
        "@SQ\tSN:chr13\tLN:115169878\n"
        "@SQ\tSN:chr14\tLN:107349540\n"
        "@SQ\tSN:chr15\tLN:102531392\n"
        "@SQ\tSN:chr16\tLN:90354753\n"
        "@SQ\tSN:chr17\tLN:81195210\n"
        "@SQ\tSN:chr18\tLN:78077248\n"
        "@SQ\tSN:chr19\tLN:59128983\n"
        "@SQ\tSN:chr20\tLN:63025520\n"
        "@SQ\tSN:chr21\tLN:48129895\n"
        "@SQ\tSN:chr22\tLN:51304566\n"
        "SRR14724454.101180596\t147\tchr22\t16060374\t60\t150M\t=\t16060196\t-"
        "328\tCTTTGCCTTCTGCCAGAATTGTGAACTTTCTGAGACCTCCCCAGAAATGGATGCCAGCATTATGCTTCCTATACAGCCTGCAGAACCATGAGCCAATTAACTCTCTTTTTTCTTTTTCTTTTTCT"
        "TTTTCTTTTTCTTTTTCTTTTTCTC\tJJJJJJJJJJJJF<"
        "JJJJJJJJJAFJAJJJJJJFFJJJJJJJJJJJJJJJJJJJJJJJFJJJJJJJJJJJJJJJJJJJJJJJFFJJJJJJJJFJJJJJJFJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJ"
        "JFFFAA\tXA:Z:chr14,+19782500,17S133M,2;\tMC:Z:150M\tMD:Z:149T0\tRG:Z:NA12878.1\tNM:i:1\tAS:i:149\tXS:i:123\n"
        "SRR14724454.126880449\t147\tchr22\t16060406\t60\t110M2D40M\t=\t16060212\t-"
        "346\tAGACCTCCCCAGAAATGGATGCCAGCATTATGCTTCCTATACAGCCTGCAGAACCATGAGCCAATTAACTCTCTTTTTTCTTTTTCTTTTTCTTTTTCTTTTTCTTTTTCTTTCTTTCTTTCTTT"
        "CTTTCTTTCTTTCTTTCTTTCTTTC\tJF<"
        "JJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJAJAJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJ"
        "JJJJJJJJJJJJFFFAA\tXA:Z:chr14,+19782451,150M,5;\tMC:Z:150M\tMD:Z:110^TT40\tRG:Z:NA12878.1\tNM:i:2\tAS:i:142\tXS:i:125\n"
        "SRR14724454.36687637\t83\tchr22\t16060414\t60\t108M2D42M\t=\t16060214\t-"
        "352\tCCAGAAATGGATGCCAGCATTATGCTTCCTATACAGCCTGCAGAACCATGAGCCAATTAACTCTCTTTTTTCTTTTTCTTTTTCTTTTTCTTTTTCTTTTTCTTTTTCTCTTTCTTTCTTTCTTT"
        "CTTTCTTTCTTTCTTTCTTTCTTTC\tJJJJJJJJJJJJJJJFJJJJJJJJF<JJJJJJJJAJJJJJJJJJFJF7JF<FFFJJJJF<"
        "JJJJFJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJFJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJFFFAA\tMC:Z:150M\tMD:Z:108^TT42\tRG:Z:"
        "NA12878.1\tNM:i:2\tAS:i:142\tXS:i:112\n"
        "SRR14724454.313097788\t147\tchr22\t16060428\t60\t52M6D98M\t=\t16060093\t-"
        "491\tCAGCATTATGCTTCCTATACTGCCTGCAGTACCATGAGCCAATTATCTCTCTTTTTTCTTTTTCTTTTTCTTTTTCTTTTTCTTTTTCTTTCTTTTTTTCTTTCTTTCTTTCTTTCTTTCTTTCT"
        "TTCTTTCTTTCTTTCTTTCTTTTTG\tJF77-<AF<AA-777--7)7)-AF<-7----<7-A--77A-<JF7-"
        "FFJJJJAFFJJJJFFJFAJJJJJJJJJJJJJFJJJJJFJJJJJJJJJJJJJJJJJJFJJJJJJJJJJJJJJJFFJJFJJJJJJJFJJJJJJJJJJJJJJFFFAA\tMC:Z:150M\tMD:Z:"
        "20A8A15A6^TTTTTC43C54\tRG:Z:NA12878.1\tNM:i:10\tAS:i:118\tXS:i:93\n"
        "SRR14724454.134730543\t147\tchr22\t16060465\t40\t55M2I46M1D14M33S\t=\t16060230\t-"
        "351\tGCCAATTAACTCTCTTTTTTCTTTTTCTTTTTCTTTTTCTTTTTCTTTTTCTTTTTCTCTTTCTTTCTTTCTTTCTTTCTTTCTTTCTTTCTTTCTTTCTTTCTTCTTTCTTTCTTTCTTTCTTT"
        "CTTTCTTTCTTTCTTTCTTTCTTTC\t7AFF<"
        "AJAAJFJFF7AJAJJJJFJFAJAJAAJJJJJJJJJFJJJJJJJJJJFJJJJFFJJFJFFJJJJJFFJJJFJJJJJJJJJFJFJJFFFAAFJJFJJJAJJFJJJJJJJJJJAFJFFFFJJFJJJFJJF7AF"
        "AAJJJAJFFAFF<AA\tMC:Z:150M\tMD:Z:101^T14\tRG:Z:NA12878.1\tNM:i:3\tAS:i:100\tXS:i:101\n"
        "SRR14724454.234428536\t99\tchr22\t16060471\t20\t9M6D97M44S\t="
        "\t16060625\t304\tTAACTCTCTTTTTTCTTTTTCTTTTTCTTTTTCTTTTTCTTTTTCTTTCTTTCTTTCTTTCTTTCTTTCTTTTTTTTTTTCTTTTTTTTTTTCTTTCTTTTTTTTTTTTTTTT"
        "TTTTTTTTTTTTTCTTCTTTTTTGTTTGTTTTTTGCT\tAAFFFJJJJJJJJJJJJJJJJJFJJJJJJJJJJJJJJJJJJJJJJJJJJJJJFJJJAJJF7-FFAJJJFJJJ-FFA-<JF-AJJ-"
        "7FJ7FJJ-FJJ-AJJ--FJJJ-AJFF-AJJ<7F77<F--7--7<-7AJFJJ-7AA<--7-7---7\tSA:Z:chr11,98877530,+,89S38M23S,0,0;\tXA:Z:chr9,+134295821,"
        "21S71M58S,2;chr19,+16386842,27S65M58S,1\tMC:Z:150M\tMD:Z:9^TTTTTC63C3C7C3C11C5\tRG:Z:NA12878.1\tNM:i:11\tAS:i:72\tXS:i:61\n"
        "SRR14724454.335292771\t83\tchr22\t16060472\t60\t8M6D142M\t=\t16060258\t-"
        "370\tAACTCTCTTTTTTCTTTTTCTTTTTCTTTTTCTTTTTCTTTTTCTTTCTTTCTTTCTTTCTTTCTTTCTTTCTTTCTTTCTTTCTTTCTTTCTTTCTTTCTTTTTGTTTTCTTTCATCTTTCCTT"
        "CTTCTTTTTTGATGGAGTCTCACTC\tAAJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJ"
        "JJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJFAFAA\tMC:Z:150M\tMD:Z:8^TTTTTC142\tRG:Z:NA12878.1\tNM:i:6\tAS:i:142\tXS:i:91\n"
        "SRR14724454.55609990\t83\tchr22\t16060486\t60\t8S142M\t=\t16060258\t-"
        "370\tTACTCTCTTTTTTCTTTTTCTTTTTCTTTTTCTTTTTCTTTTTCTTTCTTTCTTTCTTTCTTTCTTTCTTTCTTTCTTTCTTTCTTTCTTTCTTTCTTTCTTTTTGTTTTCTTTCATCTTTCCTT"
        "CTTCTTTTTTGATGGAGTCTCACTC\t-<JJFA<JJJJJFFJJJJJFJJJJJJJJJJF<JFJJF<JJJJJJJJJFFJJAJJJJJJJJJJJFJJJJAJA-"
        "FFFJJJFFJJJJJJJJJFFJJJFFJJJJJJJJJJJJJA-JJJA<JJJJJJFFJJJJJJJJJJJJFJJJJJFAJFA<AA\tMC:Z:150M\tMD:Z:142\tRG:Z:NA12878.1\tNM:i:0\tAS:i:"
        "142\tXS:i:91\n"
        "SRR14724454.142114690\t99\tchr22\t16060486\t60\t6S127M17S\t="
        "\t16060650\t314\tCTCTCTTTTTTCTTTTTCTTTTTCTTTTTCTTTTTCTTTTTCTTTCTTTTTTTCTTTCTTTCTTTCTTTCTTTCTTTCTTTTTTTCTTTCTTTCTTTCTTTTTGTTTTTTTTC"
        "ATCTTTCTTTCTTTTTTTTTTGTTTTTCTTTGCTTTG\tAAAFFJJAFAJ<<AJAFJAF-AAF<JA7<-FJ<7JJAFJJJ-JJJ7<A<-FFA---<AJFFAFFFJ-F<--FF<-<7-AJJ-FAAJ<JFF<"
        "7-<F<F-FFJJF-<<7-7FFJ--A-7-77----77-7J--A---7---7-7--------\tMC:Z:150M\tMD:Z:43C31C26C11C5C6\tRG:Z:NA12878.1\tNM:i:5\tAS:i:"
        "102\tXS:i:62\n"
        "SRR14724454.305261164\t147\tchr22\t16060492\t60\t150M\t=\t16060228\t-"
        "414\tTTTTTCTTTTTCTTTTTCTTTTTCTTTTTCTTTCTTTCTTTCTTTCTTTCTTTCTTTCTTTCTTTCTTTCTTTCTTTCTTTCTTTCTTTTTGTTTTCTTTCATCTTTCCTTCTTCTTTTTTGATG"
        "GAGTCTCACTCTGTTGCCTGGGCTG\tJJJJJJJJJJJJJJJJFJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJFJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJ"
        "JJJJJJJFFFJJJJJJJJJJJJJJJJJJJFJJJJJJJJJJJJFFFAA\tMC:Z:150M\tMD:Z:150\tRG:Z:NA12878.1\tNM:i:0\tAS:i:150\tXS:i:105\n"
        "SRR14724454.73810296\t147\tchr22\t16060519\t60\t150M\t=\t16060304\t-"
        "365\tTTCTTTCTTTCTTTCTTTCTTTCTTTCTTTCTTTCTTTCTTTCTTTCTTTCTTTCTTTCTTTTTGTTTTCTTTCATCTTTCCTTCTTCTTTTTTGATGGAGTCTCACTCTGTTGCCTGGGCTGGA"
        "GTGCAGTAAGTGGTACGATATTGGC\tJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJ"
        "JJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJFFFAA\tXA:Z:chr14,+19782303,128M1I13M2D8M,4;\tMC:Z:150M\tMD:Z:150\tRG:Z:NA12878.1\tNM:i:"
        "0\tAS:i:150\tXS:i:132\n"
        "SRR14724454.353109752\t147\tchr22\t16060523\t27\t33S117M\t=\t16060293\t-"
        "347\tCTTTCTTTCTTTCTTTCTTTCTTTCTTTCTTTCTTCTTTCTTTCTTTCTTTCTTTCTTTCTTTCTTTCTTTCTTTCTTTCTTTCTTTCTTTTTGTTTTCTTTCATCTTTCCTTCTTCTTTTTTGA"
        "TGGAGTCTCACTCTGTTGCCTGGGC\tJFAJJJJJJJAFFAJJJJFJFFJJJJJJJJFJJJJJJJJJJJJJJJJJJJJJJJJJFJJJJJJJJJJJJJJJJJJJJJFJJFJJJJFJJJJJJJJJJFJJJFJ"
        "JJJJJJJJJJJJJJFJJJJJJJJJJJJJJJJJJJJJJJJJJJFFFAA\tXA:Z:chr14,+19782332,99M1I15M1I34M,3;\tMC:Z:150M\tMD:Z:117\tRG:Z:NA12878.1\tNM:i:"
        "0\tAS:i:117\tXS:i:129\n";

    samFile* in = sam_open(sam, "r");
    bam_hdr_t* h = sam_hdr_read(in);
    bam1_t* bam = bam_init1();
    for (int i = 0; i < 100; i++) _bit->set(i, true);
    std::vector<double> expected_values({-11.731726414803024, -3.6130970503309916, -40.79483441063536});
    // memset(_buffer, 0, sizeof(uint8_t) * s_max_bufer_size);
    HcActiveBase base(_pool, 2, 21, 16060522, 16060623, _bit, 0, ref, ref_len);
    HcActiveBase compute_base(_pool, 2, 21, 16060522, 16060623, _bit, 0, ref, ref_len);
    int num = 0;
    while (sam_read1(in, h, bam) >= 0) {
        num++;
        base.process_bam_to_slot(bam, 0, 0);
        compute_base.process_bam_to_slot(bam, 0, 0);
    }
    EXPECT_EQ(num, 12);
    base.compute_genotype_PL(0);
    auto likelihoods = base.get_current_likelihood();
    EXPECT_EQ(likelihoods.size(), expected_values.size());
    for (int i = 0; i < 3; i++) {
        EXPECT_NEAR(likelihoods[i], expected_values[i], 1E-13);
    }

    compute_base.compulte_all_likelihood();
    auto result = compute_base.get_active_result(0);

    EXPECT_NEAR(result.acitve_value, 0.9999956312137954, 1E-12);
    EXPECT_EQ(result.extend_length, 33);
    EXPECT_EQ(result.state, BaseState::HC_REGION_BASES_ACTIVE_HQ_SOFT_CLIPS);
    bam_destroy1(bam);
    bam_hdr_destroy(h);
    sam_close(in);
}
