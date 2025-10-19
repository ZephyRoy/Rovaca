#include <gtest/gtest.h>

#include "bam_data_pool.hpp"
#include "cigar_builder.h"
#include "haplotype.h"
#include "read_clipper.h"
#include "read_record.h"
#include "smithwaterman_common.h"
#include "unit_test_utils.hpp"
#include "utils/alignment_utils.h"
#include "utils/base_utils.h"
#include "utils/cigar_utils.h"
#include "utils/debug_utils.h"

using namespace rovaca;

static constexpr size_t s_buffer_size = 1024 * 1024 * 1000;

struct create_read_aligned_to_ref_data
{
    pReadRecord read;
    pHaplotype hap;
    pHaplotype ref;
    uint32_t ref_start;
    uint32_t expected_start;
    pCigar expected_cigar;
} test_data;

class AlignmentUtilsUnitTest : public ::testing::Test
{
protected:
    uint8_t* _buffer{};
    pMemoryPool _pool{};
    pBamDataPool _bampool{};

    void SetUp() override
    {
        _buffer = new uint8_t[s_buffer_size]{};
        _pool = new std::pmr::monotonic_buffer_resource(_buffer, s_buffer_size, std::pmr::null_memory_resource());
        _bampool = new BamDataPool(s_buffer_size);
    }

    void TearDown() override
    {
        delete _pool;
        delete[] _buffer;
        delete _bampool;
    }

    static std::pmr::vector<create_read_aligned_to_ref_data> make_test_data(pMemoryPool pool, pBamDataPool bampool);
    static pReadRecord create_artificial_read(const char* seq, uint32_t len, pMemoryPool pool, pBamDataPool bampool);

    static std::vector<std::tuple<std::string, int64_t, int64_t>> testReadStartOnReferenceHaplotypeData();
    static std::vector<std::tuple<std::string, int64_t, int64_t, std::string>> testTrimCigarData();
    static std::vector<std::tuple<std::string, int64_t, int64_t, std::string>> testTrimCigarByBaseData();
    static std::vector<std::tuple<std::string, std::string, std::string>> testApplyCigarToCigarData();
    static std::vector<std::tuple<std::string, std::string, std::string, std::string>> testLeftAlignIndelData();
    static std::vector<std::tuple<std::string, std::string, std::string>> testAppendClippedElementsFromOriginalCigarData();
    static std::vector<std::tuple<std::string, int64_t, int64_t, std::string, std::string>> testGetBasesCoveringRefIntervalData();
    static std::vector<std::tuple<std::string, std::string, std::string, std::string>> testGetBasesAndBaseQualitiesAlignedOneToOneData();
};

TEST_F(AlignmentUtilsUnitTest, test_create_read_aligned_to_ref)
{
    std::pmr::vector<create_read_aligned_to_ref_data> all_tests = make_test_data(_pool, _bampool);
    p_lib_sw_avx sw = sw_avx_init(0);

    for (auto testcase : all_tests) {
        AlignmentUtils::create_read_aligned_to_ref(testcase.read, testcase.hap, testcase.ref, testcase.ref_start, true, sw, _pool,
                                                   _bampool);
        ASSERT_EQ(testcase.read->get_start(), testcase.expected_start);
        ASSERT_EQ(testcase.read->cigar_length(), testcase.expected_cigar->num);
        ASSERT_EQ(testcase.read->get_length_on_reference(), testcase.read->get_length_on_reference());
        for (uint32_t i = 0; i < testcase.read->cigar_length(); i++) {
            ASSERT_EQ(testcase.read->cigar_i(i), testcase.expected_cigar->data[i]);
        }
    }
    sw_avx_finit(sw);
}

TEST_F(AlignmentUtilsUnitTest, testReadStartOnReferenceHaplotype)
{
    auto data = testReadStartOnReferenceHaplotypeData();
    int64_t read_start_on_haplotype, expected_offset_in_ref;
    for (size_t i = 0, len = data.size(); i < len; ++i) {
        const auto& tup = data.at(i);
        const std::string& cp = std::get<0>(tup);
        read_start_on_haplotype = std::get<1>(tup);
        expected_offset_in_ref = std::get<2>(tup);
        pCigar cigar = DebugUtils::str2cigar(cp, _pool);
        ASSERT_EQ(expected_offset_in_ref, AlignmentUtils::read_start_on_reference_haplotype(cigar, read_start_on_haplotype)) << "i=" << i;
    }
}

TEST_F(AlignmentUtilsUnitTest, testTrimCigar)
{
    auto data = testTrimCigarData();
    for (size_t i = 0, len = data.size(); i < len; ++i) {
        const auto& tup = data.at(i);
        const std::string& cigar_str = std::get<0>(tup);
        int64_t start = std::get<1>(tup);
        int64_t length = std::get<2>(tup);
        const std::string& expected_cigar_str = std::get<3>(tup);

        pCigar cigar = DebugUtils::str2cigar(cigar_str, _pool);

        pCigar expected_cigar = DebugUtils::str2cigar(expected_cigar_str, _pool);
        if (expected_cigar->num == 1 && cigar_op_is_del(bam_cigar_op(expected_cigar->data[0]))) {
            // trimming throws error if all but deletion elements are trimmed
            continue;
        }

        pCigar result1 = CigarBuilder::create(_pool)->add_all(expected_cigar->data, expected_cigar->num)->make();
        pCigar result2 = AlignmentUtils::trim_cigar_by_reference(cigar, start, length, _pool);

        ASSERT_EQ(result1->num, result2->num) << "i=" << i;
        for (uint32_t j = 0; j < result1->num; ++j) {
            ASSERT_EQ(result1->data[j], result2->data[j]) << "i=" << i << " j=" << j;
        }
    }
}

TEST_F(AlignmentUtilsUnitTest, testTrimCigarByBase)
{
    auto data = testTrimCigarByBaseData();
    for (size_t i = 0, len = data.size(); i < len; ++i) {
        const auto& tup = data.at(i);
        const std::string& cigar_str = std::get<0>(tup);
        int64_t start = std::get<1>(tup);
        int64_t length = std::get<2>(tup);
        const std::string& expected_cigar_str = std::get<3>(tup);

        pCigar cigar = DebugUtils::str2cigar(cigar_str, _pool);

        pCigar result2 = DebugUtils::str2cigar(expected_cigar_str, _pool);
        pCigar result1 = AlignmentUtils::trim_cigar_by_bases(cigar, start, length, _pool);

        ASSERT_EQ(result1->num, result2->num) << "i=" << i;
        for (uint32_t j = 0; j < result1->num; ++j) {
            ASSERT_EQ(result1->data[j], result2->data[j]) << "i=" << i << " j=" << j;
        }
    }
}

TEST_F(AlignmentUtilsUnitTest, testApplyCigarToCigar)
{
    auto data = testApplyCigarToCigarData();
    for (size_t i = 0, len = data.size(); i < len; ++i) {
        const auto& tup = data.at(i);
        const std::string& cigar_str0 = std::get<0>(tup);
        const std::string& cigar_str1 = std::get<1>(tup);
        const std::string& cigar_str2 = std::get<2>(tup);

        pCigar cigar0 = DebugUtils::str2cigar(cigar_str0, _pool);
        pCigar cigar1 = DebugUtils::str2cigar(cigar_str1, _pool);

        pCigar result1 = DebugUtils::str2cigar(cigar_str2, _pool);
        pCigar result2 = AlignmentUtils::apply_cigar_to_cigar(cigar0, cigar1, _pool);
        ASSERT_EQ(result1->num, result2->num) << "i=" << i;
        for (uint32_t j = 0; j < result1->num; ++j) {
            ASSERT_EQ(result1->data[j], result2->data[j]) << "i=" << i << " j=" << j;
        }
    }
}

TEST_F(AlignmentUtilsUnitTest, testLeftAlignIndel)
{
    using namespace UnitTestUtils;
    std::mt19937 gen{1024};
    std::vector<int32_t> vec05{0, 5};
    std::vector<int32_t> vec010{0, 10};
    auto data = testLeftAlignIndelData();
    for (size_t i = 0, len = data.size(); i < len; ++i) {
        const auto& tup = data.at(i);
        const std::string& ref_string = std::get<0>(tup);
        const std::string& read_string = std::get<1>(tup);
        const std::string& original_cigar = std::get<2>(tup);
        const std::string& expected_cigar = std::get<3>(tup);
        for (int32_t leading_hard_clips : vec05) {
            for (int32_t trailing_hard_clips : vec05) {
                for (int32_t leading_soft_clips : vec05) {
                    for (int32_t trailing_soft_clips : vec05) {
                        for (int32_t extra_ref_in_front : vec010) {
                            for (int32_t extra_ref_in_back : vec010) {
                                uint32_t read_len = uint32_t(read_string.length() + leading_soft_clips + trailing_soft_clips);
                                uint32_t ref_len = uint32_t(ref_string.length() + extra_ref_in_front + extra_ref_in_back);
                                pBases read_bases = new ALLOC_FLEXIBLE_IN_POOL(_pool, Bases, read_len, uint8_t) Bases{read_len};
                                pBases ref_bases = new ALLOC_FLEXIBLE_IN_POOL(_pool, Bases, ref_len, uint8_t) Bases{ref_len};

                                BaseUtils::fill_with_random_bases(read_bases, 0, leading_soft_clips, gen);
                                BaseUtils::fill_with_random_bases(read_bases, int32_t(leading_soft_clips + read_string.length()),
                                                                  int32_t(read_bases->num), gen);
                                memcpy(read_bases->data + leading_soft_clips, read_string.c_str(), read_string.size() * sizeof(uint8_t));

                                BaseUtils::fill_with_random_bases(ref_bases, 0, extra_ref_in_front, gen);
                                BaseUtils::fill_with_random_bases(ref_bases, int32_t(extra_ref_in_front + ref_string.length()),
                                                                  int32_t(ref_bases->num), gen);
                                memcpy((char*)ref_bases->data + extra_ref_in_front, ref_string.c_str(),
                                       ref_string.size() * sizeof(uint8_t));

                                std::string original_cigar_with_clips =
                                    (leading_hard_clips > 0 ? std::to_string(leading_hard_clips) + "H" : "") +
                                    (leading_soft_clips > 0 ? std::to_string(leading_soft_clips) + "S" : "") + original_cigar +
                                    (trailing_soft_clips > 0 ? std::to_string(trailing_soft_clips) + "S" : "") +
                                    (trailing_hard_clips > 0 ? std::to_string(trailing_hard_clips) + "H" : "");
                                std::string expected_cigar_with_clips =
                                    (leading_hard_clips > 0 ? std::to_string(leading_hard_clips) + "H" : "") +
                                    (leading_soft_clips > 0 ? std::to_string(leading_soft_clips) + "S" : "") + expected_cigar +
                                    (trailing_soft_clips > 0 ? std::to_string(trailing_soft_clips) + "S" : "") +
                                    (trailing_hard_clips > 0 ? std::to_string(trailing_hard_clips) + "H" : "");

                                pCigar cigar1 = DebugUtils::str2cigar(original_cigar_with_clips, _pool);
                                pCigar result1 = DebugUtils::str2cigar(expected_cigar_with_clips, _pool);
                                pCigar result2 =
                                    AlignmentUtils::left_align_indels(cigar1, ref_bases, read_bases, extra_ref_in_front, _pool).cigar;
                                ASSERT_EQ(result1->num, result2->num) << "i=" << i;
                                for (uint32_t j = 0; j < result1->num; ++j) {
                                    ASSERT_EQ(result1->data[j], result2->data[j]) << "i=" << i << " j=" << j;
                                }
                            }
                        }
                    }
                }
            }
        }
    }
}

TEST_F(AlignmentUtilsUnitTest, testAppendClippedElementsFromOriginalCigar)
{
    auto data = testAppendClippedElementsFromOriginalCigarData();
    for (size_t i = 0; i < data.size(); ++i) {
        const auto& tup = data.at(i);
        const std::string& original_cigar_str = std::get<0>(tup);
        const std::string& shifted_cigar_str = std::get<1>(tup);
        const std::string& expected_str = std::get<2>(tup);

        pCigar original = DebugUtils::str2cigar(original_cigar_str, _pool);
        pCigar shifted = DebugUtils::str2cigar(shifted_cigar_str, _pool);

        pCigar result1 = DebugUtils::str2cigar(expected_str, _pool);
        pCigar result2 = AlignmentUtils::append_clipped_elements_from_cigar_to_cigar(shifted, original->data, original->num, _pool);
        ASSERT_EQ(result1->num, result2->num) << "i=" << i;
        for (uint32_t j = 0; j < result1->num; ++j) {
            ASSERT_EQ(result1->data[j], result2->data[j]) << "i=" << i << " j=" << j;
        }
    }
}

TEST_F(AlignmentUtilsUnitTest, testGetBasesCoveringRefInterval)
{
    auto data = testGetBasesCoveringRefIntervalData();
    for (size_t i = 0; i < data.size(); ++i) {
        const auto& tup = data.at(i);
        const std::string& bases_string = std::get<0>(tup);
        int64_t ref_start = std::get<1>(tup);
        int64_t ref_end = std::get<2>(tup);
        const std::string& cigar_str = std::get<3>(tup);
        const std::string& expected = std::get<4>(tup);

        pBases bases_string_b = UnitTestUtils::str2bases(bases_string, _pool);
        pCigar cigar_str_c = DebugUtils::str2cigar(cigar_str, _pool);

        pBases actual = AlignmentUtils::get_bases_covering_ref_interval(ref_start, ref_end, bases_string_b, 0, cigar_str_c, _pool);

        if (expected.empty()) {
            ASSERT_FALSE(actual) << "i=" << i;
        }
        else {
            ASSERT_EQ(expected.size(), actual->num) << "i=" << i;
            for (uint32_t j = 0; j < actual->num; ++j) {
                ASSERT_EQ(expected.at(j), actual->data[j]) << "i=" << i << " j=" << j;
            }
        }
    }
}

TEST_F(AlignmentUtilsUnitTest, testGetBasesAndBaseQualitiesAlignedOneToOne)
{
    auto data = testGetBasesAndBaseQualitiesAlignedOneToOneData();
    UnitTestUtils::pSamHeader header = UnitTestUtils::create_artificial_sam_header(10, 0, 1000);
    for (size_t i = 0; i < data.size(); ++i) {
        const auto& tup = data.at(i);
        const std::string& read_bases = std::get<0>(tup);
        const std::string& cigar = std::get<1>(tup);
        const std::string& expected_bases = std::get<2>(tup);
        const std::string& expected_quals = std::get<3>(tup);

        std::string quals_str(read_bases.size(), 10);

        pReadRecord read = UnitTestUtils::create_artificial_read(header->header, "default_read", 0, 1000, read_bases.size(), _bampool,
                                                                 _pool, cigar, read_bases, quals_str);

        std::pair<pBases, pBases> actual = AlignmentUtils::get_bases_and_base_qualities_aligned_one_to_one(read, _pool);
        pBases result_bases = actual.first;
        pBases result_quals = actual.second;

        ASSERT_EQ(result_bases->num, expected_bases.size()) << "i=" << i;
        for (uint32_t j = 0; j < result_bases->num; ++j) {
            ASSERT_EQ(result_bases->data[j], expected_bases.at(j)) << "i=" << i << " j=" << j;
        }

        ASSERT_EQ(result_quals->num, expected_quals.size()) << "i=" << i;
        for (uint32_t j = 0; j < result_quals->num; ++j) {
            ASSERT_EQ(result_quals->data[j], expected_quals.at(j)) << "i=" << i << " j=" << j;
        }
    }
    sam_header_destroy(header);
}

std::pmr::vector<create_read_aligned_to_ref_data> AlignmentUtilsUnitTest::make_test_data(pMemoryPool pool, pBamDataPool bampool)
{
    std::pmr::vector<create_read_aligned_to_ref_data> all_tests(pool);

    const char* hap = "ACTGAAGGTTCC";
    pBases hap_bases = new ALLOC_FLEXIBLE_IN_POOL(pool, Bases, 12 + 1, uint8_t) Bases{12 + 1};
    hap_bases->num = 12;
    for (int i = 0; i < 12; i++) {
        hap_bases->data[i] = hap[i];
    }

    pHaplotype all_m_ref = Haplotype::create(pool);
    pHaplotype all_m = Haplotype::create(pool);
    all_m_ref->init_haplotype(hap_bases, 1);
    all_m->init_haplotype(hap_bases, 0);

    all_m_ref->set_alignment_start_hap_wrt_ref(0);
    all_m->set_alignment_start_hap_wrt_ref(0);

    pCigarBuilder hap_cigar_builder = CigarBuilder::create(pool);
    hap_cigar_builder->add((12 << BAM_CIGAR_SHIFT | BAM_CMATCH));
    pCigar original_hap_cigar = hap_cigar_builder->make();
    all_m->set_cigar(original_hap_cigar);
    all_m_ref->set_cigar(original_hap_cigar);

    // make sure we get back a cigar of the right length
    pReadRecord original_read = create_artificial_read(hap, hap_bases->num, pool, bampool);
    create_read_aligned_to_ref_data original_testcase{original_read, all_m, all_m_ref, 10, 10, all_m->cigar()};
    all_tests.emplace_back(original_testcase);

    const char* bases = "ACTAAAGGTTCC";
    pReadRecord read = create_artificial_read(bases, hap_bases->num, pool, bampool);
    create_read_aligned_to_ref_data testcase{read, all_m, all_m_ref, 10, 10, all_m->cigar()};
    all_tests.emplace_back(testcase);

    // make sure insertions at the front are correctly handled
    const char* front_padded_bases = "NNNNNNNNACTGAAGGTTCC";
    pReadRecord read1 = create_artificial_read(front_padded_bases, 20, pool, bampool);
    pCigarBuilder pad_front_cigar_builder = CigarBuilder::create(pool);
    pad_front_cigar_builder->add((8 << BAM_CIGAR_SHIFT | BAM_CINS))->add((12 << BAM_CIGAR_SHIFT | BAM_CMATCH));
    create_read_aligned_to_ref_data testcase1{read1, all_m, all_m_ref, 10, 10, pad_front_cigar_builder->make()};
    all_tests.emplace_back(testcase1);

    // make sure insertions at the end are correctly handled
    const char* back_padded_bases = "ACTGAAGGTTCCNNNNNNNN";
    pReadRecord read2 = create_artificial_read(back_padded_bases, 20, pool, bampool);
    hap_cigar_builder->add((8 << BAM_CIGAR_SHIFT | BAM_CINS));
    create_read_aligned_to_ref_data testcase2{read2, all_m, all_m_ref, 10, 10, hap_cigar_builder->make()};
    all_tests.emplace_back(testcase2);

    // make sure refStart and hapStart are respected
    uint32_t ref_start = 5;
    for (uint32_t hap_start = ref_start; hap_start < 10 + ref_start; hap_start++) {
        pHaplotype haplo = Haplotype::create(pool);
        haplo->init_haplotype(hap_bases, 0);
        haplo->set_alignment_start_hap_wrt_ref(hap_start);
        haplo->set_cigar(original_hap_cigar);
        pReadRecord read = create_artificial_read(hap, hap_bases->num, pool, bampool);
        create_read_aligned_to_ref_data testcase{read, haplo, all_m_ref, ref_start, ref_start + hap_start, all_m->cigar()};
        all_tests.emplace_back(testcase);
    }

    // example case of bad alignment because SW doesn't necessarily left-align indels
    {
        const char* hap_bases =
            "ACTGTGGGTTCCTCTTATTTTATTTCTACATCAATGTTCATATTTAACTTATTATTTTATCTTATTTTTAAATTTCTTTTATGTTGAGCCTTGATGAAAGCCATAGGTTCTCTCATATAATTGTAT"
            "GTGT"
            "ATGTATGTATATGTACATAATATATACATATATGTATATGTATGTGTATGTACATAATATATACGTATATGTATGTGTATGTACATAATATATACGTATATGTATGTGTATGTACATAATATATAC"
            "GTAT"
            "ATGTATGTGTATGTACATAATATATACGTATATGTATGTGTATGTACATAATATATACGTATATGTATGTGTATGTGTATTACATAATATATACATATATGTATATATTATGTATATGTACATAAT"
            "ATAT"
            "ACATATATG";
        uint32_t hap_length = strlen(hap_bases);

        const char* bad_sw_read_bases =
            "ATGTACATAATATATACATATATGTATATGTATGTACATAATATATACGTATATGTATGTGTATGTACATAATATATACGTATATGTATGTGTATGTACATAATATATACGTATATGTATGTGTAT"
            "GTAC"
            "ATAATATATACGTATATGTATGTGTATGTACATAATATATACGTATATGTATGTGTATGTGTATTACATAATATATACATATATGTATATATTATGTATATGTACATAATAT";
        uint32_t read_length = strlen(bad_sw_read_bases);

        uint32_t ref_start = 10130100;
        uint32_t hap_start = 500;
        uint32_t expected_pos = 10130740;

        pBases bad_sw_hap_bases = new ALLOC_FLEXIBLE_IN_POOL(pool, Bases, hap_length + 1, uint8_t) Bases{hap_length + 1};
        bad_sw_hap_bases->num = hap_length;
        for (uint32_t i = 0; i < hap_length; i++) {
            bad_sw_hap_bases->data[i] = hap_bases[i];
        }
        pHaplotype bad_sw_hap = Haplotype::create(pool);
        bad_sw_hap->init_haplotype(bad_sw_hap_bases, 0);

        pCigarBuilder bad_sw_cigar_builder = CigarBuilder::create(pool);
        bad_sw_cigar_builder->add(hap_length << BAM_CIGAR_SHIFT | BAM_CMATCH);
        pCigar bad_sw_cigar = bad_sw_cigar_builder->make();
        bad_sw_hap->set_cigar(bad_sw_cigar);
        bad_sw_hap->set_alignment_start_hap_wrt_ref(hap_start);

        pReadRecord bad_sw_read = create_artificial_read(bad_sw_read_bases, read_length, pool, bampool);

        pCigarBuilder good_sw_cigar_builder = CigarBuilder::create(pool);
        good_sw_cigar_builder->add(28 << BAM_CIGAR_SHIFT | BAM_CMATCH)
            ->add(6 << BAM_CIGAR_SHIFT | BAM_CDEL)
            ->add(214 << BAM_CIGAR_SHIFT | BAM_CMATCH);
        pCigar good_sw_cigar = good_sw_cigar_builder->make();

        create_read_aligned_to_ref_data testcase{bad_sw_read, bad_sw_hap, bad_sw_hap, ref_start, expected_pos, good_sw_cigar};
        all_tests.emplace_back(testcase);
    }

    // example where left-align generates a leading deletion
    {
        const char* ref =
            "CTGAACGTAACCAAAATCAATATGGATACTGAGAAATACTATTTAATAAAGACATAAATTAGACTGCTAAAAAAAATTAAAGAAATTTCAAAAGAGAATCCACCTCTTTTCCTTGCCAGTGCTCAA"
            "AAGTGAGTGTGAATCTGGTGGCTGTGGGGCTGTTTTTGGTGTGGCTCTTTGGACCAGCCTGCCTGGTAATTCAAGCCTGCCTCTCATTTCTG";
        uint32_t ref_length = strlen(ref);

        const char* hap =
            "CTGAACGTAACCAAAATCAATATGGATACTGAGAAATACTATTTAATAAAGACATAAATTAGACTGCTAAAAAAAATTAAAGAAATTTCAAAAGAGAATCCACCTCTTTTCCTTGCCAGTGCTCAA"
            "AAGTGAGTGTGAATCTGGTGGCTGCGGGGCTGTTTTTGGTGTGGCTCTTTGGACCAGCCTGCCTGGTAATTCAAGCCTGCCTCTCATTTCTG";
        uint32_t hap_length = strlen(hap);

        const char* read = "GCTGCTTTTGGTGTGGCTCTTT";

        uint32_t ref_start = 215239171;
        uint32_t hap_start = 575;
        uint32_t alignment_offset = 154;
        uint32_t expected_pos = ref_start + hap_start + alignment_offset;

        pBases bad_hap_bases = new ALLOC_FLEXIBLE_IN_POOL(pool, Bases, hap_length + 1, uint8_t) Bases{hap_length + 1};
        bad_hap_bases->num = hap_length;
        for (uint32_t i = 0; i < hap_length; i++) {
            bad_hap_bases->data[i] = hap[i];
        }
        pHaplotype bad_hap = Haplotype::create(pool);
        bad_hap->init_haplotype(bad_hap_bases, 0);

        pBases ref_hap_bases = new ALLOC_FLEXIBLE_IN_POOL(pool, Bases, hap_length + 1, uint8_t) Bases{hap_length + 1};
        ref_hap_bases->num = ref_length;
        for (uint32_t i = 0; i < ref_length; i++) {
            ref_hap_bases->data[i] = ref[i];
        }
        pHaplotype ref_hap = Haplotype::create(pool);
        ref_hap->init_haplotype(ref_hap_bases, 1);

        pCigarBuilder hap_cigar_builder = CigarBuilder::create(pool);
        hap_cigar_builder->add(hap_length << BAM_CIGAR_SHIFT | BAM_CMATCH);
        pCigar hap_cigar = hap_cigar_builder->make();
        bad_hap->set_cigar(hap_cigar);
        bad_hap->set_alignment_start_hap_wrt_ref(hap_start);

        ref_hap->set_cigar(hap_cigar);
        ref_hap->set_alignment_start_hap_wrt_ref(0);

        pReadRecord read_record = create_artificial_read(read, strlen(read), pool, bampool);

        pCigarBuilder good_cigar_builder = CigarBuilder::create(pool);
        good_cigar_builder->add(22 << BAM_CIGAR_SHIFT | BAM_CMATCH);
        pCigar good_cigar = good_cigar_builder->make();

        create_read_aligned_to_ref_data testcase{read_record, bad_hap, ref_hap, ref_start, expected_pos, good_cigar};
        all_tests.emplace_back(testcase);
    }

    // example where the haplotype has an indel relative to reference
    {
        const char* ref =
            "GGGATCCTGCTACAAAGGTGAAACCCAGGAGAGTGTGGAGTCCAGAGTGTTGCCAGGACCCAGGCACAGGCATTAGTGCCCGTTGGAGAAAACAGGGGAATCCCGAAGAAATGGTGGGTCCTGGCC"
            "ATCCGTGAGATCTTCCCAGGGCAGCTCCCCTCTGTGGAATCCAATCTGTCTTCCATCCTGC";
        uint32_t ref_length = strlen(ref);

        const char* hap =
            "GGGATCCTGCTACAAAGGTGAAACCCAGGAGAGTGTGGAGTCCAGAGTGTTGCCAGGACCCAGGCACAGGCATTAGTGCCCGTTGGAGAAAACGGGAATCCCGAAGAAATGGTGGGTCCTGGCCAT"
            "CCGTGAGATCTTCCCAGGGCAGCTCCCCTCTGTGGAATCCAATCTGTCTTCCATCCTGC";
        uint32_t hap_length = strlen(hap);

        pBases hap_bases = new ALLOC_FLEXIBLE_IN_POOL(pool, Bases, hap_length + 1, uint8_t) Bases{hap_length + 1};
        hap_bases->num = hap_length;
        for (uint32_t i = 0; i < hap_bases->num; i++) {
            hap_bases->data[i] = hap[i];
        }
        pHaplotype phap = Haplotype::create(pool);
        phap->init_haplotype(hap_bases, 0);

        pBases ref_bases = new ALLOC_FLEXIBLE_IN_POOL(pool, Bases, ref_length + 1, uint8_t) Bases{ref_length + 1};
        ref_bases->num = ref_length;
        for (uint32_t i = 0; i < ref_bases->num; i++) {
            ref_bases->data[i] = ref[i];
        }
        pHaplotype ref_hap = Haplotype::create(pool);
        ref_hap->init_haplotype(ref_bases, 1);

        const char* read = "CCCATCCGTGAGATCTTCCCAGGGCAGCTCCCCTCTGTGGAATCCAATCTGTCTTCCATCCTGC";
        pReadRecord read_record = create_artificial_read(read, strlen(read), pool, bampool);

        uint32_t ref_start = 13011;
        uint32_t hap_start = 553;
        uint32_t alignment_offset = 123;
        uint32_t expected_pos = ref_start + hap_start + alignment_offset;

        pCigarBuilder hap_cigar_builder = CigarBuilder::create(pool);
        hap_cigar_builder->add(93 << BAM_CIGAR_SHIFT | BAM_CMATCH)
            ->add(2 << BAM_CIGAR_SHIFT | BAM_CDEL)
            ->add(92 << BAM_CIGAR_SHIFT | BAM_CMATCH);
        pCigar hap_cigar = hap_cigar_builder->make();
        phap->set_cigar(hap_cigar);
        phap->set_alignment_start_hap_wrt_ref(hap_start);

        pCigarBuilder ref_cigar_builder = CigarBuilder::create(pool);
        ref_cigar_builder->add(strlen(ref) << BAM_CIGAR_SHIFT | BAM_CMATCH);
        ref_hap->set_cigar(ref_cigar_builder->make());
        ref_hap->set_alignment_start_hap_wrt_ref(0);

        pCigarBuilder good_cigar_builder = CigarBuilder::create(pool);
        good_cigar_builder->add(64 << BAM_CIGAR_SHIFT | BAM_CMATCH);
        pCigar good_cigar = good_cigar_builder->make();

        create_read_aligned_to_ref_data testcase{read_record, phap, ref_hap, ref_start, expected_pos, good_cigar};
        all_tests.emplace_back(testcase);
    }
    return all_tests;
}

pReadRecord AlignmentUtilsUnitTest::create_artificial_read(const char* seq, uint32_t len, pMemoryPool pool, pBamDataPool bampool)
{
    bam1_t* new_bam = bampool->alloc_bam_struct_pool();

    const char* qname = "test_read";

    const char* qual = seq;

    pCigarBuilder cigar_builder = CigarBuilder::create(pool);
    pCigar cigar = cigar_builder->add((len << BAM_CIGAR_SHIFT | BAM_CMATCH))->make();

    // bam_set1(new_bam, qname_len, qname, flag, tid, pos, mapping_quality, cigar->num, cigar->data, mtid, mpos, isize, seq_len,
    // (char*)seq, (char *)qual, 0);
    bam_set1(new_bam, 10, qname, 0, 0, 1, 60, cigar->num, cigar->data, 0, 1, 0, len, seq, qual, 0);

    bampool->peek_bam(new_bam);
    return ReadRecord::create(pool, nullptr, new_bam);
}

std::vector<std::tuple<std::string, int64_t, int64_t>> AlignmentUtilsUnitTest::testReadStartOnReferenceHaplotypeData()
{
    return {{{"30M5D20M"}, 50, 55},      {{"30M5I20M"}, 50, 45},      {{"55M"}, 50, 50},
            {{"30M5D30M5D30M"}, 80, 90}, {{"30M5D30M5I30M"}, 80, 80}, {{"30M5D30M5I30M"}, 80, 80}};
}

std::vector<std::tuple<std::string, int64_t, int64_t, std::string>> AlignmentUtilsUnitTest::testTrimCigarData()
{
    using namespace UnitTestUtils;
    std::vector<std::tuple<std::string, int64_t, int64_t, std::string>> result;
    std::vector<int32_t> cigar_ops{BAM_CDEL, BAM_CEQUAL, BAM_CDIFF, BAM_CMATCH};
    std::vector<int32_t> pad_ops{BAM_CDEL, BAM_CMATCH};
    for (int32_t op : cigar_ops) {
        for (int32_t my_length = 1; my_length < 6; my_length++) {
            for (int32_t start = 0; start < my_length - 1; start++) {
                for (int32_t end = start; end < my_length; end++) {
                    int32_t length = end - start + 1;
                    for (int32_t padop : pad_ops) {
                        for (int32_t left_pad = 0; left_pad < 2; left_pad++) {
                            for (int32_t right_pad = 0; right_pad < 2; right_pad++) {
                                int32_t src_num = 0, expe_num = 0;
                                std::string left_str, right_str;
                                if (left_pad > 0) {
                                    left_str.append(std::to_string(left_pad)).append(1, BAM_CIGAR_STR[padop]);
                                }
                                left_str.append(std::to_string(my_length)).append(1, BAM_CIGAR_STR[op]);
                                if (right_pad > 0) {
                                    left_str.append(std::to_string(right_pad)).append(1, BAM_CIGAR_STR[padop]);
                                }
                                src_num = start + left_pad;
                                expe_num = end + left_pad;
                                right_str.append(std::to_string(length)).append(1, BAM_CIGAR_STR[op]);
                                result.push_back({left_str, src_num, expe_num, right_str});
                            }
                        }
                    }
                }
            }
        }
    }

    std::vector<int32_t> left_pads{0, 1, 2, 5};
    std::vector<int32_t> right_pads{0, 1, 2, 5};
    std::vector<int32_t> ins_sizes{1, 10};
    for (int32_t left_pad : left_pads) {
        for (int32_t right_pad : right_pads) {
            int32_t length = left_pad + right_pad;
            if (length > 0) {
                for (int32_t ins_size : ins_sizes) {
                    for (int32_t start = 0; start <= left_pad; start++) {
                        for (int32_t stop = left_pad; stop < length; stop++) {
                            int32_t left_pad_remaining = left_pad - start;
                            int32_t right_pad_remaining = stop - left_pad + 1;
                            std::string left_str, right_str;
                            left_str.append(std::to_string(left_pad))
                                .append(1, 'M')
                                .append(std::to_string(ins_size))
                                .append(1, 'I')
                                .append(std::to_string(right_pad))
                                .append(1, 'M');

                            if (left_pad_remaining > 0) {
                                right_str.append(std::to_string(left_pad_remaining)).append(1, 'M');
                            }
                            right_str.append(std::to_string(ins_size)).append(1, 'I');
                            if (right_pad_remaining > 0) {
                                right_str.append(std::to_string(right_pad_remaining)).append(1, 'M');
                            }
                            result.emplace_back(left_str, start, stop, right_str);
                        }
                    }
                }
            }
        }
    }

    result.emplace_back("3M2D4M", 0, 8, "3M2D4M");
    result.emplace_back("3M2D4M", 2, 8, "1M2D4M");
    result.emplace_back("3M2D4M", 2, 6, "1M2D2M");
    result.emplace_back("3M2D4M", 3, 6, "2D2M");
    result.emplace_back("3M2D4M", 4, 6, "1D2M");
    result.emplace_back("3M2D4M", 5, 6, "2M");
    result.emplace_back("3M2D4M", 6, 6, "1M");
    result.emplace_back("2M3I4M", 0, 5, "2M3I4M");
    result.emplace_back("2M3I4M", 1, 5, "1M3I4M");
    result.emplace_back("2M3I4M", 1, 4, "1M3I3M");
    result.emplace_back("2M3I4M", 2, 4, "3I3M");
    result.emplace_back("2M3I4M", 2, 3, "3I2M");
    result.emplace_back("2M3I4M", 2, 2, "3I1M");
    result.emplace_back("2M3I4M", 3, 4, "2M");
    result.emplace_back("2M3I4M", 3, 3, "1M");
    result.emplace_back("2M3I4M", 4, 4, "1M");

    return result;
}

std::vector<std::tuple<std::string, int64_t, int64_t, std::string>> AlignmentUtilsUnitTest::testTrimCigarByBaseData()
{
    return {{"2M3I4M", 0, 8, "2M3I4M"}, {"2M3I4M", 1, 8, "1M3I4M"}, {"2M3I4M", 2, 8, "3I4M"}, {"2M3I4M", 3, 8, "2I4M"},
            {"2M3I4M", 4, 8, "1I4M"},   {"2M3I4M", 4, 7, "1I3M"},   {"2M3I4M", 4, 6, "1I2M"}, {"2M3I4M", 4, 5, "1I1M"},
            {"2M3I4M", 4, 4, "1I"},     {"2M3I4M", 5, 5, "1M"},     {"2M2D2I", 0, 3, "2M2I"}, {"2M2D2I", 1, 3, "1M2I"},
            {"2M2D2I", 2, 3, "2I"},     {"2M2D2I", 3, 3, "1I"},     {"2M2D2I", 2, 2, "1I"},   {"2M2D2I", 1, 2, "1M1I"},
            {"2M2D2I", 0, 1, "2M"},     {"2M2D2I", 1, 1, "1M"}};
}

std::vector<std::tuple<std::string, std::string, std::string>> AlignmentUtilsUnitTest::testApplyCigarToCigarData()
{
    return {{"1M", "1M", "1M"},
            {"2M", "2M", "2M"},
            {"3M", "3M", "3M"},
            {"4M", "4M", "4M"},
            {"3M", "2M3D1M", "2M3D1M"},
            {"3M1I2M", "2M1D3M", "2M1D1M1I2M"},
            {"1M1D2M", "1M1D3M", "1M2D2M"},
            {"1M2D2M", "1M1D1M1I2M", "1M2D2M"},
            {"1M1I4M", "5M", "1M1I4M"},
            {"1M2D2M", "5M", "1M2D2M"},
            {"108M14D24M2M18I29M92M1000M", "2M1I3M", "2M1I3M"}};
}

std::vector<std::tuple<std::string, std::string, std::string, std::string>> AlignmentUtilsUnitTest::testLeftAlignIndelData()
{
    return {{"ACGT", "ACGT", "4M", "4M"},
            {"ACCT", "ACGT", "4M", "4M"},
            {"ACGT", "ACAT", "2M1X1M", "2M1X1M"},
            {"AAATTT", "AAACCCTTT", "3M3I3M", "3M3I3M"},
            {"CCCTTT", "AAACCCTTT", "3I6M", "3I6M"},
            {"AAACCC", "AAACCCTTT", "6M3I", "6M3I"},
            {"AAACCC", "AAACCGTTT", "6M3I", "6M3I"},
            {"AAACCCTTT", "AAATTT", "3M3D3M", "3M3D3M"},
            {"AAACCCTTT", "AAACCCCCCTTT", "5M3I4M", "3M3I6M"},
            {"AAACCCTTT", "AAACCCCCCTTT", "6M3I3M", "3M3I6M"},
            {"AAACCCTTT", "AAGCCCCCCTGT", "6M3I3M", "3M3I6M"},
            {"AAACGCGCGCGTTT", "AAACGCGCGCGCGCGTTT", "7M4I7M", "3M4I11M"},
            {"CCGCCG", "CCGCCGCCG", "6M3I", "3I6M"},
            {"ACCGCCG", "TCCGCCGCCG", "7M3I", "1M3I6M"},
            {"AAACCCCCCTTT", "AAACCCTTT", "5M3D4M", "3M3D6M"},
            {"AAACCCCCCTTT", "AAACCCTTT", "6M3D3M", "3M3D6M"},
            {"AAACGCGCGCGCGCGTTT", "AAACGCGCGCGTTT", "7M4D7M", "3M4D11M"},
            {"AAACCCTTTGGGAAA", "AAACCCCCCTTTGGGGGGAAA", "6M3I6M3I3M", "3M3I6M3I6M"},
            {"AAACCCTTTGGGGGGAAA", "AAACCCCCCTTTGGGAAA", "6M3I6M3D3M", "3M3I6M3D6M"},
            {"AAACCCCCTTT", "AAACCCCCTTT", "4M3I3D4M", "11M"},
            {"AAACCCCCTTT", "AAACCCCCTTT", "4M3D3I4M", "11M"},
            {"AAACCCCCTTT", "AAACCCCCTTT", "3M3I2M3D3M", "11M"},
            {"AACGCGCGCGTT", "AACGCGCGCGCGCGTT", "2M2I8M2I2M", "2M4I10M"},
            {"AACGCGCGCGCGCGTT", "AACGCGCGCGTT", "2M2D8M2D2M", "2M4D10M"}};
}

std::vector<std::tuple<std::string, std::string, std::string>> AlignmentUtilsUnitTest::testAppendClippedElementsFromOriginalCigarData()
{
    return {{"30M", "30M", "30M"},     {"30M", "15M6I15M", "15M6I15M"},     {"5S30M", "30M", "5S30M"},
            {"5H30M", "30M", "5H30M"}, {"5H5S30M", "30M", "5H5S30M"},       {"30M5H", "30M", "30M5H"},
            {"30M5S", "30M", "30M5S"}, {"10H30M5S5H", "30M", "10H30M5S5H"}, {"10H10M6D6M6D50M5S5H", "10M6I50M", "10H10M6I50M5S5H"}};
}

std::vector<std::tuple<std::string, int64_t, int64_t, std::string, std::string>>
AlignmentUtilsUnitTest::testGetBasesCoveringRefIntervalData()
{
    return {{"ACGT", 0, 3, "4M", "ACGT"},        {"ACGT", 1, 3, "4M", "CGT"},
            {"ACGT", 1, 2, "4M", "CG"},          {"ACGT", 1, 1, "4M", "C"},
            {"ACGT", 0, 5, "2M2D2M", "ACGT"},    {"ACGT", 1, 5, "2M2D2M", "CGT"},
            {"ACGT", 2, 5, "2M2D2M", ""},        {"ACGT", 3, 5, "2M2D2M", ""},
            {"ACGT", 4, 5, "2M2D2M", "GT"},      {"ACGT", 5, 5, "2M2D2M", "T"},
            {"ACGT", 0, 4, "2M2D2M", "ACG"},     {"ACGT", 0, 3, "2M2D2M", ""},
            {"ACGT", 0, 2, "2M2D2M", ""},        {"ACGT", 0, 1, "2M2D2M", "AC"},
            {"ACGT", 0, 0, "2M2D2M", "A"},       {"ACTTGT", 0, 3, "2M2I2M", "ACTTGT"},
            {"ACTTGT", 1, 3, "2M2I2M", "CTTGT"}, {"ACTTGT", 2, 3, "2M2I2M", "GT"},
            {"ACTTGT", 3, 3, "2M2I2M", "T"},     {"ACTTGT", 0, 2, "2M2I2M", "ACTTG"},
            {"ACTTGT", 0, 1, "2M2I2M", "AC"},    {"ACTTGT", 1, 2, "2M2I2M", "CTTG"},
            {"ACTTGT", 2, 2, "2M2I2M", "G"},     {"ACTTGT", 1, 1, "2M2I2M", "C"},
            {"ACTTGT", 0, 3, "2I4M", "TTGT"},    {"ACTTGT", 0, 3, "4M2I", "ACTT"},
            {"ACGT", 0, 1, "2M2I", "AC"},        {"ACGT", 1, 1, "2M2I", "C"},
            {"ACGT", 0, 0, "2M2I", "A"}};
}

std::vector<std::tuple<std::string, std::string, std::string, std::string>>
AlignmentUtilsUnitTest::testGetBasesAndBaseQualitiesAlignedOneToOneData()
{
    return {{"ATCGATCG", "8M", "ATCGATCG", {10, 10, 10, 10, 10, 10, 10, 10}},
            {"ATCGATCG", "4M4D4M", "ATCG----ATCG", {10, 10, 10, 10, 0, 0, 0, 0, 10, 10, 10, 10}},
            {"ATCGATCG", "2M3I3M", "ATTCG", {10, 10, 10, 10, 10}},
            {"ATCGATCG", "2I6M", "CGATCG", {10, 10, 10, 10, 10, 10}},
            {"ATCGATCG", "6M2I", "ATCGAT", {10, 10, 10, 10, 10, 10}},
            {"ATCGATCG", "1D8M1D", "-ATCGATCG-", {0, 10, 10, 10, 10, 10, 10, 10, 10, 0}},
            {"ATCGATCG", "2M1I1D2M2I1M2D", "AT-GAG--", {10, 10, 0, 10, 10, 10, 0, 0}}};
}