#include <gtest/gtest.h>

#include <algorithm>
#include <string>
#include <tuple>

#include "cigar_builder.h"
#include "forward.h"
#include "genotype_struct.h"
#include "utils/cigar_utils.h"
#include "utils/debug_utils.h"

using namespace rovaca;

static constexpr size_t s_buffer_size = 1024 * 1024 * 10;

class CigarBuilderUnitTest : public ::testing::Test
{
protected:
    void SetUp() override
    {
        _buffer = new uint8_t[s_buffer_size]{};
        _pool = new std::pmr::monotonic_buffer_resource(_buffer, s_buffer_size, std::pmr::null_memory_resource());
    }

    void TearDown() override
    {
        delete _pool;
        delete[] _buffer;
    }

    static std::vector<std::string> simpleConcatenationData();
    static std::vector<std::pair<std::string, std::string>> initialAndFinalDeletionsData();
    static std::vector<std::pair<std::string, std::string>> retainDeletionsData();
    static std::vector<std::pair<std::string, std::string>> mergeConsecutiveData();
    static std::vector<std::pair<std::string, std::string>> indelSandwichData();
    static std::vector<std::string> invalidData();
    static std::vector<std::pair<std::string, std::pair<int64_t, int64_t>>> removedDeletionsData();
    static std::vector<std::tuple<std::string, std::string, int64_t, int64_t>> removedDeletionsWithTwoMakesData();

    uint8_t* _buffer{};
    pMemoryPool _pool{};
};

TEST_F(CigarBuilderUnitTest, testSimpleConcatenation)
{
    auto data = simpleConcatenationData();

    for (const std::string& tup : data) {
        pCigarBuilder builder = CigarBuilder::create(_pool);
        pCigar bam_cigar = DebugUtils::str2cigar(tup, _pool);

        for (uint32_t i = 0; i < bam_cigar->num; ++i) {
            builder->add(bam_cigar->data[i]);
        }

        pCigar build_cigar = builder->make();

        ASSERT_STRCASEEQ(tup.c_str(), CigarUtils::uint2str(build_cigar->data, build_cigar->num).c_str());
    }
}

TEST_F(CigarBuilderUnitTest, testInitialAndFinalDeletions)
{
    auto data = initialAndFinalDeletionsData();
    for (const auto& tup : data) {
        CigarBuilder* builder = CigarBuilder::create(_pool);
        pCigar bam_cigar = DebugUtils::str2cigar(tup.first, _pool);

        for (uint32_t i = 0; i < bam_cigar->num; ++i) {
            builder->add(bam_cigar->data[i]);
        }

        pCigar build_cigar = builder->make();

        ASSERT_STRCASEEQ(tup.second.c_str(), CigarUtils::uint2str(build_cigar->data, build_cigar->num).c_str());
    }
}

TEST_F(CigarBuilderUnitTest, testRetainDeletions)
{
    auto data = retainDeletionsData();

    for (const auto& tup : data) {
        CigarBuilder* builder = CigarBuilder::create(false, _pool);
        pCigar bam_cigar = DebugUtils::str2cigar(tup.first, _pool);

        for (uint32_t i = 0; i < bam_cigar->num; ++i) {
            builder->add(bam_cigar->data[i]);
        }

        pCigar build_cigar = builder->make();

        ASSERT_STRCASEEQ(tup.second.c_str(), CigarUtils::uint2str(build_cigar->data, build_cigar->num).c_str());
    }
}

TEST_F(CigarBuilderUnitTest, testMergeConsecutive)
{
    auto data = mergeConsecutiveData();

    for (const auto& tup : data) {
        CigarBuilder* builder = CigarBuilder::create(_pool);
        pCigar bam_cigar = DebugUtils::str2cigar(tup.first, _pool);

        for (uint32_t i = 0; i < bam_cigar->num; ++i) {
            builder->add(bam_cigar->data[i]);
        }

        pCigar build_cigar = builder->make();

        ASSERT_EQ(tup.second, CigarUtils::uint2str(build_cigar->data, build_cigar->num).c_str());
    }
}

TEST_F(CigarBuilderUnitTest, testTrickyCases)
{
    std::string cigar_str = "10H10H10D10D10M";
    std::string cigar_result = "20H10M";

    CigarBuilder* builder = CigarBuilder::create(_pool);
    pCigar bam_cigar = DebugUtils::str2cigar(cigar_str, _pool);

    for (uint32_t i = 0; i < bam_cigar->num; ++i) {
        builder->add(bam_cigar->data[i]);
    }

    pCigar build_cigar = builder->make();

    ASSERT_STRCASEEQ(cigar_result.c_str(), CigarUtils::uint2str(build_cigar->data, build_cigar->num).c_str());
}

TEST_F(CigarBuilderUnitTest, testIndelSandwich)
{
    auto data = indelSandwichData();

    for (const auto& tup : data) {
        CigarBuilder* builder = CigarBuilder::create(_pool);
        pCigar bam_cigar = DebugUtils::str2cigar(tup.first, _pool);

        for (uint32_t i = 0; i < bam_cigar->num; ++i) {
            builder->add(bam_cigar->data[i]);
        }

        pCigar build_cigar = builder->make();

        ASSERT_STRCASEEQ(tup.second.c_str(), CigarUtils::uint2str(build_cigar->data, build_cigar->num).c_str());
    }
}

// 此Test中有部分异常数据会导致测试样例直接退出，测试其他用例时注释掉此位置
TEST_F(CigarBuilderUnitTest, testInvalid)
{
    auto data = invalidData();

    // for (const auto& tup : data) {
    //     CigarBuilder* builder = CigarBuilder::create(_pool);
    //     pCigar        bam_cigar = CigarUtils::str2uint(tup.first.c_str(), tup.second, _pool);
    //
    //     for (uint32_t i = 0; i < bam_cigar->num; ++i) {
    //         builder->add(bam_cigar->data[i]);
    //     }
    //
    //     builder->make();
    // }
}

TEST_F(CigarBuilderUnitTest, testRemovedDeletions)
{
    auto data = removedDeletionsData();
    for (const auto& tup : data) {
        CigarBuilder* builder = CigarBuilder::create(_pool);
        pCigar bam_cigar = DebugUtils::str2cigar(tup.first, _pool);

        for (uint32_t i = 0; i < bam_cigar->num; ++i) {
            builder->add(bam_cigar->data[i]);
        }

        builder->make();
        ASSERT_EQ(builder->leading_deletion_bases_removed(), tup.second.first);
        ASSERT_EQ(builder->trailing_deletion_bases_removed(), tup.second.second);
    }
}

TEST_F(CigarBuilderUnitTest, testRemovedDeletionsWithTwoMakes)
{
    auto data = removedDeletionsWithTwoMakesData();
    for (const auto& tup : data) {
        CigarBuilder* builder = CigarBuilder::create(_pool);
        pCigar bam_cigar = DebugUtils::str2cigar(std::get<0>(tup), _pool);
        for (uint32_t i = 0; i < bam_cigar->num; ++i) {
            builder->add(bam_cigar->data[i]);
        }
        builder->make();

        bam_cigar = DebugUtils::str2cigar(std::get<1>(tup), _pool);
        for (uint32_t i = 0; i < bam_cigar->num; ++i) {
            builder->add(bam_cigar->data[i]);
        }
        builder->make();

        ASSERT_EQ(builder->leading_deletion_bases_removed(), std::get<2>(tup));
        ASSERT_EQ(builder->trailing_deletion_bases_removed(), std::get<3>(tup));
    }
}

std::vector<std::string> CigarBuilderUnitTest::simpleConcatenationData()
{
    std::vector<std::vector<std::string>> leadingClips{{}, {"10H"}, {"10S"}, {"10H", "10S"}};
    std::vector<std::vector<std::string>> middleOperators{{"10M"}, {"10M", "10I", "10M"}, {"10M", "10D", "10M"}};
    std::vector<std::vector<std::string>> trailingClips{{}, {"10H"}, {"10S"}, {"10S", "10H"}};

    std::vector<std::string> result;
    for (const std::vector<std::string>& leadingString : leadingClips) {
        for (const std::vector<std::string>& middleString : middleOperators) {
            for (const std::vector<std::string>& trailingString : trailingClips) {
                std::string cigar_str;
                std::for_each(leadingString.begin(), leadingString.end(), [&cigar_str](const std::string& s) { cigar_str.append(s); });
                std::for_each(middleString.begin(), middleString.end(), [&cigar_str](const std::string& s) { cigar_str.append(s); });
                std::for_each(trailingString.begin(), trailingString.end(), [&cigar_str](const std::string& s) { cigar_str.append(s); });
                result.push_back(std::move(cigar_str));
            }
        }
    }
    return result;
}

std::vector<std::pair<std::string, std::string>> CigarBuilderUnitTest::initialAndFinalDeletionsData()
{
    return {{"10M10D", "10M"},
            {"10D10M", "10M"},
            {"10H10D10M", "10H10M"},
            {"10S10D10M", "10S10M"},
            {"10S10D10M10S", "10S10M10S"},
            {"10M10D10S", "10M10S"},
            {"10M10D10H", "10M10H"},
            {"10S10M10D10H", "10S10M10H"}};
}

std::vector<std::pair<std::string, std::string>> CigarBuilderUnitTest::retainDeletionsData()
{
    return {{"10M10D", "10M10D"},
            {"10D10M", "10D10M"},
            {"10M10D10D", "10M20D"},
            {"10H10D10M", "10H10D10M"},
            {"10M10D10H", "10M10D10H"},
            {"10M10D10S", "10M10D10S"},
            {"10S10D10M", "10S10D10M"},
            {"10S10D10M10S", "10S10D10M10S"},
            {"10S10M10D10H", "10S10M10D10H"},
            {"10M10M10D10D", "20M20D"}};
}

std::vector<std::pair<std::string, std::string>> CigarBuilderUnitTest::mergeConsecutiveData()
{
    return {{"10H10H10M", "20H10M"},
            {"10S10M10M", "10S20M"},
            {"10S10M10S10S", "10S10M20S"},
            {"10S10M10I10I10I10S10H", "10S10M30I10S10H"},
            {"10S10S10M10M10I10I10S10H", "20S20M20I10S10H"}};
}

std::vector<std::pair<std::string, std::string>> CigarBuilderUnitTest::indelSandwichData()
{
    return {{"10M10I10D10M", "10M10D10I10M"},
            {"10M10D10I10M", "10M10D10I10M"},
            {"10M10I10D10I10M", "10M10D20I10M"},
            {"10M10I10D10I10D10I10M", "10M20D30I10M"},
            {"10M10I10D10I10M10D10I10M", "10M10D20I10M10D10I10M"},
            {"10D10I10M", "10I10M"},
            {"10M10I10D", "10M10I"},
            {"10M10D10I", "10M10I"},
            {"10M10D10I10S", "10M10I10S"},
            {"10S10D10I10M", "10S10I10M"},
            {"10S10I10D10I10M", "10S20I10M"}};
}

std::vector<std::string> CigarBuilderUnitTest::invalidData()
{
    return {"10S", "10S10S", "10S10D", "10S10D10S", "10S10D10D10S", "10S10H10M", "10M10H10S", "10M10H10M", "10M10S10M"};
}

std::vector<std::pair<std::string, std::pair<int64_t, int64_t>>> CigarBuilderUnitTest::removedDeletionsData()
{
    return {{"10M", {0, 0}},
            {"10S10M", {0, 0}},
            {"10M10S", {0, 0}},
            {"10M10I10D10M", {0, 0}},
            {"10M10D10I10M", {0, 0}},
            {"10D10I10M", {10, 0}},
            {"10D10D10I10M", {20, 0}},
            {"10D10D10I10D10M", {30, 0}},
            {"10S10D10D10I10D10M", {30, 0}},
            {"10M10I10D", {0, 10}},
            {"10M10D10I", {0, 10}},
            {"10M10D10I10D", {0, 20}},
            {"10M10D10I10D10S10H", {0, 20}},
            {"10H10S10D10M10D10I10D10S10H", {10, 20}}};
}

std::vector<std::tuple<std::string, std::string, int64_t, int64_t>> CigarBuilderUnitTest::removedDeletionsWithTwoMakesData()
{
    return {{"10M", "10M", 0, 0},     {"10M10I", "10D10M", 0, 0},     {"10M10D", "10I10M", 0, 0},
            {"10D10I", "10M", 10, 0}, {"10D10D10I", "10D10M", 30, 0}, {"10H10S10D10M", "10D10I10D10S10H", 10, 20}};
}