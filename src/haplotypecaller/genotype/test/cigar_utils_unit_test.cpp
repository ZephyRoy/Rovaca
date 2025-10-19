#include <gtest/gtest.h>

#include <tuple>
#include <vector>

#include "genotype_struct.h"
#include "htslib/sam.h"
#include "utils/cigar_utils.h"

using namespace rovaca;
typedef std::pair<std::string, uint32_t> CigarPair;

static constexpr size_t s_buffer_size = 1024 * 1024 * 10;

class CigarUtilsUnitTest : public ::testing::Test
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

    static std::vector<std::tuple<CigarPair, int64_t, int64_t, std::string, std::string>> testClipCigarData();
    static std::vector<std::tuple<CigarPair, int64_t, int64_t>> testAlignmentStartShiftData();
    static std::vector<std::tuple<CigarPair, uint32_t, uint32_t, uint32_t>> testCountRefBasesBasedOnCigarData();

    uint8_t* _buffer{};
    pMemoryPool _pool{};
};

TEST_F(CigarUtilsUnitTest, testClipCigar)
{
    auto data = testClipCigarData();
    int64_t start, stop;
    pCigar bam_cigar;
    for (const auto& tup : data) {
        const CigarPair& pair = std::get<0>(tup);
        bam_cigar = CigarUtils::str2uint(pair.first.c_str(), pair.second, _pool);

        start = std::get<1>(tup), stop = std::get<2>(tup);
        const std::string& soft_result = std::get<3>(tup);
        const std::string& hard_result = std::get<4>(tup);

        pCigar soft_clip = CigarUtils::clip_cigar(bam_cigar->data, bam_cigar->num, start, stop, BAM_CSOFT_CLIP, _pool);
        ASSERT_STRCASEEQ(soft_result.c_str(), CigarUtils::uint2str(soft_clip->data, soft_clip->num).c_str());

        pCigar hard_clip = CigarUtils::clip_cigar(bam_cigar->data, bam_cigar->num, start, stop, BAM_CHARD_CLIP, _pool);
        ASSERT_STRCASEEQ(hard_result.c_str(), CigarUtils::uint2str(hard_clip->data, hard_clip->num).c_str());
    }
}

TEST_F(CigarUtilsUnitTest, testAlignmentStartShift)
{
    auto data = testAlignmentStartShiftData();
    int64_t num_clips, result;
    pCigar bam_cigar;
    for (const auto& tup : data) {
        const CigarPair& pair = std::get<0>(tup);
        num_clips = std::get<1>(tup), result = std::get<2>(tup);
        bam_cigar = CigarUtils::str2uint(pair.first.c_str(), pair.second, _pool);
        ASSERT_EQ(result, CigarUtils::alignment_start_shift(bam_cigar->data, bam_cigar->num, num_clips));
    }
}

TEST_F(CigarUtilsUnitTest, testCountRefBasesBasedOnCigar)
{
    auto data = testCountRefBasesBasedOnCigarData();
    pCigar cigar;
    uint32_t start, stop, expected;
    for (size_t i = 0, len = data.size(); i < len; ++i) {
        const auto& tup = data.at(i);
        const CigarPair& pair = std::get<0>(tup);
        start = std::get<1>(tup);
        stop = std::get<2>(tup);
        expected = std::get<3>(tup);
        cigar = CigarUtils::str2uint(pair.first.c_str(), pair.second, _pool);
        ASSERT_EQ(expected, CigarUtils::count_ref_bases_and_clips(cigar->data, cigar->num, start, stop)) << "i=" << i;
    }
}

std::vector<std::tuple<CigarPair, int64_t, int64_t, std::string, std::string>> CigarUtilsUnitTest::testClipCigarData()
{
    return {// simple cases
            {{"10M", 1}, 0, 5, "5S5M", "5H5M"},
            {{"10M", 1}, 5, 10, "5M5S", "5M5H"},
            {{"10H10M", 2}, 0, 5, "10H5S5M", "15H5M"},
            {{"10H10M", 2}, 5, 10, "10H5M5S", "10H5M5H"},
            {{"10M10H", 2}, 0, 5, "5S5M10H", "5H5M10H"},

            // clipping into insertion
            {{"10M10I10M", 3}, 0, 5, "5S5M10I10M", "5H5M10I10M"},
            {{"10M10I10M", 3}, 0, 15, "15S5I10M", "15H5I10M"},
            {{"10M10I10M", 3}, 15, 30, "10M5I15S", "10M5I15H"},

            // clipping into a soft clip
            {{"10S10M10S", 3}, 0, 5, "10S10M10S", "5H5S10M10S"},
            {{"10S10M10S", 3}, 25, 30, "10S10M10S", "10S10M5S5H"},
            {{"10S10M10S", 3}, 0, 15, "15S5M10S", "15H5M10S"},

            // clipping over a deletion
            {{"10M10D10M", 3}, 0, 10, "10S10M", "10H10M"},
            {{"10M10D10M", 3}, 0, 15, "15S5M", "15H5M"},
            {{"10M10D10M", 3}, 5, 20, "5M15S", "5M15H"},

            // removing leading deletions
            {{"10D10M", 2}, 0, 5, "5S5M", "5H5M"}};
}

std::vector<std::tuple<CigarPair, int64_t, int64_t>> CigarUtilsUnitTest::testAlignmentStartShiftData()
{
    return {{{"70M", 1}, 10, 10},       {{"70M", 1}, 0, 0},         {{"30M10D30M", 3}, 29, 29}, {{"30M10D30M", 3}, 30, 40},
            {{"30M10D30M", 3}, 31, 41}, {{"30M10D30M", 3}, 29, 29}, {{"30M10I30M", 3}, 30, 30}, {{"30M10I30M", 3}, 31, 30},
            {{"30M10I30M", 3}, 40, 30}, {{"30M10I30M", 3}, 41, 31}, {{"10H10M", 2}, 5, 5},      {{"10S10M", 2}, 5, 0},
            {{"10S10M", 2}, 5, 0}};
}

std::vector<std::tuple<CigarPair, uint32_t, uint32_t, uint32_t>> CigarUtilsUnitTest::testCountRefBasesBasedOnCigarData()
{
    return {{{"10M", 1}, 0, 1, 10},         {{"10M1D", 2}, 0, 1, 10},     {{"10M1D1S", 3}, 0, 1, 10},     {{"10M1D1N1S", 4}, 0, 1, 10},
            {{"10M1D1N1S1H", 5}, 0, 1, 10}, {{"10M1D", 2}, 0, 2, 11},     {{"10M1D1S", 3}, 0, 2, 11},     {{"10M1D1N1S", 4}, 0, 2, 11},
            {{"10M1D1N1S1H", 5}, 0, 2, 11}, {{"10M1=1X1S", 4}, 0, 2, 11}, {{"10M1=1X1S1H", 5}, 0, 2, 11}, {{"10M1D2N4S", 4}, 1, 3, 3},
            {{"10M1D2N4S8H", 5}, 1, 4, 7},  {{"10M1I2N4S", 4}, 1, 3, 2},  {{"10M1I2N4S8H", 5}, 1, 4, 6},  {{"1M1S1H1H", 4}, 1, 2, 1},
            {{"1M1I1H1H", 4}, 1, 2, 0},     {{"1M1P1H1H", 4}, 1, 2, 0},   {{"1M1H1H1H", 4}, 1, 2, 1}};
}