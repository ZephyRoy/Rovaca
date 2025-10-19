#include <gtest/gtest.h>

#include <tuple>

#include "genotype_struct.h"
#include "htslib/sam.h"
#include "utils/cigar_utils.h"
#include "utils/read_record_utils.h"

using namespace rovaca;
typedef std::pair<std::string, uint32_t> CigarPair;

static constexpr size_t s_buffer_size = 1024 * 1024 * 10;

class ReadRecordUtilsUnitTest : public ::testing::Test
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

    static std::vector<std::tuple<CigarPair, int64_t, int64_t, int64_t, uint32_t>> testGetReadCoordinateForReferenceCoordinateData();

    uint8_t* _buffer{};
    pMemoryPool _pool{};
};

TEST_F(ReadRecordUtilsUnitTest, testGetReadCoordinateForReferenceCoordinate)
{
    uint32_t op;
    int64_t start, ref_cood, expected;
    auto data = testGetReadCoordinateForReferenceCoordinateData();
    for (const auto& tup : data) {
        const CigarPair& cigar = std::get<0>(tup);
        pCigar bam_cigar = CigarUtils::str2uint(cigar.first.c_str(), cigar.second, _pool);
        start = std::get<1>(tup), ref_cood = std::get<2>(tup), expected = std::get<3>(tup);
        op = std::get<4>(tup);
        auto result = ReadRecordUtils::get_read_index_for_reference_coordinate(start, bam_cigar->num, bam_cigar->data, ref_cood);
        ASSERT_EQ(result.first, expected);
        ASSERT_EQ(result.second, op);
    }
}

std::vector<std::tuple<CigarPair, int64_t, int64_t, int64_t, uint32_t>>
ReadRecordUtilsUnitTest::testGetReadCoordinateForReferenceCoordinateData()
{
    return {{{"10M", 1}, 1, 1, 0, BAM_CMATCH},
            {{"10M", 1}, 5, 5, 0, BAM_CMATCH},
            {{"10M", 1}, 1, 10, 9, BAM_CMATCH},
            {{"10M", 1}, 1, 5, 4, BAM_CMATCH},
            {{"10M", 1}, 1, 0, ReadRecordUtils::s_read_index_not_found, BAM_CUNINITIALIZE},
            {{"10M", 1}, 1, 11, ReadRecordUtils::s_read_index_not_found, BAM_CUNINITIALIZE},
            {{"5M5D5M", 3}, 1, 1, 0, BAM_CMATCH},
            {{"5M5D5M", 3}, 1, 5, 4, BAM_CMATCH},
            {{"5M5D5M", 3}, 1, 6, 5, BAM_CDEL},
            {{"5M5D5M", 3}, 1, 10, 5, BAM_CDEL},
            {{"5M5D5M", 3}, 1, 11, 5, BAM_CMATCH},
            {{"5M5D5M", 3}, 1, 15, 9, BAM_CMATCH},
            {{"5M5D5M", 3}, 1, 16, ReadRecordUtils::s_read_index_not_found, BAM_CUNINITIALIZE},
            {{"5M5I5M", 3}, 1, 1, 0, BAM_CMATCH},
            {{"5M5I5M", 3}, 1, 5, 4, BAM_CMATCH},
            {{"5M5I5M", 3}, 1, 6, 10, BAM_CMATCH},
            {{"5M5I5M", 3}, 1, 10, 14, BAM_CMATCH},
            {{"5M5I5M", 3}, 1, 11, ReadRecordUtils::s_read_index_not_found, BAM_CUNINITIALIZE}};
}