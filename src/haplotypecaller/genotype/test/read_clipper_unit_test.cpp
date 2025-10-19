#include <gtest/gtest.h>

#include <vector>

#include "bam_data_pool.hpp"
#include "cigar_builder.h"
#include "genotype_struct.h"
#include "htslib/sam.h"
#include "read_clipper.h"
#include "read_record.h"

using namespace rovaca;
typedef std::pair<std::string, uint32_t> CigarPair;

static constexpr size_t s_buffer_size = 1024 * 1024 * 10;

class ReadClipperUnitTest : public ::testing::Test
{
protected:
    void SetUp() override
    {
        _buffer = new uint8_t[s_buffer_size]{};
        _pool = new std::pmr::monotonic_buffer_resource(_buffer, s_buffer_size, std::pmr::null_memory_resource());
        _bampool = new BamDataPool(s_buffer_size);
    }
    void TearDown() override
    {
        delete _bampool;
        delete _pool;
        delete[] _buffer;
    }

    std::vector<rovaca::pReadRecord> generateReadList();

    static void assertRefAlignmentConsistent(rovaca::pReadRecord read);
    static void assertReadLengthConsistent(rovaca::pReadRecord read);

    uint8_t *_buffer{};
    pMemoryPool _pool{};
    pBamDataPool _bampool{};
};

TEST_F(ReadClipperUnitTest, testHardClipBothEndsByReferenceCoordinates)
{
    // auto data = generateReadList();
    // for (rovaca::pReadRecord read : data) {
    //     std::cout << read->seq_i(0) << std::endl;
    // }
}

TEST_F(ReadClipperUnitTest, testHardClipByReferenceCoordinates)
{
    auto data = generateReadList();

    int64_t i, start, end;
    rovaca::pReadRecord clipLeft, clipRight;
    size_t idx;
    for (idx = 0; idx < data.size(); ++idx) {
        rovaca::pReadRecord read = data.at(idx);
        start = read->get_soft_start();
        end = read->get_soft_end();
        for (i = start; i <= end; ++i) {
            clipLeft = clipRight = nullptr;
            clipLeft = rovaca::ReadClipper(read, _pool, _bampool).hard_clip_by_reference_coordinates_left_tail(i);
            if (clipLeft) {
                ASSERT_TRUE(clipLeft->get_start() >= std::min(read->get_stop(), i + 1)) << idx;
                assertRefAlignmentConsistent(clipLeft);
                assertReadLengthConsistent(clipLeft);
            }

            clipRight = rovaca::ReadClipper(read, _pool, _bampool).hard_clip_by_reference_coordinates_right_tail(i);
            if (clipRight && clipRight->get_start() <= clipRight->get_stop()) {
                ASSERT_TRUE(clipRight->get_stop() <= std::max(read->get_start(), i - 1)) << idx;
                assertRefAlignmentConsistent(clipRight);
                assertReadLengthConsistent(clipRight);
            }
        }
    }
}

std::vector<rovaca::pReadRecord> ReadClipperUnitTest::generateReadList()
{
    // clang-format off
    const char sam[] = "data:,"
        "@SQ\tSN:chr1\tLN:2492\n"
        "rr\t0\tchr1\t10000\t0\t1M1I1M1I1M1I\t*\t0\t0\tACTGAC\t??????\n"
        "rr\t0\tchr1\t10000\t0\t1M1D1M1I1M1I\t*\t0\t0\tACTGA\t?????\n"
        "rr\t0\tchr1\t10000\t0\t1I2M1I1M1I\t*\t0\t0\tACTGAC\t??????\n"
        "rr\t0\tchr1\t10000\t0\t3M1I1M1I\t*\t0\t0\tACTGAC\t??????\n"
        "rr\t0\tchr1\t10000\t0\t1M1I1M1D1M1I\t*\t0\t0\tACTGA\t?????\n"
        "rr\t0\tchr1\t10000\t0\t1M1D1M1D1M1I\t*\t0\t0\tACTG\t????\n"
        "rr\t0\tchr1\t10000\t0\t1I2M1D1M1I\t*\t0\t0\tACTGA\t?????\n"
        "rr\t0\tchr1\t10000\t0\t3M1D1M1I\t*\t0\t0\tACTGA\t?????\n"
        "rr\t0\tchr1\t10000\t0\t1I1M1I2M1I\t*\t0\t0\tACTGAC\t??????\n"
        "rr\t0\tchr1\t10000\t0\t2M1I2M1I\t*\t0\t0\tACTGAC\t??????\n"
        "rr\t0\tchr1\t10000\t0\t1I1M1D2M1I\t*\t0\t0\tACTGA\t?????\n"
        "rr\t0\tchr1\t10000\t0\t2M1D2M1I\t*\t0\t0\tACTGA\t?????\n"
        "rr\t0\tchr1\t10000\t0\t1M1I3M1I\t*\t0\t0\tACTGAC\t??????\n"
        "rr\t0\tchr1\t10000\t0\t1M1D3M1I\t*\t0\t0\tACTGA\t?????\n"
        "rr\t0\tchr1\t10000\t0\t1I4M1I\t*\t0\t0\tACTGAC\t??????\n"
        "rr\t0\tchr1\t10000\t0\t5M1I\t*\t0\t0\tACTGAC\t??????\n"
        "rr\t0\tchr1\t10000\t0\t1I1M1I1M1I1M\t*\t0\t0\tACTGAC\t??????\n"
        "rr\t0\tchr1\t10000\t0\t2M1I1M1I1M\t*\t0\t0\tACTGAC\t??????\n"
        "rr\t0\tchr1\t10000\t0\t1I1M1D1M1I1M\t*\t0\t0\tACTGA\t?????\n"
        "rr\t0\tchr1\t10000\t0\t2M1D1M1I1M\t*\t0\t0\tACTGA\t?????\n"
        "rr\t0\tchr1\t10000\t0\t1M1I2M1I1M\t*\t0\t0\tACTGAC\t??????\n"
        "rr\t0\tchr1\t10000\t0\t1M1D2M1I1M\t*\t0\t0\tACTGA\t?????\n"
        "rr\t0\tchr1\t10000\t0\t1I3M1I1M\t*\t0\t0\tACTGAC\t??????\n"
        "rr\t0\tchr1\t10000\t0\t4M1I1M\t*\t0\t0\tACTGAC\t??????\n"
        "rr\t0\tchr1\t10000\t0\t1I1M1I1M1D1M\t*\t0\t0\tACTGA\t?????\n"
        "rr\t0\tchr1\t10000\t0\t2M1I1M1D1M\t*\t0\t0\tACTGA\t?????\n"
        "rr\t0\tchr1\t10000\t0\t1I1M1D1M1D1M\t*\t0\t0\tACTG\t????\n"
        "rr\t0\tchr1\t10000\t0\t2M1D1M1D1M\t*\t0\t0\tACTG\t????\n"
        "rr\t0\tchr1\t10000\t0\t1M1I2M1D1M\t*\t0\t0\tACTGA\t?????\n"
        "rr\t0\tchr1\t10000\t0\t1M1D2M1D1M\t*\t0\t0\tACTG\t????\n"
        "rr\t0\tchr1\t10000\t0\t1I3M1D1M\t*\t0\t0\tACTGA\t?????\n"
        "rr\t0\tchr1\t10000\t0\t4M1D1M\t*\t0\t0\tACTGA\t?????\n"
        "rr\t0\tchr1\t10000\t0\t1M1I1M1I2M\t*\t0\t0\tACTGAC\t??????\n"
        "rr\t0\tchr1\t10000\t0\t1M1D1M1I2M\t*\t0\t0\tACTGA\t?????\n"
        "rr\t0\tchr1\t10000\t0\t1I2M1I2M\t*\t0\t0\tACTGAC\t??????\n"
        "rr\t0\tchr1\t10000\t0\t3M1I2M\t*\t0\t0\tACTGAC\t??????\n"
        "rr\t0\tchr1\t10000\t0\t1M1I1M1D2M\t*\t0\t0\tACTGA\t?????\n"
        "rr\t0\tchr1\t10000\t0\t1M1D1M1D2M\t*\t0\t0\tACTG\t????\n"
        "rr\t0\tchr1\t10000\t0\t1I2M1D2M\t*\t0\t0\tACTGA\t?????\n"
        "rr\t0\tchr1\t10000\t0\t3M1D2M\t*\t0\t0\tACTGA\t?????\n"
        "rr\t0\tchr1\t10000\t0\t1I1M1I3M\t*\t0\t0\tACTGAC\t??????\n"
        "rr\t0\tchr1\t10000\t0\t2M1I3M\t*\t0\t0\tACTGAC\t??????\n"
        "rr\t0\tchr1\t10000\t0\t1I1M1D3M\t*\t0\t0\tACTGA\t?????\n"
        "rr\t0\tchr1\t10000\t0\t2M1D3M\t*\t0\t0\tACTGA\t?????\n"
        "rr\t0\tchr1\t10000\t0\t1M1I4M\t*\t0\t0\tACTGAC\t??????\n"
        "rr\t0\tchr1\t10000\t0\t1M1D4M\t*\t0\t0\tACTGA\t?????\n"
        "rr\t0\tchr1\t10000\t0\t1I5M\t*\t0\t0\tACTGAC\t??????\n"
        "rr\t0\tchr1\t10000\t0\t6M\t*\t0\t0\tACTGAC\t??????\n"
        "rr\t0\tchr1\t10000\t0\t1I1M1I1M1I1H\t*\t0\t0\tACTGA\t?????\n"
        "rr\t0\tchr1\t10000\t0\t2M1I1M1I1H\t*\t0\t0\tACTGA\t?????\n"
        "rr\t0\tchr1\t10000\t0\t1I1M1D1M1I1H\t*\t0\t0\tACTG\t????\n"
        "rr\t0\tchr1\t10000\t0\t2M1D1M1I1H\t*\t0\t0\tACTG\t????\n"
        "rr\t0\tchr1\t10000\t0\t1M1I2M1I1H\t*\t0\t0\tACTGA\t?????\n"
        "rr\t0\tchr1\t10000\t0\t1M1D2M1I1H\t*\t0\t0\tACTG\t????\n"
        "rr\t0\tchr1\t10000\t0\t1I3M1I1H\t*\t0\t0\tACTGA\t?????\n"
        "rr\t0\tchr1\t10000\t0\t4M1I1H\t*\t0\t0\tACTGA\t?????\n"
        "rr\t0\tchr1\t10000\t0\t1M1I1M1I1M1H\t*\t0\t0\tACTGA\t?????\n"
        "rr\t0\tchr1\t10000\t0\t1M1D1M1I1M1H\t*\t0\t0\tACTG\t????\n"
        "rr\t0\tchr1\t10000\t0\t1I2M1I1M1H\t*\t0\t0\tACTGA\t?????\n"
        "rr\t0\tchr1\t10000\t0\t3M1I1M1H\t*\t0\t0\tACTGA\t?????\n"
        "rr\t0\tchr1\t10000\t0\t1M1I1M1D1M1H\t*\t0\t0\tACTG\t????\n"
        "rr\t0\tchr1\t10000\t0\t1M1D1M1D1M1H\t*\t0\t0\tACT\t???\n"
        "rr\t0\tchr1\t10000\t0\t1I2M1D1M1H\t*\t0\t0\tACTG\t????\n"
        "rr\t0\tchr1\t10000\t0\t3M1D1M1H\t*\t0\t0\tACTG\t????\n"
        "rr\t0\tchr1\t10000\t0\t1I1M1I2M1H\t*\t0\t0\tACTGA\t?????\n"
        "rr\t0\tchr1\t10000\t0\t2M1I2M1H\t*\t0\t0\tACTGA\t?????\n"
        "rr\t0\tchr1\t10000\t0\t1I1M1D2M1H\t*\t0\t0\tACTG\t????\n"
        "rr\t0\tchr1\t10000\t0\t2M1D2M1H\t*\t0\t0\tACTG\t????\n"
        "rr\t0\tchr1\t10000\t0\t1M1I3M1H\t*\t0\t0\tACTGA\t?????\n"
        "rr\t0\tchr1\t10000\t0\t1M1D3M1H\t*\t0\t0\tACTG\t????\n"
        "rr\t0\tchr1\t10000\t0\t1I4M1H\t*\t0\t0\tACTGA\t?????\n"
        "rr\t0\tchr1\t10000\t0\t5M1H\t*\t0\t0\tACTGA\t?????\n"
        "rr\t0\tchr1\t10000\t0\t1I1M1I1M1I1S\t*\t0\t0\tACTGAC\t??????\n"
        "rr\t0\tchr1\t10000\t0\t2M1I1M1I1S\t*\t0\t0\tACTGAC\t??????\n"
        "rr\t0\tchr1\t10000\t0\t1I1M1D1M1I1S\t*\t0\t0\tACTGA\t?????\n"
        "rr\t0\tchr1\t10000\t0\t2M1D1M1I1S\t*\t0\t0\tACTGA\t?????\n"
        "rr\t0\tchr1\t10000\t0\t1M1I2M1I1S\t*\t0\t0\tACTGAC\t??????\n"
        "rr\t0\tchr1\t10000\t0\t1M1D2M1I1S\t*\t0\t0\tACTGA\t?????\n"
        "rr\t0\tchr1\t10000\t0\t1I3M1I1S\t*\t0\t0\tACTGAC\t??????\n"
        "rr\t0\tchr1\t10000\t0\t4M1I1S\t*\t0\t0\tACTGAC\t??????\n"
        "rr\t0\tchr1\t10000\t0\t1M1I1M1I1M1S\t*\t0\t0\tACTGAC\t??????\n"
        "rr\t0\tchr1\t10000\t0\t1M1D1M1I1M1S\t*\t0\t0\tACTGA\t?????\n"
        "rr\t0\tchr1\t10000\t0\t1I2M1I1M1S\t*\t0\t0\tACTGAC\t??????\n"
        "rr\t0\tchr1\t10000\t0\t3M1I1M1S\t*\t0\t0\tACTGAC\t??????\n"
        "rr\t0\tchr1\t10000\t0\t1M1I1M1D1M1S\t*\t0\t0\tACTGA\t?????\n"
        "rr\t0\tchr1\t10000\t0\t1M1D1M1D1M1S\t*\t0\t0\tACTG\t????\n"
        "rr\t0\tchr1\t10000\t0\t1I2M1D1M1S\t*\t0\t0\tACTGA\t?????\n"
        "rr\t0\tchr1\t10000\t0\t3M1D1M1S\t*\t0\t0\tACTGA\t?????\n"
        "rr\t0\tchr1\t10000\t0\t1I1M1I2M1S\t*\t0\t0\tACTGAC\t??????\n"
        "rr\t0\tchr1\t10000\t0\t2M1I2M1S\t*\t0\t0\tACTGAC\t??????\n"
        "rr\t0\tchr1\t10000\t0\t1I1M1D2M1S\t*\t0\t0\tACTGA\t?????\n"
        "rr\t0\tchr1\t10000\t0\t2M1D2M1S\t*\t0\t0\tACTGA\t?????\n"
        "rr\t0\tchr1\t10000\t0\t1M1I3M1S\t*\t0\t0\tACTGAC\t??????\n"
        "rr\t0\tchr1\t10000\t0\t1M1D3M1S\t*\t0\t0\tACTGA\t?????\n"
        "rr\t0\tchr1\t10000\t0\t1I4M1S\t*\t0\t0\tACTGAC\t??????\n"
        "rr\t0\tchr1\t10000\t0\t5M1S\t*\t0\t0\tACTGAC\t??????\n"
        "rr\t0\tchr1\t10000\t0\t1M1I1M1I1S1H\t*\t0\t0\tACTGA\t?????\n"
        "rr\t0\tchr1\t10000\t0\t1M1D1M1I1S1H\t*\t0\t0\tACTG\t????\n"
        "rr\t0\tchr1\t10000\t0\t1I2M1I1S1H\t*\t0\t0\tACTGA\t?????\n"
        "rr\t0\tchr1\t10000\t0\t3M1I1S1H\t*\t0\t0\tACTGA\t?????\n"
        "rr\t0\tchr1\t10000\t0\t1I1M1I1M1S1H\t*\t0\t0\tACTGA\t?????\n"
        "rr\t0\tchr1\t10000\t0\t2M1I1M1S1H\t*\t0\t0\tACTGA\t?????\n"
        "rr\t0\tchr1\t10000\t0\t1I1M1D1M1S1H\t*\t0\t0\tACTG\t????\n"
        "rr\t0\tchr1\t10000\t0\t2M1D1M1S1H\t*\t0\t0\tACTG\t????\n"
        "rr\t0\tchr1\t10000\t0\t1M1I2M1S1H\t*\t0\t0\tACTGA\t?????\n"
        "rr\t0\tchr1\t10000\t0\t1M1D2M1S1H\t*\t0\t0\tACTG\t????\n"
        "rr\t0\tchr1\t10000\t0\t1I3M1S1H\t*\t0\t0\tACTGA\t?????\n"
        "rr\t0\tchr1\t10000\t0\t4M1S1H\t*\t0\t0\tACTGA\t?????\n"
        "rr\t0\tchr1\t10000\t0\t1H1I1M1I1M1I\t*\t0\t0\tACTGA\t?????\n"
        "rr\t0\tchr1\t10000\t0\t1H2M1I1M1I\t*\t0\t0\tACTGA\t?????\n"
        "rr\t0\tchr1\t10000\t0\t1H1I1M1D1M1I\t*\t0\t0\tACTG\t????\n"
        "rr\t0\tchr1\t10000\t0\t1H2M1D1M1I\t*\t0\t0\tACTG\t????\n"
        "rr\t0\tchr1\t10000\t0\t1H1M1I2M1I\t*\t0\t0\tACTGA\t?????\n"
        "rr\t0\tchr1\t10000\t0\t1H1M1D2M1I\t*\t0\t0\tACTG\t????\n"
        "rr\t0\tchr1\t10000\t0\t1H1I3M1I\t*\t0\t0\tACTGA\t?????\n"
        "rr\t0\tchr1\t10000\t0\t1H4M1I\t*\t0\t0\tACTGA\t?????\n"
        "rr\t0\tchr1\t10000\t0\t1H1M1I1M1I1M\t*\t0\t0\tACTGA\t?????\n"
        "rr\t0\tchr1\t10000\t0\t1H1M1D1M1I1M\t*\t0\t0\tACTG\t????\n"
        "rr\t0\tchr1\t10000\t0\t1H1I2M1I1M\t*\t0\t0\tACTGA\t?????\n"
        "rr\t0\tchr1\t10000\t0\t1H3M1I1M\t*\t0\t0\tACTGA\t?????\n"
        "rr\t0\tchr1\t10000\t0\t1H1M1I1M1D1M\t*\t0\t0\tACTG\t????\n"
        "rr\t0\tchr1\t10000\t0\t1H1M1D1M1D1M\t*\t0\t0\tACT\t???\n"
        "rr\t0\tchr1\t10000\t0\t1H1I2M1D1M\t*\t0\t0\tACTG\t????\n"
        "rr\t0\tchr1\t10000\t0\t1H3M1D1M\t*\t0\t0\tACTG\t????\n"
        "rr\t0\tchr1\t10000\t0\t1H1I1M1I2M\t*\t0\t0\tACTGA\t?????\n"
        "rr\t0\tchr1\t10000\t0\t1H2M1I2M\t*\t0\t0\tACTGA\t?????\n"
        "rr\t0\tchr1\t10000\t0\t1H1I1M1D2M\t*\t0\t0\tACTG\t????\n"
        "rr\t0\tchr1\t10000\t0\t1H2M1D2M\t*\t0\t0\tACTG\t????\n"
        "rr\t0\tchr1\t10000\t0\t1H1M1I3M\t*\t0\t0\tACTGA\t?????\n"
        "rr\t0\tchr1\t10000\t0\t1H1M1D3M\t*\t0\t0\tACTG\t????\n"
        "rr\t0\tchr1\t10000\t0\t1H1I4M\t*\t0\t0\tACTGA\t?????\n"
        "rr\t0\tchr1\t10000\t0\t1H5M\t*\t0\t0\tACTGA\t?????\n"
        "rr\t0\tchr1\t10000\t0\t1H1M1I1M1I1H\t*\t0\t0\tACTG\t????\n"
        "rr\t0\tchr1\t10000\t0\t1H1M1D1M1I1H\t*\t0\t0\tACT\t???\n"
        "rr\t0\tchr1\t10000\t0\t1H1I2M1I1H\t*\t0\t0\tACTG\t????\n"
        "rr\t0\tchr1\t10000\t0\t1H3M1I1H\t*\t0\t0\tACTG\t????\n"
        "rr\t0\tchr1\t10000\t0\t1H1I1M1I1M1H\t*\t0\t0\tACTG\t????\n"
        "rr\t0\tchr1\t10000\t0\t1H2M1I1M1H\t*\t0\t0\tACTG\t????\n"
        "rr\t0\tchr1\t10000\t0\t1H1I1M1D1M1H\t*\t0\t0\tACT\t???\n"
        "rr\t0\tchr1\t10000\t0\t1H2M1D1M1H\t*\t0\t0\tACT\t???\n"
        "rr\t0\tchr1\t10000\t0\t1H1M1I2M1H\t*\t0\t0\tACTG\t????\n"
        "rr\t0\tchr1\t10000\t0\t1H1M1D2M1H\t*\t0\t0\tACT\t???\n"
        "rr\t0\tchr1\t10000\t0\t1H1I3M1H\t*\t0\t0\tACTG\t????\n"
        "rr\t0\tchr1\t10000\t0\t1H4M1H\t*\t0\t0\tACTG\t????\n"
        "rr\t0\tchr1\t10000\t0\t1H1M1I1M1I1S\t*\t0\t0\tACTGA\t?????\n"
        "rr\t0\tchr1\t10000\t0\t1H1M1D1M1I1S\t*\t0\t0\tACTG\t????\n"
        "rr\t0\tchr1\t10000\t0\t1H1I2M1I1S\t*\t0\t0\tACTGA\t?????\n"
        "rr\t0\tchr1\t10000\t0\t1H3M1I1S\t*\t0\t0\tACTGA\t?????\n"
        "rr\t0\tchr1\t10000\t0\t1H1I1M1I1M1S\t*\t0\t0\tACTGA\t?????\n"
        "rr\t0\tchr1\t10000\t0\t1H2M1I1M1S\t*\t0\t0\tACTGA\t?????\n"
        "rr\t0\tchr1\t10000\t0\t1H1I1M1D1M1S\t*\t0\t0\tACTG\t????\n"
        "rr\t0\tchr1\t10000\t0\t1H2M1D1M1S\t*\t0\t0\tACTG\t????\n"
        "rr\t0\tchr1\t10000\t0\t1H1M1I2M1S\t*\t0\t0\tACTGA\t?????\n"
        "rr\t0\tchr1\t10000\t0\t1H1M1D2M1S\t*\t0\t0\tACTG\t????\n"
        "rr\t0\tchr1\t10000\t0\t1H1I3M1S\t*\t0\t0\tACTGA\t?????\n"
        "rr\t0\tchr1\t10000\t0\t1H4M1S\t*\t0\t0\tACTGA\t?????\n"
        "rr\t0\tchr1\t10000\t0\t1H1I1M1I1S1H\t*\t0\t0\tACTG\t????\n"
        "rr\t0\tchr1\t10000\t0\t1H2M1I1S1H\t*\t0\t0\tACTG\t????\n"
        "rr\t0\tchr1\t10000\t0\t1H1M1I1M1S1H\t*\t0\t0\tACTG\t????\n"
        "rr\t0\tchr1\t10000\t0\t1H1M1D1M1S1H\t*\t0\t0\tACT\t???\n"
        "rr\t0\tchr1\t10000\t0\t1H1I2M1S1H\t*\t0\t0\tACTG\t????\n"
        "rr\t0\tchr1\t10000\t0\t1H3M1S1H\t*\t0\t0\tACTG\t????\n"
        "rr\t0\tchr1\t10000\t0\t1S1I1M1I1M1I\t*\t0\t0\tACTGAC\t??????\n"
        "rr\t0\tchr1\t10000\t0\t1S2M1I1M1I\t*\t0\t0\tACTGAC\t??????\n"
        "rr\t0\tchr1\t10000\t0\t1S1I1M1D1M1I\t*\t0\t0\tACTGA\t?????\n"
        "rr\t0\tchr1\t10000\t0\t1S2M1D1M1I\t*\t0\t0\tACTGA\t?????\n"
        "rr\t0\tchr1\t10000\t0\t1S1M1I2M1I\t*\t0\t0\tACTGAC\t??????\n"
        "rr\t0\tchr1\t10000\t0\t1S1M1D2M1I\t*\t0\t0\tACTGA\t?????\n"
        "rr\t0\tchr1\t10000\t0\t1S1I3M1I\t*\t0\t0\tACTGAC\t??????\n"
        "rr\t0\tchr1\t10000\t0\t1S4M1I\t*\t0\t0\tACTGAC\t??????\n"
        "rr\t0\tchr1\t10000\t0\t1S1M1I1M1I1M\t*\t0\t0\tACTGAC\t??????\n"
        "rr\t0\tchr1\t10000\t0\t1S1M1D1M1I1M\t*\t0\t0\tACTGA\t?????\n"
        "rr\t0\tchr1\t10000\t0\t1S1I2M1I1M\t*\t0\t0\tACTGAC\t??????\n"
        "rr\t0\tchr1\t10000\t0\t1S3M1I1M\t*\t0\t0\tACTGAC\t??????\n"
        "rr\t0\tchr1\t10000\t0\t1S1M1I1M1D1M\t*\t0\t0\tACTGA\t?????\n"
        "rr\t0\tchr1\t10000\t0\t1S1M1D1M1D1M\t*\t0\t0\tACTG\t????\n"
        "rr\t0\tchr1\t10000\t0\t1S1I2M1D1M\t*\t0\t0\tACTGA\t?????\n"
        "rr\t0\tchr1\t10000\t0\t1S3M1D1M\t*\t0\t0\tACTGA\t?????\n"
        "rr\t0\tchr1\t10000\t0\t1S1I1M1I2M\t*\t0\t0\tACTGAC\t??????\n"
        "rr\t0\tchr1\t10000\t0\t1S2M1I2M\t*\t0\t0\tACTGAC\t??????\n"
        "rr\t0\tchr1\t10000\t0\t1S1I1M1D2M\t*\t0\t0\tACTGA\t?????\n"
        "rr\t0\tchr1\t10000\t0\t1S2M1D2M\t*\t0\t0\tACTGA\t?????\n"
        "rr\t0\tchr1\t10000\t0\t1S1M1I3M\t*\t0\t0\tACTGAC\t??????\n"
        "rr\t0\tchr1\t10000\t0\t1S1M1D3M\t*\t0\t0\tACTGA\t?????\n"
        "rr\t0\tchr1\t10000\t0\t1S1I4M\t*\t0\t0\tACTGAC\t??????\n"
        "rr\t0\tchr1\t10000\t0\t1S5M\t*\t0\t0\tACTGAC\t??????\n"
        "rr\t0\tchr1\t10000\t0\t1S1M1I1M1I1H\t*\t0\t0\tACTGA\t?????\n"
        "rr\t0\tchr1\t10000\t0\t1S1M1D1M1I1H\t*\t0\t0\tACTG\t????\n"
        "rr\t0\tchr1\t10000\t0\t1S1I2M1I1H\t*\t0\t0\tACTGA\t?????\n"
        "rr\t0\tchr1\t10000\t0\t1S3M1I1H\t*\t0\t0\tACTGA\t?????\n"
        "rr\t0\tchr1\t10000\t0\t1S1I1M1I1M1H\t*\t0\t0\tACTGA\t?????\n"
        "rr\t0\tchr1\t10000\t0\t1S2M1I1M1H\t*\t0\t0\tACTGA\t?????\n"
        "rr\t0\tchr1\t10000\t0\t1S1I1M1D1M1H\t*\t0\t0\tACTG\t????\n"
        "rr\t0\tchr1\t10000\t0\t1S2M1D1M1H\t*\t0\t0\tACTG\t????\n"
        "rr\t0\tchr1\t10000\t0\t1S1M1I2M1H\t*\t0\t0\tACTGA\t?????\n"
        "rr\t0\tchr1\t10000\t0\t1S1M1D2M1H\t*\t0\t0\tACTG\t????\n"
        "rr\t0\tchr1\t10000\t0\t1S1I3M1H\t*\t0\t0\tACTGA\t?????\n"
        "rr\t0\tchr1\t10000\t0\t1S4M1H\t*\t0\t0\tACTGA\t?????\n"
        "rr\t0\tchr1\t10000\t0\t1S1M1I1M1I1S\t*\t0\t0\tACTGAC\t??????\n"
        "rr\t0\tchr1\t10000\t0\t1S1M1D1M1I1S\t*\t0\t0\tACTGA\t?????\n"
        "rr\t0\tchr1\t10000\t0\t1S1I2M1I1S\t*\t0\t0\tACTGAC\t??????\n"
        "rr\t0\tchr1\t10000\t0\t1S3M1I1S\t*\t0\t0\tACTGAC\t??????\n"
        "rr\t0\tchr1\t10000\t0\t1S1I1M1I1M1S\t*\t0\t0\tACTGAC\t??????\n"
        "rr\t0\tchr1\t10000\t0\t1S2M1I1M1S\t*\t0\t0\tACTGAC\t??????\n"
        "rr\t0\tchr1\t10000\t0\t1S1I1M1D1M1S\t*\t0\t0\tACTGA\t?????\n"
        "rr\t0\tchr1\t10000\t0\t1S2M1D1M1S\t*\t0\t0\tACTGA\t?????\n"
        "rr\t0\tchr1\t10000\t0\t1S1M1I2M1S\t*\t0\t0\tACTGAC\t??????\n"
        "rr\t0\tchr1\t10000\t0\t1S1M1D2M1S\t*\t0\t0\tACTGA\t?????\n"
        "rr\t0\tchr1\t10000\t0\t1S1I3M1S\t*\t0\t0\tACTGAC\t??????\n"
        "rr\t0\tchr1\t10000\t0\t1S4M1S\t*\t0\t0\tACTGAC\t??????\n"
        "rr\t0\tchr1\t10000\t0\t1S1I1M1I1S1H\t*\t0\t0\tACTGA\t?????\n"
        "rr\t0\tchr1\t10000\t0\t1S2M1I1S1H\t*\t0\t0\tACTGA\t?????\n"
        "rr\t0\tchr1\t10000\t0\t1S1M1I1M1S1H\t*\t0\t0\tACTGA\t?????\n"
        "rr\t0\tchr1\t10000\t0\t1S1M1D1M1S1H\t*\t0\t0\tACTG\t????\n"
        "rr\t0\tchr1\t10000\t0\t1S1I2M1S1H\t*\t0\t0\tACTGA\t?????\n"
        "rr\t0\tchr1\t10000\t0\t1S3M1S1H\t*\t0\t0\tACTGA\t?????\n"
        "rr\t0\tchr1\t10000\t0\t1H1S1M1I1M1I\t*\t0\t0\tACTGA\t?????\n"
        "rr\t0\tchr1\t10000\t0\t1H1S1M1D1M1I\t*\t0\t0\tACTG\t????\n"
        "rr\t0\tchr1\t10000\t0\t1H1S1I2M1I\t*\t0\t0\tACTGA\t?????\n"
        "rr\t0\tchr1\t10000\t0\t1H1S3M1I\t*\t0\t0\tACTGA\t?????\n"
        "rr\t0\tchr1\t10000\t0\t1H1S1I1M1I1M\t*\t0\t0\tACTGA\t?????\n"
        "rr\t0\tchr1\t10000\t0\t1H1S2M1I1M\t*\t0\t0\tACTGA\t?????\n"
        "rr\t0\tchr1\t10000\t0\t1H1S1I1M1D1M\t*\t0\t0\tACTG\t????\n"
        "rr\t0\tchr1\t10000\t0\t1H1S2M1D1M\t*\t0\t0\tACTG\t????\n"
        "rr\t0\tchr1\t10000\t0\t1H1S1M1I2M\t*\t0\t0\tACTGA\t?????\n"
        "rr\t0\tchr1\t10000\t0\t1H1S1M1D2M\t*\t0\t0\tACTG\t????\n"
        "rr\t0\tchr1\t10000\t0\t1H1S1I3M\t*\t0\t0\tACTGA\t?????\n"
        "rr\t0\tchr1\t10000\t0\t1H1S4M\t*\t0\t0\tACTGA\t?????\n"
        "rr\t0\tchr1\t10000\t0\t1H1S1I1M1I1H\t*\t0\t0\tACTG\t????\n"
        "rr\t0\tchr1\t10000\t0\t1H1S2M1I1H\t*\t0\t0\tACTG\t????\n"
        "rr\t0\tchr1\t10000\t0\t1H1S1M1I1M1H\t*\t0\t0\tACTG\t????\n"
        "rr\t0\tchr1\t10000\t0\t1H1S1M1D1M1H\t*\t0\t0\tACT\t???\n"
        "rr\t0\tchr1\t10000\t0\t1H1S1I2M1H\t*\t0\t0\tACTG\t????\n"
        "rr\t0\tchr1\t10000\t0\t1H1S3M1H\t*\t0\t0\tACTG\t????\n"
        "rr\t0\tchr1\t10000\t0\t1H1S1I1M1I1S\t*\t0\t0\tACTGA\t?????\n"
        "rr\t0\tchr1\t10000\t0\t1H1S2M1I1S\t*\t0\t0\tACTGA\t?????\n"
        "rr\t0\tchr1\t10000\t0\t1H1S1M1I1M1S\t*\t0\t0\tACTGA\t?????\n"
        "rr\t0\tchr1\t10000\t0\t1H1S1M1D1M1S\t*\t0\t0\tACTG\t????\n"
        "rr\t0\tchr1\t10000\t0\t1H1S1I2M1S\t*\t0\t0\tACTGA\t?????\n"
        "rr\t0\tchr1\t10000\t0\t1H1S3M1S\t*\t0\t0\tACTGA\t?????\n"
        "rr\t0\tchr1\t10000\t0\t1H1S1M1I1S1H\t*\t0\t0\tACTG\t????\n"
        "rr\t0\tchr1\t10000\t0\t1H1S1I1M1S1H\t*\t0\t0\tACTG\t????\n"
        "rr\t0\tchr1\t10000\t0\t1H1S2M1S1H\t*\t0\t0\tACTG\t????\n"
        "rr\t0\tchr1\t10000\t0\t2M3I5M\t*\t0\t0\tACTGACTGAC\t??????????\n";
    // clang-format on

    std::vector<bam1_t *> bam_result;
    samFile *in = sam_open(sam, "r");
    bam_hdr_t *header = sam_hdr_read(in);
    bam1_t *aln;
    while (true) {
        aln = _bampool->alloc_bam_struct_pool();
        if (sam_read1(in, header, aln) >= 0) {
            _bampool->peek_bam(aln);
            bam_result.push_back(aln);
        }
        else {
            break;
        }
    }

    bam_hdr_destroy(header);
    sam_close(in);

    rovaca::pReadRecord read;
    std::vector<rovaca::pReadRecord> result;
    result.reserve(bam_result.size());
    for (bam1_t *b : bam_result) {
        read = ReadRecord::create(_pool, nullptr, b);
        result.push_back(read);
    }

    return result;
}

void ReadClipperUnitTest::assertRefAlignmentConsistent(rovaca::pReadRecord read)
{
    uint32_t cigarRefLength = bam_cigar2rlen((int32_t)read->cigar_length(), read->cigar());
    uint32_t readRefLength = read->is_unmapped() ? 0 : (uint32_t)read->get_length_on_reference();
    ASSERT_EQ(cigarRefLength, readRefLength);
}

void ReadClipperUnitTest::assertReadLengthConsistent(rovaca::pReadRecord read)
{
    uint32_t cigarRefLength = bam_cigar2qlen((int32_t)read->cigar_length(), read->cigar());
    uint32_t readRefLength = read->seq_length();
    ASSERT_EQ(cigarRefLength, readRefLength);
}