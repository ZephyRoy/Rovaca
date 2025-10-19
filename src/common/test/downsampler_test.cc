#include <gtest/gtest.h>

#include <cstdlib>
#include <ctime>
#include <string>
#include <vector>

#include "downsampler_hc.h"
#include "htslib/sam.h"

// Create a mapped read with the specified length, reference ID, and position
bam1_t* CreateOneMappedRead(int32_t ref_id, int position)
{
    bam1_t* read = bam_init1();
    read->core.tid = ref_id;
    read->core.pos = position;
    read->core.qual = 30;
    read->core.l_qname = 5;
    read->core.n_cigar = 1;
    read->core.l_qseq = 5;
    read->core.flag = 0;
    read->core.mtid = -1;
    read->core.mpos = -1;
    read->core.isize = 0;
    uint8_t* cigar_data = new uint8_t[4];
    cigar_data[0] = 0x05;
    cigar_data[1] = BAM_CMATCH;
    cigar_data[2] = 0x00;
    cigar_data[3] = 0x00;
    read->data = cigar_data;
    return read;
}

// Create a unmapped read with the specified length
bam1_t* CreateOneUnmappedRead()
{
    bam1_t* read = bam_init1();
    read->core.tid = -1;
    read->core.pos = -1;
    read->core.qual = 30;
    read->core.l_qname = 5;
    read->core.n_cigar = 1;
    read->core.l_qseq = 5;
    read->core.flag = BAM_FUNMAP;
    read->core.mtid = -1;
    read->core.mpos = -1;
    read->core.isize = 0;
    uint8_t* cigar_data = new uint8_t[4];
    cigar_data[0] = 0x05;
    cigar_data[1] = BAM_CMATCH;
    cigar_data[2] = 0x00;
    cigar_data[3] = 0x00;
    read->data = cigar_data;
    return read;
}

// Create a unmapped read with the specified length and position
bam1_t* CreateOneUnmappedReadsWithPosition(int32_t ref_id, int position)
{
    bam1_t* read = bam_init1();
    read->core.tid = ref_id;
    read->core.pos = position;
    read->core.qual = 30;
    read->core.l_qname = 5;
    read->core.n_cigar = 1;
    read->core.l_qseq = 5;
    read->core.flag = BAM_FUNMAP;
    read->core.mtid = -1;
    read->core.mpos = -1;
    read->core.isize = 0;
    uint8_t* cigar_data = new uint8_t[4];
    cigar_data[0] = 0x05;
    cigar_data[1] = BAM_CMATCH;
    cigar_data[2] = 0x00;
    cigar_data[3] = 0x00;
    read->data = cigar_data;
    return read;
}

class HcDownsamperTest : public ::testing::Test
{
protected:
    void SetUp() override
    {
        downsampler = new HCDownsampler(100);
        for (int i = 0; i < 10; i++) {
            mapped_reads.emplace_back(CreateOneMappedRead(0, 100));
        }
    }

    void TearDown() override { delete downsampler; }

public:
    HCDownsampler* downsampler;
    std::vector<bam1_t*> mapped_reads;
    std::vector<bam1_t*> unmapped_reads;
};

TEST_F(HcDownsamperTest, TestSubmit)
{
    bam1_t* item = CreateOneMappedRead(0, 100);
    downsampler->submit(item);
    EXPECT_EQ(downsampler->has_finalized_items(), false);
    bam_destroy1(item);
}

TEST_F(HcDownsamperTest, TestHasFinalizedItems) { EXPECT_EQ(downsampler->has_finalized_items(), false); }

TEST_F(HcDownsamperTest, TestSignalEndOfInput)
{
    bam1_t* item = CreateOneMappedRead(0, 100);
    downsampler->submit(item);
    downsampler->input_end_signal();
    EXPECT_EQ(downsampler->has_finalized_items(), true);
    bam_destroy1(item);
}

TEST_F(HcDownsamperTest, TestConsumeFinalizedItems)
{
    bam1_t* item1 = CreateOneMappedRead(0, 100);
    downsampler->submit(item1);

    bam1_t* item2 = CreateOneMappedRead(0, 200);
    downsampler->submit(item2);
    downsampler->input_end_signal();
    std::vector<bam1_t*> finalized_items;
    downsampler->consume_finalized_items(finalized_items);
    EXPECT_EQ(finalized_items.size(), 2);
    EXPECT_EQ(finalized_items[0]->core.pos, 100);
    EXPECT_EQ(finalized_items[1]->core.pos, 200);
    for (auto item : finalized_items) {
        bam_destroy1(item);
    }
}

TEST_F(HcDownsamperTest, TestHandlePositionalChange)
{
    bam1_t* item1 = CreateOneMappedRead(0, 100);
    downsampler->submit(item1);

    bam1_t* item2 = CreateOneMappedRead(0, 200);
    downsampler->submit(item2);
    std::vector<bam1_t*> finalized_items;
    downsampler->consume_finalized_items(finalized_items);
    EXPECT_EQ(finalized_items.size(), 1);
    for (auto item : finalized_items) {
        bam_destroy1(item);
    }
}