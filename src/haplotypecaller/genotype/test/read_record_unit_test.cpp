#include <gtest/gtest.h>

#include "assemble_interface.h"
#include "bam_data_pool.hpp"
#include "htslib/sam.h"
#include "read_record.h"

using namespace rovaca;
typedef std::pair<std::string, uint32_t> CigarPair;

static constexpr size_t s_buffer_size = 1024 * 1024 * 10;

class ReadRecordUnitTest : public ::testing::Test
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

    std::vector<bam1_t *> createBamData();

    uint8_t *_buffer{};
    pMemoryPool _pool{};
    pBamDataPool _bampool{};
};

static void commonTestbody(const pReadRecord &read, bam1_t *b)
{
    uint32_t i;

    // read     为 1 bases, 区间为 [start, stop]
    // bam1_t   为 0 bases, 区间为 [start, stop]
    ASSERT_EQ(read->get_tid(), b->core.tid);
    ASSERT_EQ(read->mate_tid(), b->core.mtid);
    ASSERT_EQ(read->get_start(), b->core.pos);
    ASSERT_EQ(read->mate_pos(), b->core.mpos);
    ASSERT_EQ(read->get_stop(), b->core.pos + bam_cigar2rlen(b->core.n_cigar, bam_get_cigar(b)) - 1);
    ASSERT_EQ(read->insert_size(), b->core.isize);
    ASSERT_EQ(read->cigar_length(), b->core.n_cigar);
    ASSERT_EQ(read->flag(), b->core.flag);
    ASSERT_EQ(read->mapping_quality(), b->core.qual);
    ASSERT_EQ(read->qname_len(), b->core.l_qname);
    ASSERT_EQ(read->seq_length(), b->core.l_qseq);
    ASSERT_EQ(read->get_sam_flags_for_read(), b->core.flag);

    ASSERT_FALSE(strncmp(read->qname(), bam_get_qname(b), b->core.l_qname));
    for (i = 0; i < b->core.n_cigar; ++i) {
        ASSERT_EQ(read->cigar_i(i), bam_get_cigar(b)[i]);
    }
    auto *seq = bam_get_seq(b);
    auto *qual = bam_get_qual(b);
    for (i = 0; i < (uint32_u)b->core.l_qseq; ++i) {
        ASSERT_EQ(read->seq_i(i), seq_nt16_str[bam_seqi(seq, i)]);
        ASSERT_EQ(read->qual_i(i), qual[i]);
    }

    // ASSERT_FALSE(read->is_paired());
    // ASSERT_FALSE(read->is_properly_paired());
    // ASSERT_FALSE(read->is_unmapped());
    // ASSERT_FALSE(read->mate_is_unmapped());
    // ASSERT_TRUE(read->is_reverse_strand());
    // ASSERT_FALSE(read->mate_is_reverse_strand());
    // ASSERT_FALSE(read->is_second_of_pair());
    // ASSERT_FALSE(read->is_first_of_pair());
    // ASSERT_FALSE(read->is_secondary_alignment());
    // ASSERT_FALSE(read->fails_vendor_quality_check());
    // ASSERT_FALSE(read->is_duplicate());
    // ASSERT_FALSE(read->is_supplementary_alignment());
}

TEST_F(ReadRecordUnitTest, testConstructByBam)
{
    auto data = createBamData();
    for (const auto b : data) {
        pReadRecord read = ReadRecord::create(_pool, nullptr, b);
        commonTestbody(read, b);
    }
}

TEST_F(ReadRecordUnitTest, testConstructByAssemble)
{
    auto data = createBamData();
    for (const auto b : data) {
        auto *assemble_read = new (_pool->allocate(sizeof(hc_apply_one_read))) hc_apply_one_read{};
        assemble_read->read = *b;
        assemble_read->pos_start = (int)b->core.pos;
        assemble_read->cigar_len = b->core.n_cigar;
        assemble_read->read_len = b->core.l_qseq;
        assemble_read->ref_len = bam_cigar2rlen(b->core.n_cigar, bam_get_cigar(b));
        assemble_read->insert_size = (int32_t)b->core.isize;
        assemble_read->cigar = (uint32_t *)bam_get_cigar(b);
        std::string seq(b->core.l_qseq, 'a');
        for (uint32_t i = 0; i < (uint32_t)b->core.l_qseq; ++i) {
            seq[i] = seq_nt16_str[bam_seqi(bam_get_seq(b), i)];
        }
        assemble_read->seq = (uint8_t *)seq.c_str();
        assemble_read->qual = bam_get_qual(b);
        auto read = ReadRecord::create(_pool, nullptr, assemble_read);
        commonTestbody(read, b);
    }
}

std::vector<bam1_t *> ReadRecordUnitTest::createBamData()
{
    // clang-format off
    const char sam[] = "data:,"
        "@SQ\tSN:chr1\tLN:249250621\n"
        "V350035846L1C001R0501038935\t129\tchr1\t10000\t0\t100M\t*\t0\t0\tATAACCCTAACTCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAA\t&????5????????????????????????????&+?+???????????????5????????????5???????????\'???????????5??????5??\n"
        "V350035846L2C006R0040057358\t401\tchr1\t11676\t0\t64S36M\t*\t0\t0\tTCAACACCGGCCATGCAGCAAAATCATCAGTGGAAATTTAAAAAAATACACATGGCCAGGCCCCAGCCCTGGAGATTCTTATTAGTGATTTGGGCTGGGG\t?????????????????????????????????????????????????????????????5??????????????????????????????????????\n"
        "V350035846L1C003R0310165749\t99\tchr1\t13362\t0\t56M4I40M\t*\t0\t230\tCTGCTGTGTGGAAGTTCACTCCTGCCTTTTCCTTTCCCTAGAGCCTCCACCACCCCGAGAGAGATCACATTTCTCACTGCCTTTTGTCTGCCCAGTTTCA\t??????????????????????????????????????????????????????????????????????????????5??????????????+??????\n"
        "V350035846L1C006R0221258126\t83\tchr1\t13635\t6\t22M2D78M\t*\t0\t-250\tATTAGTGCCCGTTGGAGAAAACGTGAATCGCGAAGAAATGGTAGGTCCTGGCCATCCGTGAGATCTTCCCAGGGCAACTCCACTCTGTGGAATCCAATCT\t???+5???$?5??????????5?+????5+#??????5???+&?+?????+?5????5?????????5+???\"?5?'5?55+??5?5?????????????\n";
    // clang-format on

    std::vector<bam1_t *> result;
    samFile *in = sam_open(sam, "r");
    bam_hdr_t *header = sam_hdr_read(in);
    bam1_t *aln;
    while (true) {
        aln = _bampool->alloc_bam_struct_pool();
        if (sam_read1(in, header, aln) >= 0) {
            _bampool->peek_bam(aln);
            result.push_back(aln);
        }
        else {
            break;
        }
    }

    bam_hdr_destroy(header);
    sam_close(in);

    return result;
}