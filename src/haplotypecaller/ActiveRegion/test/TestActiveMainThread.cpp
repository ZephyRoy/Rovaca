#include <thread>

#include "ActiveMainThread.h"
#include "bed_loader.h"
#include "fasta_loader.h"
#include "gtest/gtest.h"
#include "htslib/faidx.h"
#include "htslib/sam.h"

static constexpr size_t s_max_bufer_size = 1024 * 1024 * 20;

typedef std::pmr::memory_resource MemoryPool, *pMemoryPool;
class TestActiveMainThreadLogic : public ::testing::Test
{
protected:
    static void SetUpTestSuite()
    {
        m_fasta_loader = new ReferenceManager(std::string("/data/pipelines/WGS_bgionline/db/GRCh37/ref/GRCh37_no_alt.fna"), 5);
        fasta_info = new contig_info_t;
        FastaLoader::get_fasta_dict("/data/pipelines/WGS_bgionline/db/GRCh37/ref/GRCh37_no_alt.fna", fasta_info);
        thread_manager = new std::thread(&ReferenceManager::run, m_fasta_loader);
        _buffer = new uint8_t[s_max_bufer_size]{};
        _pool = new std::pmr::monotonic_buffer_resource(_buffer, s_max_bufer_size, std::pmr::null_memory_resource());
    }
    static void TearDownTestSuite()
    {
        m_fasta_loader->assign_finish();
        thread_manager->join();
        delete m_fasta_loader;
        delete fasta_info;
        delete thread_manager;
        delete _pool;
        delete[] _buffer;
    }
    void SetUp() {}
    void TearDown() {}
    static ReferenceManager *m_fasta_loader;
    static contig_info_t *fasta_info;
    static std::thread *thread_manager;
    static uint8_t *_buffer;
    static pMemoryPool _pool;
};
std::thread *TestActiveMainThreadLogic::thread_manager = nullptr;
ReferenceManager *TestActiveMainThreadLogic::m_fasta_loader = nullptr;
contig_info_t *TestActiveMainThreadLogic::fasta_info = nullptr;
uint8_t *TestActiveMainThreadLogic::_buffer = nullptr;
pMemoryPool TestActiveMainThreadLogic::_pool = nullptr;
TEST_F(TestActiveMainThreadLogic, bedtest1)
{
    BedLoader *bed_loader =
        new BedLoader(std::string("/home/yinlonghui/workspace/hc/src/haplotypecaller/ActiveRegion/test/Resource/a.bed"), 0, *fasta_info);
    ActiveMainThreadDispatchTasks inst(nullptr, bed_loader, m_fasta_loader, fasta_info, nullptr, nullptr, nullptr, nullptr, nullptr);

    boost::dynamic_bitset<> *_bit = new boost::dynamic_bitset<>(1024, 0x0);
    int max = 1024;

    inst.get_target_biset(0, 0, 20, _bit, max);

    for (int i = 0; i < 2; i++) {
        EXPECT_FALSE(_bit->test(i));
    }

    for (int i = 2; i < 3; i++) {
        EXPECT_TRUE(_bit->test(i));
    }

    for (int i = 3; i < 10; i++) {
        EXPECT_FALSE(_bit->test(i));
    }

    for (int i = 10; i < 11; i++) {
        EXPECT_TRUE(_bit->test(i));
    }

    for (int i = 11; i < 20; i++) {
        EXPECT_FALSE(_bit->test(i));
    }

    inst.get_target_biset(1, 0, 20, _bit, max);
    for (int i = 0; i < 20; i++) {
        EXPECT_FALSE(_bit->test(i));
    }
    inst.get_target_biset(2, 0, 20, _bit, max);

    for (int i = 0; i < 10; i++) {
        EXPECT_FALSE(_bit->test(i));
    }

    for (int i = 10; i < 20; i++) {
        EXPECT_TRUE(_bit->test(i));
    }

    delete _bit;
    delete bed_loader;
}

TEST_F(TestActiveMainThreadLogic, split_long_reads_case1)
{
    ActiveRegionBamBlockListSource resource(10);
    ActiveMainThreadDispatchTasks inst(nullptr, nullptr, m_fasta_loader, fasta_info, &resource, nullptr, nullptr, nullptr, nullptr);
    const char sam[] =
        "data:,"
        "@SQ\tSN:chr1\tLN:249250621\n"
        "@SQ\tSN:chr2\tLN:243199373\n"
        "1\t129\tchr1\t10000\t0\t10S10M\t*\t0\t0\tAAAAAAAAAAAAAAAAAAAA\tAAAAAAAAAAAAAAAAAAAA\n"
        "1\t129\tchr1\t10010\t0\t10S5M5S\t*\t0\t0\tAAAAAAAAAAAAAAAAAAAA\tAAAAAAAAAAAAAAAAAAAA\n"
        "1\t129\tchr2\t10000\t0\t10S5M5S\t*\t0\t0\tAAAAAAAAAAAAAAAAAAAA\t!AAAAAAAAAAAAAAAAAAA\n"
        "1\t129\tchr2\t10020\t0\t10S5M5S\t*\t0\t0\tAAAAAAAAAAAAAAAAAAAA\t!AAAAAAAAAAAAA!AAAAA\n";

    samFile *in = sam_open(sam, "r");
    bam_hdr_t *h = sam_hdr_read(in);
    bam1_t *bam = bam_init1();
    for (int i = 0; i < 4; i++) {
        EXPECT_TRUE(sam_read1(in, h, bam) >= 0);
        resource.insert(bam);
    }
    resource.m_reads_buffer.back().is_out_interval = true;
    resource.m_reads_buffer.back().is_out_chrom = true;
    bam_destroy1(bam);
    EXPECT_EQ(resource.m_reads_buffer.size(), 2);

    EXPECT_EQ(resource.m_reads_buffer.front().tid, 0);
    EXPECT_EQ(resource.m_reads_buffer.front().start, 9999);
    EXPECT_EQ(resource.m_reads_buffer.front().end, 10009);

    EXPECT_EQ(resource.m_reads_buffer.back().tid, 1);
    EXPECT_EQ(resource.m_reads_buffer.back().start, 9999);
    EXPECT_EQ(resource.m_reads_buffer.back().end, 10019);
    std::pmr::vector<bam1_t *> reads(_pool);
    int current_tid = -1;
    hts_pos_t start = -1, end = -1;
    hts_pos_t actual_start = -1, actual_end = -1;
    inst.split_reads(current_tid, start, end, actual_start, actual_end, reads);
    EXPECT_EQ(current_tid, 0);
    EXPECT_EQ(start, 0);
    EXPECT_EQ(end, 249250621);
    EXPECT_EQ(current_tid, 0);
    EXPECT_EQ(actual_start, 9999);
    EXPECT_EQ(actual_end, 10014);
    EXPECT_EQ(reads.size(), 2);
    reads.clear();
    inst.split_reads(current_tid, start, end, actual_start, actual_end, reads);
    EXPECT_EQ(current_tid, 1);
    EXPECT_EQ(start, 0);
    EXPECT_EQ(end, 243199373);
    EXPECT_EQ(actual_start, 9999);
    EXPECT_EQ(actual_end, 10024);
    EXPECT_EQ(reads.size(), 2);
    reads.clear();
    for (int i = 1; i < 84; i++) {
        EXPECT_TRUE(inst.split_reads(current_tid, start, end, actual_start, actual_end, reads)) << "TID:" << current_tid << std::endl;
        EXPECT_EQ(current_tid, i + 1);
    }
    EXPECT_FALSE(inst.split_reads(current_tid, start, end, actual_start, actual_end, reads)) << "TID:" << current_tid << std::endl;
}

TEST_F(TestActiveMainThreadLogic, split_long_reads_case2)
{
    ActiveRegionBamBlockListSource resource(10);
    ActiveMainThreadDispatchTasks inst(nullptr, nullptr, m_fasta_loader, fasta_info, &resource, nullptr, nullptr, nullptr, nullptr);
    const char sam[] =
        "data:,"
        "@SQ\tSN:chr1\tLN:249250621\n"
        "@SQ\tSN:chr2\tLN:243199373\n"
        "@SQ\tSN:chr3\tLN:198022430\n"
        "1\t129\tchr1\t10000\t0\t10S30M\t*\t0\t0\tAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA\tAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA\n"
        "2\t129\tchr1\t10010\t0\t10S5M5S\t*\t0\t0\tAAAAAAAAAAAAAAAAAAAA\tAAAAAAAAAAAAAAAAAAAA\n"
        "3\t129\tchr3\t10000\t0\t10S5M5S\t*\t0\t0\tAAAAAAAAAAAAAAAAAAAA\t!AAAAAAAAAAAAAAAAAAA\n"
        "4\t129\tchr3\t11001\t0\t10S5M5S\t*\t0\t0\tAAAAAAAAAAAAAAAAAAAA\t!AAAAAAAAAAAAA!AAAAA\n";

    samFile *in = sam_open(sam, "r");
    bam_hdr_t *h = sam_hdr_read(in);
    bam1_t *bam = bam_init1();
    for (int i = 0; i < 4; i++) {
        EXPECT_TRUE(sam_read1(in, h, bam) >= 0);
        resource.insert(bam);
    }
    resource.m_reads_buffer.back().is_out_interval = true;
    resource.m_reads_buffer.back().is_out_chrom = true;
    bam_destroy1(bam);
    std::pmr::vector<bam1_t *> reads(_pool);
    int current_tid = -1;
    hts_pos_t start = -1, end = -1;
    hts_pos_t actual_start = -1, actual_end = -1;
    inst.split_reads(current_tid, start, end, actual_start, actual_end, reads);
    EXPECT_EQ(current_tid, 0);
    EXPECT_EQ(start, 0);
    EXPECT_EQ(end, 249250621);
    EXPECT_EQ(current_tid, 0);
    EXPECT_EQ(actual_start, 9999);
    EXPECT_EQ(actual_end, 10029);
    EXPECT_EQ(reads.size(), 2);

    reads.clear();
    inst.split_reads(current_tid, start, end, actual_start, actual_end, reads);
    EXPECT_EQ(current_tid, 1);
    EXPECT_EQ(start, 0);
    EXPECT_EQ(end, 243199373);
    EXPECT_EQ(actual_start, 0);
    EXPECT_EQ(actual_end, 0);

    EXPECT_EQ(reads.size(), 0);

    reads.clear();
    inst.split_reads(current_tid, start, end, actual_start, actual_end, reads);
    EXPECT_EQ(current_tid, 2);
    EXPECT_EQ(start, 0);
    EXPECT_EQ(end, 10004);
    EXPECT_EQ(actual_start, 9999);
    EXPECT_EQ(actual_end, 10004);
    EXPECT_EQ(reads.size(), 1);
    reads.clear();
    inst.split_reads(current_tid, start, end, actual_start, actual_end, reads);
    EXPECT_EQ(current_tid, 2);
    EXPECT_EQ(start, 10004);
    EXPECT_EQ(end, 198022430);
    EXPECT_EQ(actual_start, 11000);
    EXPECT_EQ(actual_end, 11005);
    EXPECT_EQ(reads.size(), 1);
}

TEST_F(TestActiveMainThreadLogic, split_long_reads_case3)
{
    ActiveRegionBamBlockListSource resource(10);
    ActiveMainThreadDispatchTasks inst(nullptr, nullptr, m_fasta_loader, fasta_info, &resource, nullptr, nullptr, nullptr, nullptr);
    const char sam[] =
        "data:,"
        "@SQ\tSN:chr1\tLN:249250621\n"
        "@SQ\tSN:chr2\tLN:243199373\n"
        "@SQ\tSN:chr3\tLN:198022430\n"
        "1\t129\tchr2\t10000\t0\t10S10M\t*\t0\t0\tAAAAAAAAAAAAAAAAAAAA\tAAAAAAAAAAAAAAAAAAAA\n"
        "2\t129\tchr2\t10001\t0\t10S5M5S\t*\t0\t0\tAAAAAAAAAAAAAAAAAAAA\tAAAAAAAAAAAAAAAAAAAA\n"
        "3\t129\tchr2\t10002\t0\t10S5M5S\t*\t0\t0\tAAAAAAAAAAAAAAAAAAAA\t!AAAAAAAAAAAAAAAAAAA\n"
        "4\t129\tchr2\t10003\t0\t10S5M5S\t*\t0\t0\tAAAAAAAAAAAAAAAAAAAA\t!AAAAAAAAAAAAA!AAAAA\n"
        "5\t129\tchr2\t10004\t0\t10S5M5S\t*\t0\t0\tAAAAAAAAAAAAAAAAAAAA\t!AAAAAAAAAAAAAAAAAAA\n"
        "6\t129\tchr2\t10005\t0\t10S5M5S\t*\t0\t0\tAAAAAAAAAAAAAAAAAAAA\t!AAAAAAAAAAAAAAAAAAA\n"
        "7\t129\tchr2\t10006\t0\t10S5M5S\t*\t0\t0\tAAAAAAAAAAAAAAAAAAAA\t!AAAAAAAAAAAAAAAAAAA\n"
        "8\t129\tchr2\t10007\t0\t10S5M5S\t*\t0\t0\tAAAAAAAAAAAAAAAAAAAA\t!AAAAAAAAAAAAAAAAAAA\n"
        "9\t129\tchr2\t10008\t0\t10S5M5S\t*\t0\t0\tAAAAAAAAAAAAAAAAAAAA\t!AAAAAAAAAAAAAAAAAAA\n"
        "10\t129\tchr2\t10009\t0\t10S5M5S\t*\t0\t0\tAAAAAAAAAAAAAAAAAAAA\t!AAAAAAAAAAAAAAAAAAA\n"
        "11\t129\tchr2\t10010\t0\t10S5M5S\t*\t0\t0\tAAAAAAAAAAAAAAAAAAAA\t!AAAAAAAAAAAAAAAAAAA\n";

    samFile *in = sam_open(sam, "r");
    bam_hdr_t *h = sam_hdr_read(in);
    bam1_t *bam = bam_init1();
    for (int i = 0; i < 4; i++) {
        EXPECT_TRUE(sam_read1(in, h, bam) >= 0);
        resource.insert(bam);
    }
    resource.m_finalize = true;
    for (int i = 0; i < 4; i++) {
        EXPECT_TRUE(sam_read1(in, h, bam) >= 0);
        resource.insert(bam);
    }
    resource.m_finalize = true;
    for (int i = 0; i < 3; i++) {
        EXPECT_TRUE(sam_read1(in, h, bam) >= 0);
        resource.insert(bam);
    }
    resource.m_finalize = true;

    resource.m_reads_buffer.back().is_out_interval = true;
    resource.m_reads_buffer.back().is_out_chrom = true;

    bam_destroy1(bam);
    std::pmr::vector<bam1_t *> reads(_pool);
    int current_tid = -1;
    hts_pos_t start = -1, end = -1;
    hts_pos_t actual_start = -1, actual_end = -1;
    inst.split_reads(current_tid, start, end, actual_start, actual_end, reads);
    EXPECT_EQ(current_tid, 0);
    EXPECT_EQ(start, 0);
    EXPECT_EQ(end, 249250621);
    EXPECT_EQ(actual_start, 0);
    EXPECT_EQ(actual_end, 0);
    EXPECT_EQ(reads.size(), 0);

    inst.split_reads(current_tid, start, end, actual_start, actual_end, reads);
    EXPECT_EQ(current_tid, 1);
    EXPECT_EQ(start, 0);
    EXPECT_EQ(end, 10003);
    EXPECT_EQ(actual_start, 9999);
    EXPECT_EQ(actual_end, 10003);
    EXPECT_EQ(reads.size(), 4);
    reads.clear();
    inst.split_reads(current_tid, start, end, actual_start, actual_end, reads);

    EXPECT_EQ(current_tid, 1);
    EXPECT_EQ(start, 10003);
    EXPECT_EQ(end, 10007);
    EXPECT_EQ(reads.size(), 8);
    reads.clear();
    inst.split_reads(current_tid, start, end, actual_start, actual_end, reads);
    EXPECT_EQ(current_tid, 1);
    EXPECT_EQ(start, 10007);
    EXPECT_EQ(end, 243199373);
    EXPECT_EQ(reads.size(), 8);
}

TEST_F(TestActiveMainThreadLogic, get_region_reads_block_case1)
{
    ActiveRegionBamBlockListSource resource(10);
    rovaca::BamDataPool pool(1024);

    ActiveMainThreadReduce inst(nullptr, nullptr, nullptr, nullptr, nullptr, &resource, nullptr, nullptr, true);
    const char sam[] =
        "data:,"
        "@SQ\tSN:chr1\tLN:249250621\n"
        "@SQ\tSN:chr2\tLN:243199373\n"
        "@SQ\tSN:chr3\tLN:198022430\n"
        "1\t129\tchr1\t10000\t0\t10S10M\t*\t0\t0\tAAAAAAAAAAAAAAAAAAAA\tAAAAAAAAAAAAAAAAAAAA\n"
        "2\t129\tchr1\t10100\t0\t10S5M5S\t*\t0\t0\tAAAAAAAAAAAAAAAAAAAA\tAAAAAAAAAAAAAAAAAAAA\n"
        "3\t129\tchr1\t10200\t0\t10S5M5S\t*\t0\t0\tAAAAAAAAAAAAAAAAAAAA\t!AAAAAAAAAAAAAAAAAAA\n"
        "4\t129\tchr1\t10300\t0\t10S5M5S\t*\t0\t0\tAAAAAAAAAAAAAAAAAAAA\t!AAAAAAAAAAAAA!AAAAA\n"
        "5\t129\tchr1\t10400\t0\t10S5M5S\t*\t0\t0\tAAAAAAAAAAAAAAAAAAAA\t!AAAAAAAAAAAAAAAAAAA\n"
        "6\t129\tchr1\t10500\t0\t10S5M5S\t*\t0\t0\tAAAAAAAAAAAAAAAAAAAA\t!AAAAAAAAAAAAAAAAAAA\n"
        "7\t129\tchr1\t10600\t0\t10S5M5S\t*\t0\t0\tAAAAAAAAAAAAAAAAAAAA\t!AAAAAAAAAAAAAAAAAAA\n"
        "8\t129\tchr2\t10700\t0\t10S5M5S\t*\t0\t0\tAAAAAAAAAAAAAAAAAAAA\t!AAAAAAAAAAAAAAAAAAA\n"
        "9\t129\tchr2\t10800\t0\t10S5M5S\t*\t0\t0\tAAAAAAAAAAAAAAAAAAAA\t!AAAAAAAAAAAAAAAAAAA\n"
        "10\t129\tchr2\t10900\t0\t10S5M5S\t*\t0\t0\tAAAAAAAAAAAAAAAAAAAA\t!AAAAAAAAAAAAAAAAAAA\n";

    samFile *in = sam_open(sam, "r");
    bam_hdr_t *h = sam_hdr_read(in);

    std::pmr::vector<bam1_t *> reads(_pool);
    p_hc_region_active_storage region = new hc_region_active_storage;

    region->tid = 0;
    region->start_index = 10109;
    region->end_index = 10400;
    bam1_t *bam = bam_init1();
    for (int i = 0; i < 10; i++) {
        EXPECT_TRUE(sam_read1(in, h, bam) >= 0);
        resource.insert(bam);
    }
    inst.get_region_reads_block(&pool, region, reads);
    EXPECT_EQ(reads.size(), 5);
    reads.clear();
    region->tid = 0;
    region->start_index = 10109;
    region->end_index = 10109;

    inst.get_region_reads_block(&pool, region, reads);
    EXPECT_EQ(reads.size(), 2);
    reads.clear();
    region->tid = 1;
    region->start_index = 10800;
    region->end_index = 10800;
    inst.get_region_reads_block(&pool, region, reads);
    EXPECT_EQ(reads.size(), 3);
    bam_destroy1(bam);
    bam_hdr_destroy(h);
    sam_close(in);
}

TEST_F(TestActiveMainThreadLogic, remove_non_overlapping_reads_region_long_case1)
{
    ActiveRegionBamBlockListSource resource(10);
    rovaca::BamDataPool pool(1024);

    ActiveMainThreadReduce inst(nullptr, nullptr, nullptr, nullptr, nullptr, &resource, nullptr, nullptr, true);
    const char sam[] =
        "data:,"
        "@SQ\tSN:chr1\tLN:249250621\n"
        "@SQ\tSN:chr2\tLN:243199373\n"
        "@SQ\tSN:chr3\tLN:198022430\n"
        "1\t129\tchr1\t10000\t0\t10S10M\t*\t0\t0\tAAAAAAAAAAAAAAAAAAAA\tAAAAAAAAAAAAAAAAAAAA\n"
        "2\t129\tchr1\t10003\t0\t10S5M5S\t*\t0\t0\tAAAAAAAAAAAAAAAAAAAA\tAAAAAAAAAAAAAAAAAAAA\n"
        "3\t129\tchr1\t10200\t0\t10S5M5S\t*\t0\t0\tAAAAAAAAAAAAAAAAAAAA\t!AAAAAAAAAAAAAAAAAAA\n"
        "4\t129\tchr1\t10300\t0\t10S5M5S\t*\t0\t0\tAAAAAAAAAAAAAAAAAAAA\t!AAAAAAAAAAAAA!AAAAA\n"
        "5\t129\tchr1\t10400\t0\t10S5M5S\t*\t0\t0\tAAAAAAAAAAAAAAAAAAAA\t!AAAAAAAAAAAAAAAAAAA\n"
        "6\t129\tchr1\t10500\t0\t10S5M5S\t*\t0\t0\tAAAAAAAAAAAAAAAAAAAA\t!AAAAAAAAAAAAAAAAAAA\n"
        "7\t129\tchr1\t10600\t0\t10S5M5S\t*\t0\t0\tAAAAAAAAAAAAAAAAAAAA\t!AAAAAAAAAAAAAAAAAAA\n"
        "8\t129\tchr2\t10700\t0\t10S5M5S\t*\t0\t0\tAAAAAAAAAAAAAAAAAAAA\t!AAAAAAAAAAAAAAAAAAA\n"
        "9\t129\tchr2\t10800\t0\t10S5M5S\t*\t0\t0\tAAAAAAAAAAAAAAAAAAAA\t!AAAAAAAAAAAAAAAAAAA\n"
        "10\t129\tchr2\t10900\t0\t10S5M5S\t*\t0\t0\tAAAAAAAAAAAAAAAAAAAA\t!AAAAAAAAAAAAAAAAAAA\n";

    samFile *in = sam_open(sam, "r");
    bam_hdr_t *h = sam_hdr_read(in);
    bam1_t *bam = bam_init1();

    for (int i = 0; i < 2; i++) {
        EXPECT_TRUE(sam_read1(in, h, bam) >= 0);
        resource.insert(bam);
    }
    resource.m_finalize = true;
    for (int i = 0; i < 4; i++) {
        EXPECT_TRUE(sam_read1(in, h, bam) >= 0);
        resource.insert(bam);
    }
    resource.m_finalize = true;
    for (int i = 0; i < 4; i++) {
        EXPECT_TRUE(sam_read1(in, h, bam) >= 0);
        resource.insert(bam);
    }
    EXPECT_EQ(resource.m_reads_buffer.size(), 4);
    EXPECT_EQ(resource.m_reads_idle.size(), 0);
    p_hc_region_active_storage region = new hc_region_active_storage;
    region->tid = 0;
    region->start_index = 10102;
    region->end_index = 10109;
    inst.remove_non_overlapping_reads_region_block(region);
    EXPECT_EQ(resource.m_reads_buffer.size(), 4);
    EXPECT_EQ(resource.m_reads_idle.size(), 0);

    region->tid = 0;
    region->start_index = 10108;
    region->end_index = 10210;
    inst.remove_non_overlapping_reads_region_block(region);
    EXPECT_EQ(resource.m_reads_buffer.size(), 4);
    EXPECT_EQ(resource.m_reads_idle.size(), 0);

    region->tid = 0;
    region->start_index = 10209;
    region->end_index = 10210;
    inst.remove_non_overlapping_reads_region_block(region);
    EXPECT_EQ(resource.m_reads_buffer.size(), 3);
    EXPECT_EQ(resource.m_reads_idle.size(), 2);

    region->tid = 0;
    region->start_index = 10500;
    region->end_index = 10600;
    inst.remove_non_overlapping_reads_region_block(region);
    EXPECT_EQ(resource.m_reads_buffer.size(), 3);
    EXPECT_EQ(resource.m_reads_idle.size(), 2);

    region->tid = 0;
    region->start_index = 10700;
    region->end_index = 10800;
    inst.remove_non_overlapping_reads_region_block(region);
    EXPECT_EQ(resource.m_reads_buffer.size(), 2);
    EXPECT_EQ(resource.m_reads_idle.size(), 6);

    region->tid = 1;
    region->start_index = 0;
    region->end_index = 200;
    inst.remove_non_overlapping_reads_region_block(region);
    EXPECT_EQ(resource.m_reads_buffer.size(), 1);
    EXPECT_EQ(resource.m_reads_idle.size(), 7);

    region->tid = 1;
    region->start_index = 10800;
    region->end_index = 10805;
    inst.remove_non_overlapping_reads_region_block(region);

    EXPECT_EQ(resource.m_reads_buffer.size(), 1);
    EXPECT_EQ(resource.m_reads_idle.size(), 7);

    region->tid = 1;
    region->start_index = 110086;
    region->end_index = 110086;
    inst.remove_non_overlapping_reads_region_block(region);
    EXPECT_EQ(resource.m_reads_buffer.size(), 0);
    EXPECT_EQ(resource.m_reads_idle.size(), 10);
    bam_hdr_destroy(h);
    sam_close(in);
    bam_destroy1(bam);
}