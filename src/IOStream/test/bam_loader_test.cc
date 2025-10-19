#include "../bam_loader.h"

#include <gtest/gtest.h>

#include <chrono>
#include <string>

#include "../fasta_loader.h"
#include "htslib/thread_pool.h"

void print_bam_record(const bam1_t* bam_record)
{
    int32_t tid = bam_record->core.tid;
    std::cout << bam_get_qname(bam_record) << "\t" << bam_record->core.flag << "\t" << tid << "\t" << bam_record->core.pos + 1 << "\n";
}

std::string path1 = "/data/rovaca-dev/sungeng/Datas/NA12878.sort.dup.bam";
std::string path3 = "/data/rovaca-dev/sungeng/Datas/mutect2/normal.bam";
std::string path4 = "/data/rovaca-dev/sungeng/Datas/mutect2/tumor.bam";
std::string path5 = "/data/rovaca-dev/wanghaoling/workdir/work/21S18093206/bwa/21S18093206.final.bam";
std::string path6 = "/data/rovaca-dev/wanghaoling/workdir/work/21B20366817/bwa/21B20366817.final.bam";
std::string path7 = "/data/rovaca-dev/wanghaoling/workdir/work/21B20366816/bwa/21B20366816.final.bam";
std::string path8 = "/data/rovaca-dev/wanghaoling/workdir/work/21S17435345/bwa/21S17435345.final.bam";

TEST(BamLoaderTest, LoadSingleBamFile)
{
    htsThreadPool tpool{nullptr, 0};
    tpool.pool = hts_tpool_init(10);
    BamLoader loader(path1.c_str(), &tpool, false);
    bool success = loader.has_next();
    EXPECT_TRUE(success);
}

TEST(BamLoaderTest, LoadMultipleBamFiles)
{
    htsThreadPool tpool{nullptr, 0};
    tpool.pool = hts_tpool_init(10);
    std::vector<const char*> bam_files = {path3.c_str(), path4.c_str()};
    BamLoader* loader = new BamLoader(bam_files, &tpool, false);
    bool success = loader->has_next();
    EXPECT_TRUE(success);
    delete loader;
}

TEST(BamLoaderTest, GetNextReadFromSingle)
{
    htsThreadPool tpool{nullptr, 0};
    tpool.pool = hts_tpool_init(10);
    BamLoader* loader = new BamLoader(path1.c_str(), &tpool, true);
    const char* region = "chr4:15597707-15597830";
    const char* regarray[] = {region};

    bed_intervals region_array{regarry : const_cast<char**>(regarray), start : nullptr, end : nullptr, m : 10, n : 1};
    loader->set_target(&region_array);
    while (loader->has_next()) {
        int file_index;
        auto read = loader->get_next_read(file_index);
        if (!read) continue;
        loader->read_recovery(read);
        file_index++;
    }
    delete loader;
}

TEST(BamLoaderTest, GetNextReadFromSingleFileMultiRegion)
{
    htsThreadPool tpool{nullptr, 0};
    tpool.pool = hts_tpool_init(10);
    BamLoader loader(path1.c_str(), &tpool, true);
    const char* region1 = "chr1:24126876-24127034";
    const char* region2 = "chr1:169509531-169512352";
    const char* region3 = "chr2:152422229-152422334";
    const char* region4 = "chr3:48621730-48621808";
    const char* region5 = "chr4:15597707-15597830";
    const char* regarray[] = {region1, region2, region3, region4, region5};

    bed_intervals jj{regarry : const_cast<char**>(regarray), start : nullptr, end : nullptr, m : 10, n : 5};
    loader.set_target(&jj);
    int i = 0;
    while (loader.has_next()) {
        int file_index;
        auto read = loader.get_next_read(file_index);
        if (!read) continue;
        loader.read_recovery(read);
        i++;
    }
    EXPECT_EQ(i, 873);
}

TEST(BamLoaderTest, GetNextReadFromMultiFileMultiRegion)
{
    htsThreadPool tpool{nullptr, 0};
    tpool.pool = hts_tpool_init(10);
    std::vector<const char*> bam_files = {path5.c_str(), path6.c_str()};
    BamLoader loader(bam_files, &tpool, true);
    const char* region1 = "chr1:24126876-24127034";
    const char* region2 = "chr1:169509531-169512352";
    const char* region3 = "chr2:152422229-152422334";
    const char* region4 = "chr3:48621730-48621808";
    const char* region5 = "chr4:15597707-15597830";
    const char* regarray[] = {region1, region2, region3, region4, region5};
    bed_intervals jj{regarry : const_cast<char**>(regarray), start : nullptr, end : nullptr, m : 10, n : 5};
    loader.set_target(&jj);

    int i = 1;
    while (loader.has_next()) {
        int file_index;
        auto read = loader.get_next_read(file_index);
        if (!read) continue;
        loader.read_recovery(read);
        i++;
    }
    EXPECT_EQ(i, 28815);  // 未排序结果为28814。
}

TEST(BamLoaderTest, 30x_WGS_64G)
{
    htsThreadPool tpool{nullptr, 0};
    tpool.pool = hts_tpool_init(10);
    BamLoader loader(path1.c_str(), &tpool, false);
    int64_t count = 0;
    loader.set_target("chrY");
    while (loader.has_next()) {
        int file_index;
        auto read = loader.get_next_read(file_index);
        if (!read) continue;
        loader.read_recovery(read);
        count++;
    }
    EXPECT_EQ(count, 1035745);
}