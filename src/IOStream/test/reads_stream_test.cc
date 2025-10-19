#include "reads_stream.h"

#include <gtest/gtest.h>

#include "bam_loader.h"
#include "bed_loader.h"
#include "downsampler_hc.h"
#include "fasta_loader.h"
#include "reads_filter_hc.h"
#include "reads_stream_builder.h"
#include "transformer.h"

std::string bam_path = "/data/rovaca-dev/sungeng/Datas/NA12878.sort.dup.bam";
std::string bed_path = "/data/rovaca-dev/wanghaoling/cWES/wes-annotation/src/db/PP100.gene.info.bed";
std::string ref_path = "/data/rovaca-dev/sungeng/Datas/hg19.fasta";

TEST(ReadsStreamTest, WES)
{
    contig_info_t fasta_dict;
    FastaLoader::get_fasta_dict(ref_path, &fasta_dict);
    BedLoader bd{bed_path, 0, fasta_dict};
    BamLoader bam_info{bam_path.c_str(), true};
    ReadFilter* reads_filter = new HCReadFilter(bam_info.get_sam_hdr()[0]);
    Transformer* reads_transformer = new Transformer();
    Downsampler* reads_downsampler = new HCDownsampler(50);
    ReadStreamBuilder builder(&bam_info);
    ReadStream reads = builder.set_downsampler(reads_downsampler).set_filter(reads_filter).set_transformer(reads_transformer);
    auto bed_intervals = bd.get_all_intervals();
    auto bed_key = bd.get_bed_keys();
    int32_t count = 0;
    reads = builder.set_target(bed_intervals);
    reads.init_reads_streamr();
    while (reads.has_next()) {
        auto item = reads.next();
        reads.read_recovery(item);
        count++;
    }
    delete reads_downsampler;
    delete reads_transformer;
    delete reads_filter;
}

TEST(ReadsStreamTest, WGSNODOWNSAMPLE)
{
    std::string line;
    htsThreadPool tpool{nullptr, 0};
    tpool.pool = hts_tpool_init(10);
    BamLoader bam_info{bam_path.c_str(), &tpool, false};
    ReadFilter* reads_filter = new HCReadFilter(bam_info.get_sam_hdr()[0]);
    Transformer* reads_transformer = new Transformer();
    ReadStreamBuilder builder(&bam_info);
    builder = builder.set_filter(reads_filter).set_transformer(reads_transformer);
    int32_t count = 0;
    ReadStream reads = builder.set_target("chrY");
    if (reads.init_reads_streamr() < 0) {
        std::cerr << "Init streamer error!" << std::endl;
        goto err;
    }
    while (reads.has_next()) {
        auto item = reads.next();
        reads.read_recovery(item);
        count++;
    }
    EXPECT_EQ(count, 280149);  // from gatk, raw reads 1035745 in total.
err:
    delete reads_transformer;
    delete reads_filter;
    hts_tpool_destroy(tpool.pool);
}

TEST(ReadsStreamTest, WGS)
{
    std::string line;
    htsThreadPool tpool{nullptr, 0};
    tpool.pool = hts_tpool_init(10);
    BamLoader bam_info{bam_path.c_str(), &tpool, false};
    ReadFilter* reads_filter = new HCReadFilter(bam_info.get_sam_hdr()[0]);
    Transformer* reads_transformer = new Transformer();
    Downsampler* reads_downsampler = new HCDownsampler(50);
    ReadStreamBuilder builder(&bam_info);
    builder = builder.set_downsampler(reads_downsampler).set_filter(reads_filter).set_transformer(reads_transformer);
    int32_t count = 0;
    ReadStream reads = builder.set_target("chrY");
    if (reads.init_reads_streamr() < 0) {
        std::cerr << "Init streamer error!" << std::endl;
        goto err;
    }
    while (reads.has_next()) {
        auto item = reads.next();
        reads.read_recovery(item);
        count++;
    }
    EXPECT_EQ(count, 266412);  // from gatk.
err:
    delete reads_downsampler;
    delete reads_transformer;
    delete reads_filter;
    hts_tpool_destroy(tpool.pool);
}