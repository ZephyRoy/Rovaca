#include <gtest/gtest.h>

#include <string>

#include "htslib/sam.h"
#include "reads_filter_hc.h"
#include "reads_filter_lib.h"

TEST(ReadFilterTest, WellFormed)
{
    const char sam[] =
        "data:,"
        "@SQ\tSN:chr1\tLN:249250621\n"
        "read0\t0\tchr1\t10\t0\t10M\t*\t0\t0\tTACGTACGTA\t!@#$%^&*()\tRG:Z:NA12878\n"
        "read1\t0\tchr1\t80\t60\t10M\t*\t0\t0\tTACGTACGTA\t!@#$%^&*()\tRG:Z:NA12878\n"
        "read2\t0\tchr1\t10\t60\t10S\t*\t0\t0\tGTACGTACGT\t!@#$%^&*()\tRG:Z:NA12878\n"
        "read3\t0\tchr1\t30\t60\t10M\t*\t0\t0\tACGTACGTAC\t!@#$%^&*()\tRG:Z:NA12878\n"
        "read4\t0\tchr1\t56\t60\t10M\t*\t0\t0\tGTACGTACGT\t!@#$%^&*()\tRG:Z:NA12878\n"
        "read5\t0\tchr1\t70\t60\t10M\t*\t0\t0\tTACGTACGTA\t!@#$%^&*()\tRG:Z:NA12878\n"
        "read6\t0\tchr1\t40\t60\t10M\t*\t0\t0\tACGTACGTAC\t!@#$%^&*()\tRG:Z:NA12878\n"
        "read7\t0\tchr1\t20\t60\t10M\t*\t0\t0\tGTACGTACGT\t!@#$%^&*()\tRG:Z:NA12878\n"
        "read8\t0\tchr1\t60\t60\t10M\t*\t0\t0\tTACGTACGTA\t!@#$%^&*()\tRG:Z:NA12878\n"
        "read9\t0\tchr1\t5\t60\t10M\t*\t0\t0\tACGTACGTAC\t!@#$%^&*()\tRG:Z:NA12878\n";
    samFile* fp = sam_open(sam, "r");
    bam_hdr_t* hdr = sam_hdr_read(fp);
    bam1_t* read = bam_init1();
    EXPECT_TRUE(sam_read1(fp, hdr, read) >= 0);  // read0
    EXPECT_TRUE(ReadFilterLib::valid_alignment_start(read));
    EXPECT_TRUE(ReadFilterLib::is_well_formed(read, hdr));
    EXPECT_TRUE(sam_read1(fp, hdr, read) >= 0);  // read1
    EXPECT_TRUE(ReadFilterLib::valid_alignment_start(read));
    EXPECT_TRUE(ReadFilterLib::is_well_formed(read, hdr));
    EXPECT_TRUE(sam_read1(fp, hdr, read) >= 0);  // read2
    EXPECT_TRUE(ReadFilterLib::valid_alignment_end(read));
    EXPECT_TRUE(ReadFilterLib::is_well_formed(read, hdr));
    EXPECT_TRUE(sam_read1(fp, hdr, read) >= 0);  // read3
    EXPECT_TRUE(ReadFilterLib::valid_alignment_end(read));
    EXPECT_TRUE(ReadFilterLib::is_well_formed(read, hdr));
    EXPECT_TRUE(sam_read1(fp, hdr, read) >= 0);  // read4
    EXPECT_TRUE(ReadFilterLib::alignment_agree_with_hdr(read, hdr));
    EXPECT_TRUE(ReadFilterLib::is_well_formed(read, hdr));
    EXPECT_TRUE(sam_read1(fp, hdr, read) >= 0);  // read5
    EXPECT_TRUE(ReadFilterLib::alignment_agree_with_hdr(read, hdr));
    EXPECT_TRUE(ReadFilterLib::is_well_formed(read, hdr));
    EXPECT_TRUE(sam_read1(fp, hdr, read) >= 0);  // read6
    EXPECT_TRUE(ReadFilterLib::read_len_equal_cigar_len(read));
    EXPECT_TRUE(ReadFilterLib::is_well_formed(read, hdr));
    EXPECT_TRUE(sam_read1(fp, hdr, read) >= 0);  // read7
    EXPECT_TRUE(ReadFilterLib::read_len_equal_cigar_len(read));
    EXPECT_TRUE(ReadFilterLib::is_well_formed(read, hdr));
    EXPECT_TRUE(sam_read1(fp, hdr, read) >= 0);  // read8
    EXPECT_TRUE(ReadFilterLib::seq_is_stored(read));
    EXPECT_TRUE(ReadFilterLib::is_well_formed(read, hdr));
    EXPECT_TRUE(sam_read1(fp, hdr, read) >= 0);  // read9
    EXPECT_TRUE(ReadFilterLib::seq_is_stored(read));
    EXPECT_TRUE(ReadFilterLib::is_well_formed(read, hdr));
    bam_destroy1(read);
    bam_hdr_destroy(hdr);
    sam_close(fp);
}

TEST(ReadFilterTest, WellFormed1)
{
    const char sam[] =
        "data:,"
        "@SQ\tSN:chr1\tLN:249250621\n";
    samFile* fp = sam_open(sam, "r");
    bam_hdr_t* hdr = sam_hdr_read(fp);
    uint32_t cigar[] = {160};
    uint32_t cigar1[] = {128};
    uint32_t cigar2[] = {240};
    bam1_t* read0 = bam_init1();
    bam1_t* read1 = bam_init1();
    bam1_t* read2 = bam_init1();
    bam1_t* read3 = bam_init1();
    bam1_t* read4 = bam_init1();
    bam1_t* read5 = bam_init1();
    bam1_t* read6 = bam_init1();
    bam1_t* read7 = bam_init1();
    bam1_t* read8 = bam_init1();
    bam1_t* read9 = bam_init1();
    // invalid start.
    bam_set1(read0, 10, "read0", 0, 0, -1, 99, 1, cigar, 0, 0, 0, 10, "TACGTACGTA", "!@#$%^&*()", 0);
    // valid end.
    bam_set1(read1, 10, "read1", 0, 0, -1000, 99, 1, cigar, 0, 0, 0, 10, "TACGTACGTA", "!@#$%^&*()", 0);
    // agree with hdr
    bam_set1(read2, 10, "read2", 0, 0, 1, 99, 1, cigar, 0, 0, 0, 10, "TACGTACGTA", "!@#$%^&*()", 0);
    // disagree with hdr
    bam_set1(read3, 10, "read3", 0, 5, 1, 99, 1, cigar, 0, 0, 0, 10, "TACGTACGTA", "!@#$%^&*()", 0);
    // match baselen & cigarlen
    bam_set1(read4, 10, "read4", 0, 0, 1, 99, 1, cigar, 0, 0, 0, 10, "TACGTACGTA", "!@#$%^&*()", 0);
    // unmatch baselen & cigarlen
    bam_set1(read5, 10, "read5", 0, 0, 1, 99, 1, cigar1, 0, 0, 0, 8, "TACGTACGTA", "!@#$%^&*()", 0);
    // unmatch baselen & cigarlen
    bam_set1(read6, 10, "read6", 0, 0, 1, 99, 1, cigar2, 0, 0, 0, 15, "TACGTACGTA", "!@#$%^&*()", 0);
    // has seq
    bam_set1(read7, 10, "read7", 0, 0, 1, 99, 1, cigar, 0, 0, 0, 10, "TACGTACGTA", "!@#$%^&*()", 0);
    // has no seq
    bam_set1(read8, 10, "read8", 0, 0, 1, 99, 1, cigar, 0, 0, 0, 0, "TACGTACGTA", "!@#$%^&*()", 0);
    // normal
    bam_set1(read9, 10, "read9", 0, 0, 1, 99, 1, cigar, 0, 0, 0, 10, "TACGTACGTA", "!@#$%^&*()", 0);
    bam_aux_append(read9, "RG", 'Z', strlen("RovacaTEST") + 1, (uint8_t*)"RovacaTEST");
    EXPECT_FALSE(ReadFilterLib::valid_alignment_start(read0));
    EXPECT_TRUE(ReadFilterLib::valid_alignment_end(read1));
    EXPECT_TRUE(ReadFilterLib::alignment_agree_with_hdr(read2, hdr));
    EXPECT_FALSE(ReadFilterLib::alignment_agree_with_hdr(read3, hdr));
    EXPECT_TRUE(ReadFilterLib::read_len_equal_cigar_len(read4));
    EXPECT_TRUE(ReadFilterLib::read_len_equal_cigar_len(read5));
    EXPECT_TRUE(ReadFilterLib::read_len_equal_cigar_len(read6));
    EXPECT_TRUE(ReadFilterLib::seq_is_stored(read7));
    EXPECT_FALSE(ReadFilterLib::seq_is_stored(read8));
    EXPECT_TRUE(ReadFilterLib::is_well_formed(read9, hdr));
    bam_destroy1(read9);
    bam_destroy1(read8);
    bam_destroy1(read7);
    bam_destroy1(read6);
    bam_destroy1(read5);
    bam_destroy1(read4);
    bam_destroy1(read3);
    bam_destroy1(read2);
    bam_destroy1(read1);
    bam_destroy1(read0);
    bam_hdr_destroy(hdr);
    sam_close(fp);
}

TEST(ReadFilterTest, HasReadGroup)
{
    const char sam[] =
        "data:,"
        "@SQ\tSN:chr1\tLN:249250621\n"
        "read0\t0\tchr1\t10\t0\t10M\t*\t0\t0\tTACGTACGTA\t!@#$%^&*()\tRG:Z:NA12878\n"
        "read1\t0\tchr1\t80\t60\t10M\t*\t0\t0\tTACGTACGTA\t!@#$%^&*()\n";
    samFile* fp = sam_open(sam, "r");
    bam_hdr_t* hdr = sam_hdr_read(fp);
    bam1_t* read = bam_init1();
    EXPECT_TRUE(sam_read1(fp, hdr, read) >= 0);
    EXPECT_TRUE(ReadFilterLib::has_read_group(read));
    EXPECT_TRUE(sam_read1(fp, hdr, read) >= 0);
    EXPECT_FALSE(ReadFilterLib::has_read_group(read));
    bam_destroy1(read);
    bam_hdr_destroy(hdr);
    sam_close(fp);
}

TEST(ReadFilterTest, GoodCigar)
{
    const char sam[] =
        "data:,"
        "@SQ\tSN:chr1\tLN:249250621\n"
        "@RG\tID:LIBAI\tPL:COMPLETE\tPU:NA12878\tLB:NA12878\tSM:SF3\tCN:BGI\n"
        "read0\t0\tchr1\t10\t0\t1D10M\t*\t0\t0\tTACGTACGTA\t!@#$%^&*()\n"        // start with deletion
        "read1\t0\tchr1\t10\t0\t5S2D5M\t*\t0\t0\tTACGTACGTA\t!@#$%^&*()\n"       // start with deletion
        "read2\t0\tchr1\t80\t60\t10M1D\t*\t0\t0\tTACGTACGTA\t!@#$%^&*()\n"       // end with deletion
        "read3\t0\tchr1\t80\t60\t5M1D5S\t*\t0\t0\tTACGTACGTA\t!@#$%^&*()\n"      // end with deletion
        "read4\t0\tchr1\t10\t60\t3S3M3D3I1S\t*\t0\t0\tGTACGTACGT\t!@#$%^&*()\n"  // consective indel
        "read5\t0\tchr1\t30\t60\t3S3M3I3D1S\t*\t0\t0\tACGTACGTAC\t!@#$%^&*()\n"  // consective indel
        "read6\t0\tchr1\t56\t60\t5M5N5S\t*\t0\t0\tGTACGTACGT\t!@#$%^&*()\n"      // has N
        "read7\t0\tchr1\t70\t60\t10S5N\t*\t0\t0\tTACGTACGTA\t!@#$%^&*()\n"       // has no consume ref
        "read8\t0\tchr1\t40\t60\t5S5S\t*\t0\t0\tACGTACGTAC\t!@#$%^&*()\n"        // has no consume ref
        "read9\t0\tchr1\t20\t60\t10M\t*\t0\t0\tGTACGTACGT\t!@#$%^&*()\n"         // normal
        "read10\t0\tchr1\t20\t60\t15M\t*\t0\t0\tGTACGTACGTGTACG\tGTACGTACGTGTACG\n";
    samFile* fp = sam_open(sam, "r");
    bam_hdr_t* hdr = sam_hdr_read(fp);
    bam1_t* read = bam_init1();
    EXPECT_TRUE(sam_read1(fp, hdr, read) >= 0);  // read0
    EXPECT_FALSE(ReadFilterLib::is_good_cigar(read));
    EXPECT_TRUE(sam_read1(fp, hdr, read) >= 0);  // read1
    EXPECT_FALSE(ReadFilterLib::is_good_cigar(read));
    EXPECT_TRUE(sam_read1(fp, hdr, read) >= 0);  // read2
    EXPECT_FALSE(ReadFilterLib::is_good_cigar(read));
    EXPECT_TRUE(sam_read1(fp, hdr, read) >= 0);  // read3
    EXPECT_FALSE(ReadFilterLib::is_good_cigar(read));
    EXPECT_TRUE(sam_read1(fp, hdr, read) >= 0);  // read4
    EXPECT_FALSE(ReadFilterLib::is_good_cigar(read));
    EXPECT_TRUE(sam_read1(fp, hdr, read) >= 0);  // read5
    EXPECT_FALSE(ReadFilterLib::is_good_cigar(read));
    EXPECT_TRUE(sam_read1(fp, hdr, read) >= 0);  // read6
    EXPECT_FALSE(ReadFilterLib::is_good_cigar(read));
    EXPECT_TRUE(sam_read1(fp, hdr, read) >= 0);  // read7
    EXPECT_FALSE(ReadFilterLib::is_good_cigar(read));
    EXPECT_TRUE(sam_read1(fp, hdr, read) >= 0);  // read8
    EXPECT_FALSE(ReadFilterLib::is_good_cigar(read));
    EXPECT_TRUE(sam_read1(fp, hdr, read) >= 0);  // read9
    EXPECT_TRUE(ReadFilterLib::is_good_cigar(read));
    EXPECT_TRUE(sam_read1(fp, hdr, read) >= 0);  // read10
    EXPECT_TRUE(ReadFilterLib::is_good_cigar(read));
    bam_destroy1(read);
    bam_hdr_destroy(hdr);
    sam_close(fp);
}