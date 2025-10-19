#ifndef ROVACA_HC_UNIT_TEST_UTILS_H_
#define ROVACA_HC_UNIT_TEST_UTILS_H_
#include <random>
#include <string>

#include "allele.h"
#include "allele_likelihoods.hpp"
#include "bam_data_pool.hpp"
#include "forward.h"
#include "genotype_macors.h"
#include "genotype_struct.h"
#include "htslib/sam.h"
#include "indexed_allele_list.hpp"
#include "indexed_sample_list.hpp"
#include "rovaca_logger.h"
#include "read_record.h"
#include "utils/cigar_utils.h"
#include "utils/debug_utils.h"

using namespace rovaca;

namespace UnitTestUtils
{

#define GATK_RANDOM_SEED (47382911L)

static const size_t s_ref_len = 1300;
const char* ref =
    "GGTTGCAAAAATTTTCTCCCATTTTCTGGGTTGCCTGTTCACTCTGATGGTAGTTTCTTTTGCTGTGCAGAAGCTCTTTAGTTTAATTAGATCCCATTTGTCAATTTTGTCTTTTGTTGCCATTGCTTTT"
    "CACCGATCCCACAGAAATACAAACTACCATCAGAGAATACTACAAACACCTCTACGCAAATAAACTAGAAAATCTAGAAGAAATGGATAAATTCCTGGACACATACACTCTCCCAAGCCTAAACCAGGAA"
    "TTGAATCTCTGAATAGACCAATAACAGAAGCTGAAATTGTGGCAATAATCAATAGCTTACCAACCAAAAAGAGTCCAGGACCAGATGGATTCACAGCCGAATTCTACCAGAGGTACAAGGAGGAACTGGT"
    "TTCCTTCTGAAACTATTCCAATCAATAGAAAAAGAGGGAGTCCTCCCTAACTCATTTTATGAGGCCAGCATCATTCTGATACCAAAGCCAGGCAGAGACACAACAAAAAAAGAGAATTTTAGACCAATAT"
    "GATGAACATTGATGCAAAAATCCTCAATAAAATACTGGCAAAACGAATCCAGCAGCACATCAAAAAGCTTATCCACCAAGATCAAGTGGGCTTCATCCCTGGGATGCAAGGCTGGTTCAATATACGCAAA"
    "TAAATGTAATCCAGCATATAAACAGAGCCAAAGACAAAAACCACATGATTATCTCAATAGATGCAGAAAAGGCCTTTGACAAAATTCAACAACCCTTCATGCTAAAAACTCTCAATAAATTAGGTATTGA"
    "ACGTATTTCAAAATAATAAGAGCTATCTATGACAAACCCACAGCCAATATCATACTGAATGGGCAAAAACTGGAAGCATTCCCTTTGAAAACTGGCACAAGACAGGGATGCCCTCTCTCACCACTCCTAT"
    "CATAGTGTTGGAAGTTCTGGCCAGGGCAATTAGGCAGGAGAAGGAAATAAAGGGTATTCAGTTAGGAAAAGAGGAAGTCAAATTGTCCCTGTTTGCAGACGACATGATTGTATATCTAGAAAACCCCATT"
    "CAGCCCAAAATCTTCCTAAGCTGATAAGCAACTTCAGCAAAGTCTCAGGATACAAAATCAATGTACAAAAATCACAAGCATTCTTATACACCAACAACAGACAAACAGAGAGCCAAACCATGAGTGAACT"
    "TTCACAATTGTTTCAAAGAGAATAAAATACCTAGGAATCCAACTTACAAGGGACGTGAAGGACCTCTTCAAGGAGAACTACAAATCACTGCTCAAGGAAATAAAAGAGGATACAAAGAAATGGAAGAACA";

typedef struct SamHeader
{
    samFile* file;
    sam_hdr_t* header;
} SamHeader, *pSamHeader;

void sam_header_destroy(pSamHeader h)
{
    sam_hdr_destroy(h->header);
    sam_close(h->file);
    delete h;
}

pSamHeader create_artificial_sam_header(const char* file)
{
    pSamHeader h = new SamHeader{};
    h->file = sam_open(file, "r");
    if (!h->file) {
        RovacaLogger::error("couldn't open file");
        exit(EXIT_FAILURE);
    }

    h->header = sam_hdr_read(h->file);
    if (!h->header) {
        RovacaLogger::error("couldn't reading header from file");
        goto error;
    }
    return h;

error:
    sam_hdr_destroy(h->header);
    if (h->file) sam_close(h->file);
    exit(EXIT_FAILURE);
}

/*! @brief 创建 bam_hdr_t */
pSamHeader create_artificial_sam_header(int32_t number_of_chromosomes, int32_t starting_chromosome, uint32_t chromosome_size)
{
    static const char header_text[] =
        "data:,"
        "@HD\tVN:1.4\tGO:group\tSS:coordinate:queryname\n"
        "@CO\tThis line is good\n"
        "@CO\tThis line below will be updated\n";

    pSamHeader h = new SamHeader{};
    h->file = sam_open(header_text, "r");
    if (!h->file) {
        RovacaLogger::error("couldn't open file");
        exit(EXIT_FAILURE);
    }

    h->header = sam_hdr_read(h->file);
    if (!h->header) {
        RovacaLogger::error("couldn't reading header from file");
        goto error;
    }

    int32_t r;
    for (int32_t x = starting_chromosome; x < starting_chromosome + number_of_chromosomes; ++x) {
        r = sam_hdr_add_line(h->header, "SQ", "SN", std::to_string(x).c_str(), "LN", std::to_string(chromosome_size).c_str(), NULL);
        if (r < 0) {
            RovacaLogger::error("sam_hdr_add_line");
            goto error;
        }
    }

    return h;

error:
    sam_hdr_destroy(h->header);
    if (h->file) sam_close(h->file);
    exit(EXIT_FAILURE);
}

/*! @brief 使用bam1_t创建一条简易ReadRecord */
pReadRecord create_artificial_read(sam_hdr_t* header, const char* qname, int32_t tid, int64_t start, uint32_t length, pBamDataPool bpool,
                                   pMemoryPool mpool)
{
    bam1_t* b = bpool->alloc_bam_struct_pool();

    size_t l_qname = strlen(qname);
    uint16_t flag = 0;
    uint8_t mapq = 0;
    uint32_t cigar = length << BAM_CIGAR_SHIFT | BAM_CMATCH;
    int32_t mtid = -1;
    hts_pos_t mpos = 0;
    hts_pos_t isize = 0;
    pBases seq = new ALLOC_FLEXIBLE_IN_POOL(mpool, Bases, length, uint8_t) Bases{length};
    pBases qual = new ALLOC_FLEXIBLE_IN_POOL(mpool, Bases, length, uint8_t) Bases{length};

    for (uint32_t i = 0; i < length; ++i) {
        seq->data[i] = 'A';
        qual->data[i] = 'A';
    }

    int32_t r;
    r = bam_set1(b, l_qname, qname, flag, tid, start, mapq, 1, &cigar, mtid, mpos, isize, length, (char*)seq->data, (char*)qual->data, 0);
    if (ROVACA_UNLIKELY(r < 0)) {
        RovacaLogger::error("bam_set1 error");
        exit(EXIT_FAILURE);
    }

    pReadRecord result = ReadRecord::create(mpool, header, b);
    if (tid == -1) {
        result->set_read_unmapped_flag(true);
    }
    result->set_start(start);
    result->set_stop(start + length - 1);
    return result;
}

/*! @brief 使用bam1_t创建一条简易ReadRecord */
pReadRecord create_artificial_read(sam_hdr_t* header, const char* qname, int32_t tid, int64_t start, uint32_t length, pBamDataPool bpool,
                                   pMemoryPool mpool, const std::string& cigar, const std::string& seqs, const std::string& quals)
{
    bam1_t* b = bpool->alloc_bam_struct_pool();

    size_t l_qname = strlen(qname);
    uint16_t flag = 0;
    uint8_t mapq = 0;
    int32_t mtid = -1;
    hts_pos_t mpos = 0;
    hts_pos_t isize = 0;
    pBases seq = DebugUtils::str2bases(seqs, mpool);
    pBases qual = DebugUtils::str2bases(quals, mpool);
    pCigar c = DebugUtils::str2cigar(cigar, mpool);

    int32_t r;
    r = bam_set1(b, l_qname, qname, flag, tid, start, mapq, int32_t(c->num), c->data, mtid, mpos, isize, length, (char*)seq->data,
                 (char*)qual->data, 0);
    if (ROVACA_UNLIKELY(r < 0)) {
        RovacaLogger::error("bam_set1 error");
        exit(EXIT_FAILURE);
    }

    pReadRecord result = ReadRecord::create(mpool, header, b);
    if (tid == -1) {
        result->set_read_unmapped_flag(true);
    }
    result->set_start(start);
    result->set_stop(start + length - 1);
    return result;
}

/*! @brief 构建read_count个ReadRecord的集合 */
ReadVector read_list(pSamHeader h, int32_t sample_index, int32_t read_count, pBamDataPool bpool, pMemoryPool mpool)
{
    ReadVector reads{mpool};
    reads.reserve(read_count);
    for (int32_t j = 0, read_index = 0; j < read_count; ++j) {
        std::string read_name = "READ_" + std::to_string(sample_index) + "_" + std::to_string(read_index++);
        reads.push_back(create_artificial_read(h->header, read_name.c_str(), 1, 1, 100, bpool, mpool));
    }
    return reads;
}

/*! @brief 创建包含sample_count个样本的IndexedSampleList */
InterfaceSampleList* sample_list(int32_t sample_count)
{
    std::string sample{"SAMPLE_"};
    std::vector<std::string> samples{};
    samples.reserve(sample_count);
    for (int32_t i = 0; i < sample_count; ++i) {
        samples.push_back(sample + std::to_string(i));
    }
    return IndexedSampleList::create(samples);
}

/*! @brief 创建一个简易的IndexedAlleleList，内部等位基因为ATGC(伪)随机排列 */
InterfaceAlleleList<pAllele>* allele_list(int32_t allele_count, int32_t max_allele_length, bool skip_if_repeats, pMemoryPool pool)
{
    AlleleVector alleles{pool};
    alleles.reserve(allele_count);

    pBases bases;
    pAllele allele;
    uint32_t start, len;

    std::mt19937 gen{GATK_RANDOM_SEED};
    std::uniform_int_distribution<> len_dis(1, max_allele_length - 1);
    std::uniform_int_distribution<> start_dis(1, 1300 - max_allele_length - 10);
    for (int32_t a = 0; a < allele_count; ++a) {
        len = len_dis(gen);
        start = start_dis(gen);
        if (len == 1) {
            allele = Allele::create_allele(ref[start], a == 0);
        }
        else {
            bases = new ALLOC_FLEXIBLE_IN_POOL(pool, Bases, len, uint8_t) Bases{len};
            memcpy(bases->data, ref + start, len * sizeof(uint8_t));
            allele = Allele::create_allele(bases, a == 0, pool);
        }
        alleles.emplace_back(allele);
    }

    auto* result = IndexedAlleleList<pAllele>::create(alleles, pool);
    if (skip_if_repeats && result->number_of_alleles() != alleles.size()) {
        RovacaLogger::error("repeated alleles, should be infrequent");
        exit(EXIT_FAILURE);
    }
    return result;
}

/*! @brief produces a test likelihood depending on the sample, read and allele index. */
double test_likelihood(int32_t sample_index, int32_t allele_index, int32_t read_index)
{
    return -std::abs(3 * (sample_index + 1) + 7 * (allele_index + 1) + 11 * (read_index + 1));
}

/*!
 * @param h header生命周期贯穿测试
 * @param sample_list sample_list设计初衷希望全局仅一份，不会由内存池构建，使用完毕需要delete
 */
pRALikelihoods read_likelihoods(pSamHeader h, int32_t allele_count, const std::vector<int32_t>& read_count,
                                InterfaceSampleList* sample_list, pMemoryPool mpool, pBamDataPool bpool)
{
    int32_t sample_count = (int32_t)read_count.size();
    InterfaceAlleleList<pAllele>* al = allele_list(allele_count, 100, true, mpool);
    Int32ToReadVectorMap sample_to_reads{mpool};
    for (int32_t i = 0; i < sample_count; ++i) {
        sample_to_reads.insert({i, read_list(h, i, read_count.at(i), bpool, mpool)});
    }

    pRALikelihoods likelihoods = RALikelihoods::create<pReadRecord, pAllele>(mpool, sample_list, al, sample_to_reads);

    for (int32_t s = 0; s < sample_count; ++s) {
        auto* sample_likelihoods = likelihoods->sample_matrix((size_t)s);
        for (int32_t a = 0; a < allele_count; a++) {
            for (int32_t r = 0; r < read_count[s]; r++) {
                sample_likelihoods->set((size_t)a, (size_t)r, test_likelihood(s, a, r));
            }
        }
    }
    return likelihoods;
}

/*!
 * @brief std::string 转为 pBases形式
 */
pBases str2bases(const std::string& str, pMemoryPool pool)
{
    uint32_t len = uint32_t(str.size());
    pBases result = new ALLOC_FLEXIBLE_IN_POOL(pool, Bases, len, uint8_t) Bases{len};
    memcpy(result->data, str.data(), len * sizeof(uint8_t));
    return result;
}

}  // namespace UnitTestUtils

#endif  // ROVACA_HC_UNIT_TEST_UTILS_H_
