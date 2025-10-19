#include <gtest/gtest.h>

#include <cstring>
#include <fstream>
#include <iostream>

#include "bam_data_pool.hpp"
#include "block_combiner.h"
#include "forward.h"
#include "genotype_argument.h"
#include "germline_genotying_engine.h"
#include "homogeneous_ploidy_model.hpp"
#include "indexed_sample_list.hpp"
#include "rovaca_logger.h"
#include "unit_test_header.h"
#include "unit_test_utils.hpp"
#include "utils/adapter_utils.h"
#include "utils/debug_utils.h"
#include "variant.h"
#include "writer/writer.h"

using namespace std;
using namespace rovaca;

static constexpr size_t s_buffer_size = 1024 * 1024 * 500;

typedef struct TestData TestData, *pTestData;

static const std::string vcfInFile = "/data/rovaca-dev/zhangtiefeng/data/hc-genotype-test/source/NA12878-chr1-9960-30585-vcf-finalized.data";
static const std::string gvcfInFile =
    "/data/rovaca-dev/zhangtiefeng/data/hc-genotype-test/source/NA12878-chr1-9960-30585-gvcf-finalized.data";
static const std::string vcfOutFile =
    "/data/rovaca-dev/zhangtiefeng/data/hc-genotype-test/result/tmp/rovaca-NA12878-chr1-9960-30585-finalized.vcf";
static const std::string gvcfOutFile =
    "/data/rovaca-dev/zhangtiefeng/data/hc-genotype-test/result/tmp/rovaca-NA12878-chr1-9960-30585-finalized.g.vcf";

class GermlineGenotyingEngineUnitTest : public ::testing::Test
{
public:
    pTestData runData();

    static LineType parse_str(std::string &line);

protected:
    void SetUp() override;
    void TearDown() override;

    uint8_t *_buffer{};
    pMemoryPool _pool{};
    pBamDataPool _bampool{};

    UnitTestUtils::pSamHeader _header{};
    ifstream ifs;
    pHCArgs args{};
};

TEST_F(GermlineGenotyingEngineUnitTest, trimRegionTest)
{
    args->init_reference_confidence_mode(ReferenceConfidenceMode::NONE);
    ifs.open(vcfInFile.c_str(), ios_base::in);
    if (!ifs.is_open()) {
        std::cout << "open error" << endl;
        exit(EXIT_FAILURE);
    }

    while (true) {
        std::pmr::monotonic_buffer_resource my_pool(_buffer, s_buffer_size, std::pmr::null_memory_resource());
        _pool = std::addressof(my_pool);
        _bampool->finalize();
        pTestData data = runData();
        if (data == nullptr) {
            break;
        }

        IntervalPair trim_result = AdapterUtils::trim_region(data->original_haplotype, &data->ref_bases, data->ref_loc, data->original,
                                                             data->original_padded, args, _pool);
        std::string gatk_trim_variant = DebugUtils::interval2str(data->variant, _header->header);
        std::string rovaca_trim_variant = DebugUtils::interval2str(trim_result.first, _header->header);
        ASSERT_EQ(gatk_trim_variant, rovaca_trim_variant) << gatk_trim_variant << " / " << rovaca_trim_variant;

        std::string gatk_trim_variant_padded = DebugUtils::interval2str(data->variant_padded, _header->header);
        std::string rovaca_trim_variant_padded = DebugUtils::interval2str(trim_result.second, _header->header);
        ASSERT_EQ(gatk_trim_variant_padded, rovaca_trim_variant_padded) << gatk_trim_variant_padded << " / " << rovaca_trim_variant_padded;
    }
}

TEST_F(GermlineGenotyingEngineUnitTest, testFilterNonPassingReads1)
{
    args->init_reference_confidence_mode(ReferenceConfidenceMode::NONE);
    ifs.open(vcfInFile.c_str(), ios_base::in);
    if (!ifs.is_open()) {
        std::cout << "open error" << endl;
        exit(EXIT_FAILURE);
    }

    while (true) {
        std::pmr::monotonic_buffer_resource my_pool(_buffer, s_buffer_size, std::pmr::null_memory_resource());
        _pool = std::addressof(my_pool);
        _bampool->finalize();
        pTestData data = runData();
        if (data == nullptr) {
            break;
        }
        ReadHashSet ff{{data->trimed_reads.begin(), data->trimed_reads.end()}, _pool};
        ReadList filter = AdapterUtils::filter_non_passing_reads1(ff, args->minimum_read_length_after_trimming, _pool);
        ASSERT_EQ(filter.size(), data->filtered_reads1.size());

        std::unordered_set<std::string> result1;
        std::unordered_set<std::string> result2;
        for (pReadRecord read : filter) {
            result1.insert(DebugUtils::read2str(read, _header->header, _pool));
        }
        for (pReadRecord read : data->filtered_reads1) {
            result2.insert(DebugUtils::read2str(read, _header->header, _pool));
        }

        for (const std::string &s : result1) {
            ASSERT_TRUE(result2.count(s));
        }
    }
}

TEST_F(GermlineGenotyingEngineUnitTest, testFilterNonPassingReads2)
{
    args->init_reference_confidence_mode(ReferenceConfidenceMode::NONE);
    ifs.open(vcfInFile.c_str(), ios_base::in);
    if (!ifs.is_open()) {
        std::cout << "open error" << endl;
        exit(EXIT_FAILURE);
    }

    while (true) {
        std::pmr::monotonic_buffer_resource my_pool(_buffer, s_buffer_size, std::pmr::null_memory_resource());
        _pool = std::addressof(my_pool);
        _bampool->finalize();
        pTestData data = runData();
        if (data == nullptr) {
            break;
        }
        ReadHashSet ff{{data->trimed_reads.begin(), data->trimed_reads.end()}, _pool};
        ReadList filter = AdapterUtils::filter_non_passing_reads1(ff, args->minimum_read_length_after_trimming, _pool);
        ASSERT_EQ(filter.size(), data->filtered_reads1.size());

        ReadList filter2 = AdapterUtils::filter_non_passing_reads2(ff, args->mapping_quality_threshold, _pool);
        ASSERT_EQ(filter2.size(), data->filtered_reads2.size());
        ASSERT_EQ(ff.size(), data->pairhmm_reads.size());

        std::unordered_set<std::string> result1;
        std::unordered_set<std::string> result2;
        for (pReadRecord read : filter2) {
            result1.insert(DebugUtils::read2str(read, _header->header, _pool));
        }
        for (pReadRecord read : data->filtered_reads2) {
            result2.insert(DebugUtils::read2str(read, _header->header, _pool));
        }

        for (const std::string &s : result1) {
            ASSERT_TRUE(result2.count(s));
        }
    }
}

TEST_F(GermlineGenotyingEngineUnitTest, trimReadsTest)
{
    args->init_reference_confidence_mode(ReferenceConfidenceMode::NONE);
    ifs.open(vcfInFile.c_str(), ios_base::in);
    if (!ifs.is_open()) {
        std::cout << "open error" << endl;
        exit(EXIT_FAILURE);
    }

    while (true) {
        std::pmr::monotonic_buffer_resource my_pool(_buffer, s_buffer_size, std::pmr::null_memory_resource());
        _pool = std::addressof(my_pool);
        _bampool->finalize();
        std::pmr::list<bam1_t *> extra_memory_reads{_pool};
        pTestData data = runData();
        if (data == nullptr) {
            break;
        }

        ReadHashSet ff{{data->original_reads.begin(), data->original_reads.end()}, _pool};
        ReadHashSet trimed_reads = AdapterUtils::trim_reads_by_region(ff, data->variant_padded, _pool, _bampool, extra_memory_reads);
        ASSERT_EQ(trimed_reads.size(), data->trimed_reads.size());

        std::set<std::string> result1;
        std::set<std::string> result2;
        for (pReadRecord read : trimed_reads) {
            result1.insert(DebugUtils::read2str(read, _header->header, _pool));
        }
        for (pReadRecord read : data->trimed_reads) {
            result2.insert(DebugUtils::read2str(read, _header->header, _pool));
        }

        for (const std::string &s : result1) {
            ASSERT_TRUE(result2.count(s)) << s;
        }
    }
}

TEST_F(GermlineGenotyingEngineUnitTest, trimHaplotypeTest)
{
    args->init_reference_confidence_mode(ReferenceConfidenceMode::NONE);
    ifs.open(vcfInFile.c_str(), ios_base::in);
    if (!ifs.is_open()) {
        std::cout << "open error" << endl;
        exit(EXIT_FAILURE);
    }

    while (true) {
        std::pmr::monotonic_buffer_resource my_pool(_buffer, s_buffer_size, std::pmr::null_memory_resource());
        _pool = std::addressof(my_pool);
        _bampool->finalize();
        pTestData data = runData();
        if (data == nullptr) {
            break;
        }

        HaplotypeVector trimed_haps = AdapterUtils::trim_haplotype_by_region(data->original_haplotype, data->variant_padded, _pool);
        ASSERT_EQ(trimed_haps.size(), data->trimed_haplotype.size());

        std::set<std::string> result1;
        std::set<std::string> result2;
        for (pHaplotype h : trimed_haps) {
            result1.insert(DebugUtils::haplotype2str(h, _header->header));
        }
        for (pHaplotype h : data->trimed_haplotype) {
            result2.insert(DebugUtils::haplotype2str(h, _header->header));
        }

        for (const std::string &s : result1) {
            ASSERT_TRUE(result2.count(s)) << s;
        }
    }
}

TEST_F(GermlineGenotyingEngineUnitTest, GvcfModeTest)
{
    auto task_start = std::chrono::high_resolution_clock::now();
    args->init_reference_confidence_mode(ReferenceConfidenceMode::GVCF);
    bool is_gvcf = args->reference_confidence_mode == ReferenceConfidenceMode::GVCF;
    args->tool_name = "LUSHVariantCaller";
    args->command_line = "LUSHVariantCaller haha haha";
    for (int32_t i = 1; i <= 60; ++i) {
        args->gvcf_gq_bands.push_back(i);
    }
    args->gvcf_gq_bands.push_back(70);
    args->gvcf_gq_bands.push_back(80);
    args->gvcf_gq_bands.push_back(90);
    args->gvcf_gq_bands.push_back(99);

    GermlineGenotyingEngine engine{};

    bcf_hdr_t *vcf_hdr = Writer::init_gvcf_header(args, _header->header);
    CHECK_CONDITION_EXIT(vcf_hdr == nullptr, "init_gvcf_hdr");

    htsFile *out_file = hts_open(gvcfOutFile.c_str(), "w");
    CHECK_CONDITION_EXIT(out_file == nullptr, "hts_open");
    ifs.open(gvcfInFile.c_str(), ios_base::in);
    if (!ifs.is_open()) {
        std::cout << "open error: " << std::strerror(errno) << endl;
        exit(EXIT_FAILURE);
    }

    BlockCombiner blocker{args->gvcf_gq_bands, vcf_hdr, _header->header, false, is_gvcf};

    CHECK_CONDITION_EXIT(bcf_hdr_write(out_file, vcf_hdr) != 0, "Failed to write header");

    pInterfaceSampleList samples = IndexedSampleList::create({"SF3"});
    pInterfacePloidyModel pm = HomogeneousPloidyModel::create(args->sample_ploidy, samples);
    VariantVector ss;

    kstring_t s{0, 0, nullptr};

    while (true) {
        std::pmr::monotonic_buffer_resource my_pool(_buffer, s_buffer_size, std::pmr::new_delete_resource());
        _pool = std::addressof(my_pool);
        _bampool->finalize();
        std::pmr::list<bam1_t *> extra_memory_reads{_pool};

        pTestData data = runData();
        if (data == nullptr) {
            break;
        }

        engine.init_engine_per_loop(args, _pool, _bampool, _header->header, vcf_hdr, samples, pm);

        if (data->variant == nullptr) {
            ReadHashSet original_reads{{data->original_reads.begin(), data->original_reads.end()}, _pool};
            ss = engine.reference_model_for_no_variation(&data->ref_bases, data->ref_loc, data->original, data->original_padded, 2,
                                                         original_reads);
            for (pVariant vc : ss) {
                blocker.submit(vc, &s);
            }
        }
        else {
            DoubleVector3D values{_pool};
            values.push_back(std::move(data->likelihoods));

            ReadVector reads{data->genotype_reads, _pool};
            Int32ToReadVectorMap evidence{_pool};
            evidence.insert({0, std::move(reads)});

            Int32ToReadVectorMap filtered_read{_pool};
            ReadVector f_read{data->filtered_reads2.begin(), data->filtered_reads2.end(), _pool};
            filtered_read.insert({0, std::move(f_read)});

            auto *allele_list = IndexedAlleleList<pHaplotype>::create(data->trimed_haplotype, _pool);

            auto *rh_likelihoods = RHLikelihoods::create<pReadRecord, pHaplotype>(_pool, samples, allele_list, evidence, std::move(values));

            auto result = engine.assign_genotype_likelihoods(rh_likelihoods, &data->ref_bases, data->ref_loc, data->variant, filtered_read);

            if (is_gvcf) {
                if (!GermlineGenotyingEngine::contains_calls(result.first)) {
                    ReadHashSet original_reads{{data->original_reads.begin(), data->original_reads.end()}, _pool};
                    ss = engine.reference_model_for_no_variation(&data->ref_bases, data->ref_loc, data->original, data->original_padded, 2,
                                                                 original_reads);
                }
                else {
                    ReadHashSet original_reads{{data->original_reads.begin(), data->original_reads.end()}, _pool};
                    ReadHashSet genotype_reads{{data->genotype_reads.begin(), data->genotype_reads.end()}, _pool};
                    pHaplotype ref_haplotype = allele_list->get_allele(size_t(allele_list->index_of_reference()));
                    pSimpleInterval o = data->original;
                    pSimpleInterval op = data->original_padded;
                    pSimpleInterval v = data->variant;
                    pSimpleInterval vp = data->variant_padded;
                    pRefFragment ref = &data->ref_bases;
                    ss = engine.call_non_active_site(ref_haplotype, ref, data->ref_loc, o, op, v, vp, result, original_reads,
                                                     genotype_reads, extra_memory_reads);
                }
            }
            const VariantVector &res = is_gvcf ? ss : result.first;
            for (pVariant vc : res) {
                blocker.submit(vc, &s);
            }
        }
    }

    auto write_start = std::chrono::high_resolution_clock::now();
    blocker.force_output(&s);
    CHECK_CONDITION_EXIT(vcf_write_line(out_file, &s) != 0, "vcf_write_line");
    ks_free(&s);

    auto end = std::chrono::high_resolution_clock::now();
    auto task_time = std::chrono::duration_cast<std::chrono::milliseconds>(end - task_start);
    auto write_time = std::chrono::duration_cast<std::chrono::milliseconds>(end - write_start);

    std::cout << "Total taking time [" << task_time.count() << "] milliseconds" << endl;
    std::cout << "Write taking time [" << write_time.count() << "] milliseconds" << endl;

    delete pm;
    delete samples;

    hts_close(out_file);
    bcf_hdr_destroy(vcf_hdr);
}

TEST_F(GermlineGenotyingEngineUnitTest, VcfModeTest)
{
    auto task_start = std::chrono::high_resolution_clock::now();

    args->init_reference_confidence_mode(ReferenceConfidenceMode::NONE);
    bool is_gvcf = args->reference_confidence_mode == ReferenceConfidenceMode::GVCF;
    args->tool_name = "LUSHVariantCaller";
    args->command_line = "LUSHVariantCaller haha haha";

    GermlineGenotyingEngine engine{};

    bcf_hdr_t *vcf_hdr = Writer::init_vcf_header(args, _header->header);
    CHECK_CONDITION_EXIT(vcf_hdr == nullptr, "init_vcf_hdr");
    bcf_hdr_add_sample(vcf_hdr, "SF3");

    htsFile *out_file = hts_open(vcfOutFile.c_str(), "w");
    CHECK_CONDITION_EXIT(out_file == nullptr, "hts_open");
    ifs.open(vcfInFile.c_str(), ios_base::in);
    if (!ifs.is_open()) {
        std::cout << "open error: " << std::strerror(errno) << endl;
        exit(EXIT_FAILURE);
    }

    BlockCombiner blocker{args->gvcf_gq_bands, vcf_hdr, _header->header, false, is_gvcf};

    CHECK_CONDITION_EXIT(bcf_hdr_write(out_file, vcf_hdr) != 0, "Failed to write header");

    pInterfaceSampleList samples = IndexedSampleList::create({"SF3"});
    pInterfacePloidyModel pm = HomogeneousPloidyModel::create(args->sample_ploidy, samples);
    VariantVector ss;
    kstring_t s{0, 0, nullptr};

    while (true) {
        std::pmr::monotonic_buffer_resource my_pool(_buffer, s_buffer_size, std::pmr::new_delete_resource());
        _pool = std::addressof(my_pool);
        _bampool->finalize();
        std::pmr::list<bam1_t *> extra_memory_reads{_pool};
        pTestData data = runData();
        if (data == nullptr) {
            break;
        }

        engine.init_engine_per_loop(args, _pool, _bampool, _header->header, vcf_hdr, samples, pm);

        if (data->variant == nullptr) {
            ReadHashSet original_reads{{data->original_reads.begin(), data->original_reads.end()}, _pool};
            ss = engine.reference_model_for_no_variation(&data->ref_bases, data->ref_loc, data->original, data->original_padded, 2,
                                                         original_reads);
            for (pVariant vc : ss) {
                blocker.submit(vc, &s);
            }
        }
        else {
            DoubleVector3D values{_pool};
            values.push_back(std::move(data->likelihoods));

            ReadVector reads{data->genotype_reads, _pool};
            Int32ToReadVectorMap evidence{_pool};
            evidence.insert({0, std::move(reads)});

            Int32ToReadVectorMap filtered_read{_pool};
            ReadVector f_read{data->filtered_reads2.begin(), data->filtered_reads2.end(), _pool};
            filtered_read.insert({0, std::move(f_read)});

            auto *allele_list = IndexedAlleleList<pHaplotype>::create(data->trimed_haplotype, _pool);

            auto *rh_likelihoods = RHLikelihoods::create<pReadRecord, pHaplotype>(_pool, samples, allele_list, evidence, std::move(values));

            auto result = engine.assign_genotype_likelihoods(rh_likelihoods, &data->ref_bases, data->ref_loc, data->variant, filtered_read);

            if (is_gvcf) {
                if (!GermlineGenotyingEngine::contains_calls(result.first)) {
                    ReadHashSet original_reads{{data->original_reads.begin(), data->original_reads.end()}, _pool};
                    ss = engine.reference_model_for_no_variation(&data->ref_bases, data->ref_loc, data->original, data->original_padded, 2,
                                                                 original_reads);
                }
                else {
                    ReadHashSet original_reads{{data->original_reads.begin(), data->original_reads.end()}, _pool};
                    ReadHashSet genotype_reads{{data->genotype_reads.begin(), data->genotype_reads.end()}, _pool};
                    pHaplotype ref_haplotype = allele_list->get_allele(size_t(allele_list->index_of_reference()));
                    pSimpleInterval o = data->original;
                    pSimpleInterval op = data->original_padded;
                    pSimpleInterval v = data->variant;
                    pSimpleInterval vp = data->variant_padded;
                    pRefFragment ref = &data->ref_bases;
                    ss = engine.call_non_active_site(ref_haplotype, ref, data->ref_loc, o, op, v, vp, result, original_reads,
                                                     genotype_reads, extra_memory_reads);
                }
            }
            const VariantVector &res = is_gvcf ? ss : result.first;
            for (pVariant vc : res) {
                blocker.submit(vc, &s);
            }
        }
    }

    auto write_start = std::chrono::high_resolution_clock::now();
    blocker.force_output(&s);
    CHECK_CONDITION_EXIT(vcf_write_line(out_file, &s) != 0, "vcf_write_line");
    ks_free(&s);

    auto end = std::chrono::high_resolution_clock::now();
    auto task_time = std::chrono::duration_cast<std::chrono::milliseconds>(end - task_start);
    auto write_time = std::chrono::duration_cast<std::chrono::milliseconds>(end - write_start);

    std::cout << "Total taking time [" << task_time.count() << "] milliseconds" << endl;
    std::cout << "Write taking time [" << write_time.count() << "] milliseconds" << endl;

    delete pm;
    delete samples;

    hts_close(out_file);
    bcf_hdr_destroy(vcf_hdr);
}

void GermlineGenotyingEngineUnitTest::SetUp()
{
    args = new GenotypeArgument{};
    _buffer = new uint8_t[s_buffer_size]{};
    _bampool = new BamDataPool(s_buffer_size);
    _header = UnitTestUtils::create_artificial_sam_header("/data/rovaca-dev/zhangtiefeng/data//NA12878/NA12878.sort.dup.bam");
}

void GermlineGenotyingEngineUnitTest::TearDown()
{
    delete _bampool;
    delete[] _buffer;
    UnitTestUtils::sam_header_destroy(_header);

    ifs.close();
    delete args;
}

// clang-format off
/*!
 * @brief 输出文件格式
 * num original_reads_num trimed_reads_num filtered_reads1_num pairhmm_reads_num filtered_reads2_num genotype_reads_num original_haplotype_num trimed_haplotype_num
 * original region
 * original_padded region
 * variant region
 * variant_padded region
 * ref_loc region
 * ref_bases ref_bases
 * original_reads       // 组装出来的reads
 * trimed_reads         // 经过第一次trim
 * filtered_reads1      // 第一次过滤，过滤出来的read直接删除，无需保留
 * pairhmm_reads        // 对trimed_reads进行过滤，得到的reads进行pairhmm
 * filtered_reads2      // 第二次过滤，保留过滤的read，后续会使用
 * genotype_reads       // pairhmm对低likelihoods得reads过滤，得到的reads进行genotype
 * original_haplotype
 * trimed_haplotype
 * likelihoods likelihoods
 */
// clang-format on
pTestData GermlineGenotyingEngineUnitTest::runData()
{
    int32_t original_reads_num{INVALID_INT};
    int32_t trimed_reads_num{INVALID_INT};
    int32_t filtered_reads1_num{INVALID_INT};
    int32_t pairhmm_reads_num{INVALID_INT};
    int32_t filtered_reads2_num{INVALID_INT};
    int32_t genotype_reads_num{INVALID_INT};
    int32_t original_haplotype_num{INVALID_INT};
    int32_t trimed_haplotype_num{INVALID_INT};

    std::string line;
    pTestData result{nullptr};
    while (std::getline(ifs, line)) {
        if (line.empty()) {
            break;
        }
        LineType type = parse_str(line);
        switch (type) {
            case k_num: {
                std::istringstream is{line};

                is >> original_reads_num >> trimed_reads_num >> filtered_reads1_num >> pairhmm_reads_num >> filtered_reads2_num >>
                    genotype_reads_num >> original_haplotype_num >> trimed_haplotype_num;

                if (original_reads_num == INVALID_INT || trimed_reads_num == INVALID_INT || filtered_reads1_num == INVALID_INT ||
                    pairhmm_reads_num == INVALID_INT || filtered_reads2_num == INVALID_INT || genotype_reads_num == INVALID_INT ||
                    original_haplotype_num == INVALID_INT || trimed_haplotype_num == INVALID_INT) {
                    RovacaLogger::error("num error");
                    exit(EXIT_FAILURE);
                }

                result = new ALLOC_TYPE_IN_POOL(_pool, TestData) TestData{original_reads_num,     trimed_reads_num,     filtered_reads1_num,
                                                                          pairhmm_reads_num,      filtered_reads2_num,  genotype_reads_num,
                                                                          original_haplotype_num, trimed_haplotype_num, _pool};

                break;
            }
            case k_original: {
                result->original = DebugUtils::str2interval(line, _header->header, _pool);
                break;
            }
            case k_original_padded: {
                result->original_padded = DebugUtils::str2interval(line, _header->header, _pool);
                break;
            }
            case k_variant: {
                result->variant = DebugUtils::str2interval(line, _header->header, _pool);
                break;
            }
            case k_variant_padded: {
                result->variant_padded = DebugUtils::str2interval(line, _header->header, _pool);
                break;
            }
            case k_ref_loc: {
                result->ref_loc = DebugUtils::str2interval(line, _header->header, _pool);
                break;
            }
            case k_ref_bases: {
                result->ref_bases.len = uint32_t(line.size());
                result->ref_bases.data = (uint8_t *)_pool->allocate(line.size());
                memcpy(result->ref_bases.data, line.c_str(), result->ref_bases.len * sizeof(uint8_t));
                break;
            }
            case k_original_reads: {
                result->original_reads.push_back(DebugUtils::str2read(line, _header->header, _pool, _bampool));
                break;
            }
            case k_trimed_reads: {
                result->trimed_reads.push_back(DebugUtils::str2read(line, _header->header, _pool, _bampool));
                break;
            }
            case k_filtered_reads1: {
                result->filtered_reads1.push_back(DebugUtils::str2read(line, _header->header, _pool, _bampool));
                break;
            }
            case k_filtered_reads2: {
                result->filtered_reads2.push_back(DebugUtils::str2read(line, _header->header, _pool, _bampool));
                break;
            }
            case k_pairhmm_reads: {
                result->pairhmm_reads.push_back(DebugUtils::str2read(line, _header->header, _pool, _bampool));
                break;
            }
            case k_genotype_reads: {
                result->genotype_reads.push_back(DebugUtils::str2read(line, _header->header, _pool, _bampool));
                break;
            }
            case k_original_haplotype: {
                result->original_haplotype.push_back(DebugUtils::str2haplotype(line, _header->header, _pool));
                break;
            }
            case k_trimed_haplotype: {
                result->trimed_haplotype.push_back(DebugUtils::str2haplotype(line, _header->header, _pool));
                break;
            }
            case k_likelihoods: {
                std::istringstream is{line};
                for (int32_t row = 0; row < original_haplotype_num; ++row) {
                    for (int32_t col = 0; col < genotype_reads_num; ++col) {
                        is >> result->likelihoods[row][col];
                    }
                }
                break;
            }
            default: {
                RovacaLogger::error("invalid line: {}", line.c_str());
                exit(EXIT_FAILURE);
            }
        }
    }

    return result;
}

LineType GermlineGenotyingEngineUnitTest::parse_str(std::string &line)
{
    std::string::size_type pos = line.find_first_of('\t');
    CHECK_CONDITION_EXIT(pos == std::string::npos, "not found");
    std::string first_word = line.substr(0, pos);

    if (!s_str2type.count(first_word)) {
        RovacaLogger::error("not found type: {}", first_word.c_str());
        exit(EXIT_FAILURE);
    }

    line = line.substr(pos + 1);
    return s_str2type.at(first_word);
}
