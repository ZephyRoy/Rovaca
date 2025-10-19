#include <gtest/gtest.h>

#include <random>
#include <tuple>

#include "allele.h"
#include "allele_likelihoods.hpp"
#include "bam_data_pool.hpp"
#include "forward.h"
#include "math_utils.h"
#include "unit_test_utils.hpp"

using namespace std;
using namespace rovaca;

static constexpr double s_epsilon = 1e-6;
static constexpr int64_t s_odd_read_start = 101;
static constexpr int64_t s_even_read_start = 1;
static constexpr size_t s_buffer_size = 1024 * 1024 * 10;

vector<vector<string>> s_samples{{"A", "B", "C"}, {"A"}, {"C", "A", "D", "E", "Salsa", "Gazpacho"}};

typedef tuple<const vector<string>&, const AlleleVector&, const AlleleVector&, Int32ToReadVectorMap, pmr::map<pAllele, AlleleVector>> Stup;
typedef pmr::vector<Stup> TupleVector;

class AlleleLikelihoodsUnitTest : public ::testing::Test
{
protected:
    void SetUp() override
    {
        _buffer = new uint8_t[s_buffer_size]{};
        _pool = new std::pmr::monotonic_buffer_resource(_buffer, s_buffer_size, std::pmr::null_memory_resource());
        _bampool = new BamDataPool(s_buffer_size);
        _header = UnitTestUtils::create_artificial_sam_header(10, 0, 1000);
    }

    void TearDown() override
    {
        delete _bampool;
        delete _pool;
        delete[] _buffer;
        sam_header_destroy(_header);
    }

    pmr::vector<pmr::vector<pAllele>> allele_sets();
    pmr::map<int32_t, pmr::vector<pReadRecord>> data_set_reads(const vector<string>& samples);

    void test_sample_queries(const vector<string>& samples, const pmr::map<int32_t, pmr::vector<pReadRecord>>& reads,
                             pRALikelihoods result);
    void test_allele_queries(const pmr::vector<pAllele>& alleles, pRALikelihoods result);
    static void test_likelihood_matrix_queries(const vector<string>& samples, pRALikelihoods result, const DoubleVector3D& likelihoods);

    DoubleVector3D fill_with_random_likelihoods(const vector<string>& samples, size_t allele_count, pRALikelihoods result);

    static void check_evidence_to_index_map_is_correct(pRALikelihoods likelihoods);

    TupleVector marginalization_data_sets();
    pmr::map<pAllele, AlleleVector> random_allele_map(const AlleleVector& from_alleles, const AlleleVector& to_alleles);

    uint8_t* _buffer{};
    pMemoryPool _pool{};
    pBamDataPool _bampool{};

    pInterfaceSampleList _samples{};
    UnitTestUtils::pSamHeader _header{};
};

TEST_F(AlleleLikelihoodsUnitTest, testInstantiationAndQuery)
{
    pmr::vector<pmr::vector<pAllele>> all_alleles = allele_sets();
    for (const vector<string>& samples : s_samples) {
        _samples = IndexedSampleList::create(samples);
        for (const pmr::vector<pAllele>& alleles : all_alleles) {
            pmr::map<int32_t, pmr::vector<pReadRecord>> reads = data_set_reads(samples);
            auto* allele_list = IndexedAlleleList<pAllele>::create(alleles, _pool);
            pRALikelihoods result = RALikelihoods::create<pReadRecord, pAllele>(_pool, _samples, allele_list, reads);

            ASSERT_EQ(result->number_of_samples(), samples.size());
            ASSERT_EQ(result->number_of_alleles(), alleles.size());
            ASSERT_EQ(result->get_alleles().size(), alleles.size());

            test_sample_queries(samples, reads, result);
            test_allele_queries(alleles, result);
        }
        delete _samples;
    }
}

TEST_F(AlleleLikelihoodsUnitTest, testLikelihoodFillingAndQuery)
{
    pmr::vector<pmr::vector<pAllele>> all_alleles = allele_sets();
    for (const vector<string>& samples : s_samples) {
        _samples = IndexedSampleList::create(samples);
        for (const pmr::vector<pAllele>& alleles : all_alleles) {
            pmr::map<int32_t, pmr::vector<pReadRecord>> reads = data_set_reads(samples);
            auto* allele_list = IndexedAlleleList<pAllele>::create(alleles, _pool);
            pRALikelihoods result = RALikelihoods::create<pReadRecord, pAllele>(_pool, _samples, allele_list, reads);
            DoubleVector3D likelihoods = fill_with_random_likelihoods(samples, alleles.size(), result);

            test_likelihood_matrix_queries(samples, result, likelihoods);
        }
        delete _samples;
    }
}

TEST_F(AlleleLikelihoodsUnitTest, testBestAlleles)
{
    pmr::vector<pmr::vector<pAllele>> all_alleles = allele_sets();
    for (const vector<string>& samples : s_samples) {
        _samples = IndexedSampleList::create(samples);
        for (const pmr::vector<pAllele>& alleles : all_alleles) {
            pmr::map<int32_t, pmr::vector<pReadRecord>> reads = data_set_reads(samples);
            auto* allele_list = IndexedAlleleList<pAllele>::create(alleles, _pool);
            pRALikelihoods original = RALikelihoods::create<pReadRecord, pAllele>(_pool, _samples, allele_list, reads);
            fill_with_random_likelihoods(samples, alleles.size(), original);

            size_t sample_count = samples.size();
            size_t number_of_alleles = alleles.size();
            int32_t ref_index = original->index_of_reference();
            pAllele ref_allele = ref_index >= 0 ? original->get_allele((size_t)ref_index) : nullptr;
            for (size_t s = 0; s < sample_count; ++s) {
                uint32_t sample_read_count = original->sample_evidence_count(s);
                auto* sample_matrix = original->sample_matrix(s);
                DoubleVector best_lk_array{sample_read_count, _pool};
                Int32Vector best_index_array{sample_read_count, _pool};
                DoubleVector confidence_array{sample_read_count, _pool};
                for (uint32_t r = 0; r < sample_read_count; ++r) {
                    int32_t bestindex_of_allele = -1;
                    double best_allele_lk = NEGATIVE_INFINITY;
                    double second_best_allele_lk = NEGATIVE_INFINITY;
                    for (int32_t a = 0; a < (int32_t)number_of_alleles; ++a) {
                        double lk = sample_matrix->get((size_t)a, (size_t)r);
                        if (lk > best_allele_lk) {
                            second_best_allele_lk = best_allele_lk;
                            best_allele_lk = lk;
                            bestindex_of_allele = a;
                        }
                        else if (lk > second_best_allele_lk) {
                            second_best_allele_lk = lk;
                        }
                    }
                    best_lk_array[r] = best_allele_lk;
                    confidence_array[r] = best_allele_lk - second_best_allele_lk;
                    best_index_array[r] = bestindex_of_allele;
                }
                auto best_alleles = original->best_alleles_breaking_ties();
                for (auto* best_allele : best_alleles) {
                    int32_t read_index = original->evidence_index(s, best_allele->evidence);
                    if (read_index == -1) {
                        continue;
                    }
                    double ref_likelihood = ref_index >= 0 ? sample_matrix->get((size_t)ref_index, (size_t)read_index) : NEGATIVE_INFINITY;
                    bool ref_override = ref_index >= 0 && ref_index != best_index_array[read_index] &&
                                        best_lk_array[read_index] - ref_likelihood < RALikelihoods::s_log_10_informative_threshold;
                    ASSERT_EQ(ref_override ? ref_likelihood : best_lk_array[read_index], best_allele->likelihood) << "s=" << s;
                    ASSERT_EQ(best_allele->best_allele, ref_override ? ref_allele : alleles[best_index_array[read_index]]) << "s=" << s;
                    if (!std::isnan(best_allele->confidence) && !std::isinf(best_allele->confidence)) {
                        ASSERT_NEAR(best_allele->confidence,
                                    ref_override ? ref_likelihood - best_lk_array[read_index] : confidence_array[read_index], s_epsilon)
                            << "s=" << s;
                    }
                }
            }
        }
        delete _samples;
    }
}

TEST_F(AlleleLikelihoodsUnitTest, testFilterReadsToOverlap)
{
    pmr::vector<pmr::vector<pAllele>> all_alleles = allele_sets();
    for (const vector<string>& samples : s_samples) {
        _samples = IndexedSampleList::create(samples);
        for (const pmr::vector<pAllele>& alleles : all_alleles) {
            pmr::map<int32_t, pmr::vector<pReadRecord>> reads = data_set_reads(samples);
            auto* allele_list = IndexedAlleleList<pAllele>::create(alleles, _pool);

            pRALikelihoods original = RALikelihoods::create<pReadRecord, pAllele>(_pool, _samples, allele_list, reads);
            fill_with_random_likelihoods(samples, alleles.size(), original);
            pRALikelihoods result = RALikelihoods::create<pReadRecord, pAllele>(_pool, _samples, allele_list, reads);
            fill_with_random_likelihoods(samples, alleles.size(), result);

            pSimpleInterval even_read_overlap = SimpleInterval::create(0, 1, 1, _pool);

            check_evidence_to_index_map_is_correct(result);
            result->retain_evidence(even_read_overlap);
            check_evidence_to_index_map_is_correct(result);

            DoubleVector3D new_likelihoods{samples.size(), _pool};
            std::for_each(new_likelihoods.begin(), new_likelihoods.end(), [&](DoubleVector2D& d) { d.resize(alleles.size()); });
            for (int32_t s = 0; s < (int32_t)samples.size(); s++) {
                for (int32_t a = 0; a < (int32_t)alleles.size(); a++) {
                    new_likelihoods[s][a].resize((original->sample_evidence_count(size_t(s)) + 1) / 2);
                    auto* sample_matrix = original->sample_matrix(size_t(s));
                    for (int32_t r = 0; r < (int32_t)new_likelihoods[s][a].size(); r++) {
                        ASSERT_EQ(result->evidence_index(size_t(s), sample_matrix->get_evidence(size_t(r << 1))), r)
                            << "s=" << s << " r=" << r;
                        new_likelihoods[s][a][r] = sample_matrix->get(size_t(a), size_t(r) << 1);
                    }
                }
            }
            test_likelihood_matrix_queries(samples, result, new_likelihoods);
        }
        delete _samples;
    }
}

TEST_F(AlleleLikelihoodsUnitTest, testMarginalization)
{
    TupleVector data = marginalization_data_sets();
    for (const auto& tup : data) {
        const vector<string>& samples = get<0>(tup);
        const AlleleVector& all_old = get<1>(tup);
        const AlleleVector& new_alleles = get<2>(tup);
        const Int32ToReadVectorMap& reads = get<3>(tup);
        const pmr::map<pAllele, AlleleVector>& allele_map = get<4>(tup);

        _samples = IndexedSampleList::create(samples);

        auto* allele_list = IndexedAlleleList<pAllele>::create(all_old, _pool);
        pRALikelihoods original = RALikelihoods::create<pReadRecord, pAllele>(_pool, _samples, allele_list, reads);
        fill_with_random_likelihoods(samples, all_old.size(), original);
        auto* marginalized = original->marginalize(new_alleles, allele_map);

        size_t num_alleles = marginalized->number_of_alleles();
        size_t num_samples = samples.size();
        ASSERT_EQ(num_alleles, allele_map.size());
        for (size_t a = 0; a < num_alleles; ++a) {
            const AlleleVector& old_alleles = allele_map.at(marginalized->get_allele(a));
            ASSERT_FALSE(old_alleles.empty());
            for (size_t s = 0; s < num_samples; ++s) {
                auto* old_smaple_likelihoods = original->sample_matrix(s);
                auto* sample_likelihoods = marginalized->sample_matrix(s);
                size_t sample_read_count = sample_likelihoods->evidence_count();
                size_t old_sample_read_count = old_smaple_likelihoods->evidence_count();
                ASSERT_EQ(sample_read_count, old_sample_read_count) << "a=" << a << " s=" << s;
                for (size_t r = 0; r < sample_read_count; ++r) {
                    double old_best_lk = NEGATIVE_INFINITY;
                    for (pAllele old_allele : old_alleles) {
                        size_t a_idx = size_t(original->index_of_allele(old_allele));
                        old_best_lk = std::max(old_smaple_likelihoods->get(a_idx, r), old_best_lk);
                    }
                    ASSERT_EQ(sample_likelihoods->get(a, r), old_best_lk) << "a=" << a << " s=" << s << " r=" << r;
                }
            }
        }

        delete _samples;
    }
}

TEST_F(AlleleLikelihoodsUnitTest, testAddNonRefAllele)
{
    pmr::vector<pmr::vector<pAllele>> all_alleles = allele_sets();
    for (const vector<string>& samples : s_samples) {
        _samples = IndexedSampleList::create(samples);
        for (const pmr::vector<pAllele>& alleles : all_alleles) {
            auto* allele_list = IndexedAlleleList<pAllele>::create(alleles, _pool);
            pmr::map<int32_t, pmr::vector<pReadRecord>> reads = data_set_reads(samples);

            pRALikelihoods original = RALikelihoods::create<pReadRecord, pAllele>(_pool, _samples, allele_list, reads);
            pRALikelihoods result = RALikelihoods::create<pReadRecord, pAllele>(_pool, _samples, allele_list, reads);
            fill_with_random_likelihoods(samples, alleles.size(), original);
            DoubleVector3D original_likelihoods = fill_with_random_likelihoods(samples, alleles.size(), result);

            result->add_non_reference_allele();
            ASSERT_EQ(result->number_of_alleles(), original->number_of_alleles() + 1);
            ASSERT_EQ(result->index_of_allele(StaticAllele::get_instance()->_non_ref_allele.get()),
                      int32_t(result->number_of_alleles() - 1));

            DoubleVector3D new_likelihoods{original_likelihoods.size(), _pool};

            size_t sample_count = samples.size();
            for (size_t s = 0; s < sample_count; ++s) {
                new_likelihoods[s] = original_likelihoods.at(s);
                uint32_t sample_read_count = original->sample_evidence_count(s);
                size_t ordinarynumber_of_alleles = original_likelihoods[s].size();
                new_likelihoods[s].emplace_back(sample_read_count, 0.0);
                for (uint32_t r = 0; r < sample_read_count; ++r) {
                    double best_lk = new_likelihoods[s][0][r];
                    double second_best_lk = NEGATIVE_INFINITY;
                    for (size_t a = 0; a < ordinarynumber_of_alleles; ++a) {
                        double lk = original_likelihoods[s][a][r];
                        if (lk > best_lk) {
                            second_best_lk = best_lk;
                            best_lk = lk;
                        }
                        else if (lk > second_best_lk) {
                            second_best_lk = lk;
                        }
                    }

                    DoubleVector qualifyling_likelihoods{_pool};
                    for (size_t a = 0; a < ordinarynumber_of_alleles; a++) {
                        if (original_likelihoods[s][a][r] >= best_lk) continue;
                        qualifyling_likelihoods.push_back(original_likelihoods[s][a][r]);
                    }

                    double median_likelihood = RHLikelihoods::evaluate(qualifyling_likelihoods, long(qualifyling_likelihoods.size()));
                    double expected_non_ref_lk = !isnan(median_likelihood)        ? median_likelihood
                                                 : ordinarynumber_of_alleles <= 1 ? NAN
                                                                                  : best_lk;
                    new_likelihoods[s][ordinarynumber_of_alleles][r] = expected_non_ref_lk;
                }
            }
            test_likelihood_matrix_queries(samples, result, new_likelihoods);
        }
        delete _samples;
    }
}

pmr::vector<pmr::vector<pAllele>> AlleleLikelihoodsUnitTest::allele_sets()
{
    return {
        {
            {Allele::create_allele("A", true, _pool), Allele::create_allele("T", false, _pool), Allele::create_allele("C", false, _pool)},
            {Allele::create_allele("A", true, _pool)},
            {Allele::create_allele("ATTTA", false, _pool), Allele::create_allele("A", true, _pool)},
            {Allele::create_allele("A", false, _pool), Allele::create_allele("AT", true, _pool)},
            {Allele::create_allele("A", false, _pool), Allele::create_allele("AT", false, _pool)},
        },
        _pool};
}

pmr::map<int32_t, pmr::vector<pReadRecord>> AlleleLikelihoodsUnitTest::data_set_reads(const vector<string>& samples)
{
    std::mt19937 gen{GATK_RANDOM_SEED};
    std::uniform_int_distribution<int32_t> dis(50, 100);
    pmr::map<int32_t, pmr::vector<pReadRecord>> result{_pool};
    for (int32_t i = 0, len = (int32_t)samples.size(); i < len; ++i) {
        result.insert({i, {}});
        int32_t read_count = dis(gen);
        for (int32_t r = 0; r < read_count; ++r) {
            int64_t alignment_start = (r & 1) == 0 ? s_even_read_start : s_odd_read_start;
            string name = "RRR" + std::to_string(i) + "00" + std::to_string(r);
            pReadRecord read = UnitTestUtils::create_artificial_read(_header->header, name.c_str(), 0, alignment_start, 5, _bampool, _pool);
            result.at(i).push_back(read);
        }
    }
    return result;
}

void AlleleLikelihoodsUnitTest::test_sample_queries(const vector<string>& samples, const pmr::map<int32_t, pmr::vector<pReadRecord>>& reads,
                                                    pRALikelihoods result)
{
    pmr::set<int32_t> sample_ids{_pool};
    for (const string& sample : samples) {
        int32_t index_of_sample = result->index_of_sample(sample);
        ASSERT_TRUE(index_of_sample >= 0);
        ASSERT_FALSE(sample_ids.count(index_of_sample));
        sample_ids.insert(index_of_sample);

        const auto& sample_reads = result->sample_evidence((size_t)index_of_sample);
        const auto& expected_sample_read_array = reads.at(index_of_sample);
        ASSERT_EQ(sample_reads.size(), expected_sample_read_array.size());

        int32_t sample_read_count = (int32_t)sample_reads.size();
        for (int32_t r = 0; r < sample_read_count; ++r) {
            // 此处因为并未对数据做更改，仅比较指针
            ASSERT_EQ(sample_reads.at(r), expected_sample_read_array.at(r)) << sample << " " << r;
            int32_t read_index = result->evidence_index((size_t)index_of_sample, sample_reads.at(r));
            ASSERT_EQ(read_index, r) << sample << " " << r;
        }
    }
}

void AlleleLikelihoodsUnitTest::test_allele_queries(const pmr::vector<pAllele>& alleles, pRALikelihoods result)
{
    pmr::set<int32_t> allele_indices{_pool};
    for (pAllele a : alleles) {
        int32_t index_of_allele = result->index_of_allele(a);
        ASSERT_TRUE(index_of_allele >= 0);
        ASSERT_FALSE(allele_indices.count(index_of_allele));
        allele_indices.insert(index_of_allele);
        // 仅比较指针
        ASSERT_EQ(a, alleles.at(index_of_allele));
    }
}

void AlleleLikelihoodsUnitTest::test_likelihood_matrix_queries(const vector<string>& samples, pRALikelihoods result,
                                                               const DoubleVector3D& likelihoods)
{
    for (const string& sample : samples) {
        int32_t index_of_sample = result->index_of_sample(sample);
        uint32_t sample_read_count = result->sample_evidence_count((size_t)index_of_sample);
        size_t number_of_alleles = result->number_of_alleles();

        auto* sample_matrix = result->sample_matrix((size_t)index_of_sample);
        for (size_t a = 0; a < number_of_alleles; ++a) {
            for (size_t r = 0; r < (size_t)sample_read_count; ++r) {
                if (std::isnan(sample_matrix->get(a, r))) {
                    ASSERT_TRUE(std::isnan(likelihoods.at(index_of_sample).at(a).at(r)));
                }
                else {
                    ASSERT_NEAR(sample_matrix->get(a, r), likelihoods.at(index_of_sample).at(a).at(r), s_epsilon);
                }
            }
        }
    }
}

DoubleVector3D AlleleLikelihoodsUnitTest::fill_with_random_likelihoods(const vector<string>& samples, size_t allele_count,
                                                                       pRALikelihoods result)
{
    std::mt19937 gen{GATK_RANDOM_SEED};
    std::uniform_real_distribution<> dis(0.0, 5.0);
    DoubleVector3D ret{_pool};
    ret.resize(samples.size());
    std::for_each(ret.begin(), ret.end(), [&](DoubleVector2D& d2) { d2.resize(allele_count); });

    for (size_t s = 0, s_len = samples.size(); s < s_len; ++s) {
        size_t evidence_count = result->sample_evidence_count(s);
        auto* sample_likelihoods = result->sample_matrix(s);
        for (size_t a = 0; a < allele_count; ++a) {
            DoubleVector& ll = ret.at(s).at(a);
            ll.resize(evidence_count);
            for (size_t r = 0; r < evidence_count; ++r) {
                double value = -std::abs(dis(gen) + 0.0);
                ll.at(r) = value;
                sample_likelihoods->set(a, r, value);
            }
        }
    }
    return ret;
}

void AlleleLikelihoodsUnitTest::check_evidence_to_index_map_is_correct(pRALikelihoods likelihoods)
{
    for (size_t s = 0, slen = likelihoods->number_of_samples(); s < slen; ++s) {
        for (uint32_t r = 0, rln = likelihoods->sample_evidence_count(s); r < rln; ++r) {
            ASSERT_EQ(likelihoods->evidence_index(s, likelihoods->sample_evidence(s).at(r)), r) << "s=" << s << " r=" << r;
        }
    }
}

TupleVector AlleleLikelihoodsUnitTest::marginalization_data_sets()
{
    pmr::vector<pmr::vector<pAllele>> alleles_set = allele_sets();
    TupleVector result{_pool};
    for (const vector<string>& sample : s_samples) {
        for (const pmr::vector<pAllele>& as1 : alleles_set) {
            for (const pmr::vector<pAllele>& as2 : alleles_set) {
                if (as2.size() < as1.size()) {
                    Int32ToReadVectorMap reads = data_set_reads(sample);
                    pmr::map<pAllele, AlleleVector> allele_map = random_allele_map(as1, as2);
                    result.emplace_back(sample, as1, as2, reads, allele_map);
                }
            }
        }
    }
    return result;
}

pmr::map<pAllele, AlleleVector> AlleleLikelihoodsUnitTest::random_allele_map(const AlleleVector& from_alleles,
                                                                             const AlleleVector& to_alleles)
{
    size_t next_to_index = 0;
    pmr::vector<size_t> indexs{_pool};
    std::mt19937 gen{GATK_RANDOM_SEED};
    pmr::map<pAllele, AlleleVector> result{_pool};
    std::for_each(to_alleles.begin(), to_alleles.end(), [&](pAllele a) { result.insert({a, {}}); });
    std::for_each(from_alleles.begin(), from_alleles.end(), [&](pAllele) { indexs.push_back(next_to_index++); });
    std::shuffle(indexs.begin(), indexs.end(), gen);

    next_to_index = 0;
    std::for_each(indexs.begin(), indexs.end(), [&](size_t i) {
        result.at(to_alleles.at(next_to_index)).push_back(from_alleles.at(i));
        next_to_index = (next_to_index + 1) % to_alleles.size();
    });
    return result;
}
