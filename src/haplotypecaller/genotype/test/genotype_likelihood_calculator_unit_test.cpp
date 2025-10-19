#include <gtest/gtest.h>

#include "bam_data_pool.hpp"
#include "forward.h"
#include "genotype_allele_counts.h"
#include "genotype_allele_counts_manger.hpp"
#include "genotype_likelihood_calculator.h"
#include "genotype_likelihoods.h"
#include "math_utils.h"
#include "unit_test_utils.hpp"

using namespace std;
using namespace rovaca;
using namespace UnitTestUtils;

static constexpr size_t s_buffer_size = 1024 * 1024 * 10;

vector<int32_t> s_maximum_allele_arr{1, 2, 5, 6};

vector<vector<int32_t>> s_read_counts{
    {10, 100, 50},
    {0, 100, 10, 1, 50},
    {1, 2, 3, 4, 20},
    {10, 0},
};

/*!
 * @note 仅测试 2 倍体
 */
class GenotypeLikelihoodCalculatorUnitTest : public ::testing::Test
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

    static int32_t calculate_genotype_count(int32_t allele_count) { return ((allele_count) * (allele_count + 1)) >> 1; }

    uint8_t* _buffer{};
    pMemoryPool _pool{};
    pBamDataPool _bampool{};
};

TEST_F(GenotypeLikelihoodCalculatorUnitTest, testPloidyAndMaximumAllele)
{
    int32_t index, i, j, len, ploidy = 2;
    pGenotypeAlleleCounts gac;
    pGenotypeLikelihoodCalculator calculator;

    Int32Vector allele_array(ploidy, _pool);

    for (int32_t allele_count : s_maximum_allele_arr) {
        calculator = GenotypeLikelihoodCalculator::create(ploidy, allele_count, _pool);

        ASSERT_EQ(calculator->ploidy(), ploidy);
        ASSERT_EQ(calculator->allele_count(), allele_count);
        ASSERT_EQ(calculator->genotype_count(), calculate_genotype_count(allele_count));
        int32_t test_genotype_count = std::min(30000, calculator->genotype_count());
        for (i = 0; i < test_genotype_count; ++i) {
            gac = calculator->genotype_allele_counts_at(i);
            index = 0;
            std::fill_n(allele_array.begin(), ploidy, 0);
            for (j = 0, len = gac->distinct_allele_count(); j < len; ++j) {
                std::fill_n(allele_array.begin() + (long)index, gac->allele_count_at(j), gac->allele_index_at(j));
                index += gac->allele_count_at(j);
            }
            Int32Vector allele_count_array(gac->distinct_allele_count() << 1, _pool);
            gac->copy_allele_counts(0, allele_count_array);
            ASSERT_EQ(index, ploidy);
            ASSERT_EQ(calculator->alleles_to_index(allele_array), i) << "allele_count=" << allele_count << " i=" << i;
            ASSERT_EQ(calculator->allele_counts_to_index(allele_count_array), i);
        }
    }
}

TEST_F(GenotypeLikelihoodCalculatorUnitTest, testLikelihoodCalculation)
{
    int32_t ploidy = 2;
    pSamHeader h = create_artificial_sam_header(10, 0, 1000);
    InterfaceSampleList* sl;
    pRALikelihoods read_likelihoods;
    pGenotypeLikelihoodCalculator calculator;

    for (int32_t allele_count : s_maximum_allele_arr) {
        for (const vector<int32_t>& read_count : s_read_counts) {
            sl = sample_list((int32_t)read_count.size());
            read_likelihoods = UnitTestUtils::read_likelihoods(h, allele_count, read_count, sl, _pool, _bampool);
            calculator = GenotypeLikelihoodCalculator::create(ploidy, allele_count, _pool);

            int32_t genotype_count = calculator->genotype_count();
            int32_t test_genotype_count = std::min(30000, genotype_count);
            int32_t sample_count = (int32_t)read_count.size();

            for (int32_t s = 0; s < sample_count; ++s) {
                auto* sample_likelihoods = read_likelihoods->sample_matrix((size_t)s);
                pGenotypeLikelihoods genotype_likelihoods = calculator->genotype_likelihoods(sample_likelihoods);
                ASSERT_TRUE(genotype_likelihoods != nullptr);
                const DoubleVector& genotype_likelihoods_doubles = genotype_likelihoods->_log10likelihoods;
                ASSERT_EQ(genotype_likelihoods_doubles.size(), genotype_count);

                for (int32_t i = 0; i < test_genotype_count; ++i) {
                    pGenotypeAlleleCounts gac = calculator->genotype_allele_counts_at(i);
                    size_t evidence_count = sample_likelihoods->evidence_count();
                    DoubleVector read_genotype_likelihoods{evidence_count, _pool};
                    for (size_t r = 0; r < evidence_count; ++r) {
                        int32_t distinct_allele_count = gac->distinct_allele_count();
                        DoubleVector compoments{(size_t)distinct_allele_count, _pool};
                        for (int ar = 0; ar < distinct_allele_count; ar++) {
                            int32_t a = gac->allele_index_at(ar);
                            int32_t a_count = gac->allele_count_at(ar);
                            double read_lk = sample_likelihoods->get(a, r);
                            compoments[ar] = read_lk + std::log10(a_count);
                        }
                        read_genotype_likelihoods[r] = MathUtils::approximate_log10sum_log10(compoments) - std::log10(ploidy);
                    }
                    double genotype_likelihood = std::accumulate(read_genotype_likelihoods.begin(), read_genotype_likelihoods.end(), 0.0);

                    ASSERT_NEAR(genotype_likelihood, genotype_likelihoods_doubles.at(i), 0.0001);
                }
            }

            delete sl;
        }
    }

    sam_header_destroy(h);
}

TEST_F(GenotypeLikelihoodCalculatorUnitTest, testGenotypeIndexMap)
{
    // final int old_allele_count, final int new_allele_count
    int32_t ploidy = 2, rand;
    std::mt19937 gen{1024};
    for (int32_t old_allele_count : s_maximum_allele_arr) {
        for (int32_t new_allele_count = 0; new_allele_count < old_allele_count * 2; ++new_allele_count) {
            std::uniform_int_distribution<> rnd(0, old_allele_count - 1);
            pmr::map<int32_t, pmr::set<int32_t>> reverse_map{_pool};
            Int32Vector allele_map{(size_t)new_allele_count, _pool};
            for (int32_t n = 0; n < new_allele_count; ++n) {
                rand = rnd(gen);
                allele_map[n] = rand;
                if (!reverse_map.count(rand)) {
                    reverse_map.insert({rand, {}});
                }
                reverse_map.at(rand).insert(n);
            }

            int32_t max_allele_count = std::max(old_allele_count, new_allele_count);
            pGenotypeLikelihoodCalculator calculator = GenotypeLikelihoodCalculator::create(ploidy, max_allele_count, _pool);

            Int32Vector genotype_index_map = calculator->genotype_index_map(allele_map);

            ASSERT_EQ(genotype_index_map.size(),
                      GenotypeAlleleCountsManger::get_instance()->allele_first_genotype_offset_by_ploidy(ploidy).at(new_allele_count));

            pGenotypeLikelihoodCalculator old_calculator = GenotypeLikelihoodCalculator::create(ploidy, old_allele_count, _pool);
            pGenotypeLikelihoodCalculator new_calculator = GenotypeLikelihoodCalculator::create(ploidy, new_allele_count, _pool);

            for (int32_t i = 0, len = (int32_t)genotype_index_map.size(); i < len; ++i) {
                pGenotypeAlleleCounts old_counts = old_calculator->genotype_allele_counts_at(genotype_index_map[i]);
                pGenotypeAlleleCounts new_counts = new_calculator->genotype_allele_counts_at(i);

                Int32Vector reverse_counts{(size_t)old_allele_count, _pool};
                for (int32_t j = 0, j_len = new_counts->distinct_allele_count(); j < j_len; ++j) {
                    int32_t new_index = new_counts->allele_index_at(j);
                    int32_t new_repeats = new_counts->allele_count_at(j);
                    int32_t expected_old_index = allele_map[new_index];
                    int32_t old_index_rank = old_counts->allele_rank_for(expected_old_index);
                    ASSERT_NE(-1, old_index_rank);
                    int32_t old_index = old_counts->allele_index_at(old_index_rank);
                    int32_t old_repeats = old_counts->allele_count_at(old_index_rank);
                    ASSERT_EQ(old_index, expected_old_index);
                    ASSERT_TRUE(old_repeats >= new_repeats);
                    reverse_counts[old_index] += new_repeats;
                }
                for (int j = 0; j < old_allele_count; j++) {
                    ASSERT_EQ(old_counts->allele_count_for(j), reverse_counts[j]);
                }
            }
        }
    }
}