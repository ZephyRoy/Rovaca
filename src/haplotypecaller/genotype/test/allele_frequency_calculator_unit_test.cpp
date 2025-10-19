#include <gtest/gtest.h>

#include "af_calculation_result.h"
#include "allele.h"
#include "allele_frequency_calculator.h"
#include "forward.h"
#include "genotype.h"
#include "genotype_likelihood_calculator.h"
#include "genotypes_context.hpp"
#include "variant.h"

using namespace std;
using namespace rovaca;

static constexpr int32_t s_diploid = 2;
static constexpr int32_t s_triploid = 3;
static constexpr int32_t s_biallelic = 2;
static constexpr int32_t s_triallelic = 3;
static constexpr int32_t s_extremely_confident_pl = 1000;
static constexpr int32_t s_fairly_confident_pl = 20;
static constexpr size_t s_buffer_size = 1024 * 1024 * 10;

class AlleleFrequencyCalculatorUnitTest : public ::testing::Test
{
protected:
    void SetUp() override
    {
        _buffer = new uint8_t[s_buffer_size]{};
        _pool = new std::pmr::monotonic_buffer_resource(_buffer, s_buffer_size, std::pmr::null_memory_resource());
        A = Allele::create_allele('A', 1);
        B = Allele::create_allele('C', 0);
        C = Allele::create_allele('G', 0);
        sample_name_counter = 0;
    }

    void TearDown() override
    {
        delete _pool;
        delete[] _buffer;
    }

    pGenotype make_genotype(int32_t ploidy, Int32Vector&& pls);
    pVariant make_vc(const AlleleVector& alleles, const std::pmr::vector<pGenotype>& genotypes);
    Int32Vector pls_for_obvious_call(int32_t ploidy, int32_t num_alleles, const Int32Vector& alleles, int32_t pl);
    pGenotype genotype_with_obvious_call(int32_t ploidy, int32_t num_alleles, const Int32Vector& alleles, int32_t pl);

    uint8_t* _buffer{};
    pMemoryPool _pool{};

    pAllele A{}, B{}, C{};
    int32_t sample_name_counter{};
};

TEST_F(AlleleFrequencyCalculatorUnitTest, testMLECounts)
{
    AlleleVector alleles{{A, B, C}, _pool};
    pGenotype aa = genotype_with_obvious_call(s_diploid, s_triallelic, {{0, 2}, _pool}, s_fairly_confident_pl);
    pGenotype bb = genotype_with_obvious_call(s_diploid, s_triallelic, {{1, 2}, _pool}, s_fairly_confident_pl);
    pGenotype ab = genotype_with_obvious_call(s_diploid, s_triallelic, {{0, 1, 1, 1}, _pool}, s_fairly_confident_pl);
    pGenotype ac = genotype_with_obvious_call(s_diploid, s_triallelic, {{0, 1, 2, 1}, _pool}, s_fairly_confident_pl);

    pmr::vector<pmr::vector<int32_t>> target{{{2, 0}, {1, 0}, {2, 0}, {2, 0}, {2, 0}, {1, 1}}, _pool};
    pmr::vector<pmr::vector<pGenotype>> genotypes{{{aa, bb}, {aa, ab}, {ab, ab}, {aa, aa, bb}, {aa, ab, ab}, {aa, ab, ac}}, _pool};

    pVariant vc;
    pAlleleFrequencyCalculator calculator = AlleleFrequencyCalculator::create(1, 1, 1, 2);
    for (size_t i = 0, len = target.size(); i < len; ++i) {
        const pmr::vector<pGenotype>& gs = genotypes.at(i);
        vc = make_vc(alleles, gs);
        const Int32Vector& actual = calculator->calculate(vc)->get_allele_counts_of_mle();
        const Int32Vector& expected = target.at(i);
        ASSERT_EQ(expected, actual) << "i = " << i;
    }
    delete calculator;
}

pGenotype AlleleFrequencyCalculatorUnitTest::genotype_with_obvious_call(int32_t ploidy, int32_t num_alleles, const Int32Vector& alleles,
                                                                        int32_t pl)
{
    return make_genotype(ploidy, pls_for_obvious_call(ploidy, num_alleles, alleles, pl));
}

Int32Vector AlleleFrequencyCalculatorUnitTest::pls_for_obvious_call(int32_t ploidy, int32_t num_alleles, const Int32Vector& alleles,
                                                                    int32_t pl)
{
    pGenotypeLikelihoodCalculator calculator = GenotypeLikelihoodCalculator::create(ploidy, num_alleles, _pool);
    Int32Vector result{(size_t)calculator->genotype_count(), pl, _pool};
    result[calculator->allele_counts_to_index(alleles)] = 0;
    return result;
}

pGenotype AlleleFrequencyCalculatorUnitTest::make_genotype(int32_t ploidy, Int32Vector&& pls)
{
    pGenotype result = Genotype::create(_pool);
    AlleleVector alleles{(size_t)ploidy, StaticAllele::get_instance()->_no_call.get(), _pool};
    result->set_id(sample_name_counter++);
    result->set_alleles(std::move(alleles));
    result->set_pl(std::forward<Int32Vector>(pls));
    return result;
}

pVariant AlleleFrequencyCalculatorUnitTest::make_vc(const AlleleVector& alleles, const std::pmr::vector<pGenotype>& genotypes)
{
    pVariant result = Variant::create(_pool);
    pGenotypesContext gc = GenotypesContext::create(_pool);
    std::for_each(genotypes.begin(), genotypes.end(), [&](pGenotype g) { gc->add(g); });
    result->set_tid(1);
    result->set_alleles(alleles);
    result->set_genotype(gc);
    return result;
}
