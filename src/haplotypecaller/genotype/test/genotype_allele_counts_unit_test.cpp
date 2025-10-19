#include <gtest/gtest.h>

#include "allele.h"
#include "forward.h"
#include "genotype_allele_counts.h"
#include "genotype_struct.h"

using namespace std;
using namespace rovaca;

static constexpr size_t s_buffer_size = 1024 * 1024 * 100;

#define MAXIMUM_ALLELE_INDEX (10)

class GenotypeAlleleCountsUnitTest : public ::testing::Test
{
protected:
    void SetUp() override
    {
        _buffer = new uint8_t[s_buffer_size]{};
        _pool = new std::pmr::monotonic_buffer_resource(_buffer, s_buffer_size, std::pmr::null_memory_resource());
    }

    void TearDown() override
    {
        delete _pool;
        delete[] _buffer;
    }

    AlleleVector getTestAlleles()
    {
        AlleleVector t(_pool);
        pBases bases;
        t.push_back(Allele::create_allele('A', 1));
        for (uint32_t i = 1; i <= 50; ++i) {
            bases = new ALLOC_FLEXIBLE_IN_POOL(_pool, Bases, i + 1, uint8_t) Bases{i + 1};
            for (uint32_t j = 0; j < i + 1; ++j) {
                bases->data[j] = 'A';
            }
            t.push_back(Allele::create_allele(bases, i == 0, _pool));
        }
        return t;
    }

    uint8_t* _buffer{};
    pMemoryPool _pool{};
};

TEST_F(GenotypeAlleleCountsUnitTest, testFirst)
{
    int32_t ploidy = 2;
    pGenotypeAlleleCounts subject = GenotypeAlleleCounts::first(ploidy, _pool);
    ASSERT_EQ(subject->ploidy(), ploidy);
    ASSERT_EQ(subject->index(), 0);
    ASSERT_EQ(subject->distinct_allele_count(), 1);
    ASSERT_EQ(subject->log10combination_count(), 0.0);
    ASSERT_EQ(subject->allele_index_at(0), 0);
    ASSERT_EQ(subject->allele_count_at(0), ploidy);
    ASSERT_EQ(subject->allele_count_for(0), ploidy);
    ASSERT_EQ(subject->allele_rank_for(0), 0);
    ASSERT_EQ(subject->allele_rank_for(1), -2);
    ASSERT_EQ(subject->maximum_allele_index(), 0);
    ASSERT_EQ(subject->minimum_allele_index(), 0);
    ASSERT_TRUE(subject->contains_allele(0));
    ASSERT_FALSE(subject->contains_allele(1));

    AlleleVector alleles = getTestAlleles();
    AlleleVector allele_list = subject->as_allele_list(alleles);
    for (int32_t i = 0; i < ploidy; ++i) {
        ASSERT_TRUE(allele_list.at(i)->equals(*alleles.at(0)));
    }
}

TEST_F(GenotypeAlleleCountsUnitTest, testNext_PloidyTwo)
{
    int32_t ploidy = 2;
    AlleleVector testAlleles = getTestAlleles();

    pGenotypeAlleleCounts subject = GenotypeAlleleCounts::first(ploidy, _pool);
    pGenotypeAlleleCounts next;

    while (true) {
        next = subject->next(_pool);
        if (next->contains_allele(MAXIMUM_ALLELE_INDEX + 1)) {
            break;
        }

        // test log10CombinationCount
        ASSERT_DOUBLE_EQ(next->log10combination_count(), next->distinct_allele_count() == 2 ? std::log10(2) : 0.0);

        // test forEach
        pmr::vector<int32_t> alleleCountsAsList(_pool);
        alleleCountsAsList.reserve(next->distinct_allele_count() * 2);
        pmr::set<int32_t> absentAlleles(_pool);
        next->for_each_allele_index_and_count([&](int32_t alleleIndex, int32_t alleleCount) {
            alleleCountsAsList.push_back(alleleIndex);
            alleleCountsAsList.push_back(alleleCount);
        });
        next->for_each_absent_allele_index([&](int32_t i) { absentAlleles.insert(i); }, MAXIMUM_ALLELE_INDEX + 1);

        Int32Vector actualAlleleCounts(next->distinct_allele_count() * 2, _pool);
        next->copy_allele_counts(0, actualAlleleCounts);
        ASSERT_EQ(alleleCountsAsList, actualAlleleCounts);

        ASSERT_EQ(absentAlleles.size(), MAXIMUM_ALLELE_INDEX + 1 - next->distinct_allele_count());
        next->for_each_allele_index_and_count(
            [&](int32_t index, [[maybe_unused]] int32_t count) { ASSERT_FALSE(absentAlleles.count(index)); });

        if (1 == subject->distinct_allele_count()) {
            ASSERT_EQ(next->maximum_allele_index(), subject->maximum_allele_index() + 1);
            ASSERT_EQ(next->distinct_allele_count(), 2);
            ASSERT_EQ(next->minimum_allele_index(), 0);
        }
        else {
            ASSERT_EQ(next->maximum_allele_index(), subject->maximum_allele_index());
            ASSERT_EQ(next->minimum_allele_index(), subject->allele_count_at(0) > 1    ? 0
                                                    : subject->allele_count_at(0) == 1 ? subject->minimum_allele_index() + 1
                                                                                       : subject->minimum_allele_index());
        }

        // checking on 0's new count and subject->min_allele + 1 alleles.
        ASSERT_EQ(next->allele_count_for(0), subject->allele_count_for(subject->minimum_allele_index()) - 1);
        ASSERT_EQ(next->allele_count_for(subject->minimum_allele_index() + 1),
                  subject->allele_count_for(subject->minimum_allele_index() + 1) + 1);

        // checks subject->min_allele count
        ASSERT_EQ(next->allele_count_for(subject->minimum_allele_index()),
                  subject->minimum_allele_index() == 0 ? subject->allele_count_at(0) - 1 : 0);

        int32_t totalCountSum = 0;
        Int32Vector expectedAlleleCountsByIndex(std::max(MAXIMUM_ALLELE_INDEX, next->maximum_allele_index()) + 1, _pool);
        for (int32_t i = 0, len = next->distinct_allele_count(); i < len; ++i) {
            int32_t count = next->allele_count_at(i);
            int32_t index = next->allele_index_at(i);
            expectedAlleleCountsByIndex[index] = count;

            // Check consistency of alleleCountAt(x) and alleleCountFor(alleleIndexAt(x))
            ASSERT_EQ(next->allele_count_for(index), count);
            totalCountSum += count;

            // check on counts of, in theory, unaffected allele counts.
            if (index > subject->minimum_allele_index() + 1) {
                ASSERT_EQ(next->allele_count_for(index), subject->allele_count_for(index));
            }
        }

        ASSERT_EQ(totalCountSum, ploidy);
        ASSERT_EQ(next->index(), subject->index() + 1);
        ASSERT_EQ(next->ploidy(), ploidy);

        AlleleVector expectedList(_pool);
        expectedList.reserve(ploidy);
        for (int i = 0; i < next->distinct_allele_count(); i++) {
            for (int j = 0; j < next->allele_count_at(i); j++) {
                expectedList.push_back(testAlleles.at(next->allele_index_at(i)));
            }
        }
        ASSERT_EQ(next->as_allele_list(testAlleles), expectedList);
        subject = next;
    }
}

TEST_F(GenotypeAlleleCountsUnitTest, testIncrease_PloidyTwo)
{
    int32_t ploidy = 2;
    pGenotypeAlleleCounts next = GenotypeAlleleCounts::first(ploidy, _pool);

    while (!next->contains_allele(MAXIMUM_ALLELE_INDEX + 1)) {
        pGenotypeAlleleCounts current = next->copy(_pool);
        next->increase();

        if (current->distinct_allele_count() == 1) {
            ASSERT_EQ(next->maximum_allele_index(), current->maximum_allele_index() + 1);
            ASSERT_EQ(next->distinct_allele_count(), 2);
            ASSERT_EQ(next->minimum_allele_index(), 0);
        }
        else {
            ASSERT_EQ(next->maximum_allele_index(), current->maximum_allele_index()) << "index = " << current->index() << endl;
            ASSERT_EQ(next->minimum_allele_index(), current->allele_count_at(0) > 1    ? 0
                                                    : current->allele_count_at(0) == 1 ? current->minimum_allele_index() + 1
                                                                                       : current->minimum_allele_index());
        }

        // checking on 0's new count and current.min_allele + 1 alleles.
        ASSERT_EQ(next->allele_count_for(0), current->allele_count_for(current->minimum_allele_index()) - 1);
        ASSERT_EQ(next->allele_count_for(current->minimum_allele_index() + 1),
                  current->allele_count_for(current->minimum_allele_index() + 1) + 1);

        // checks current.min_allele count
        ASSERT_EQ(next->allele_count_for(current->minimum_allele_index()),
                  current->minimum_allele_index() == 0 ? current->allele_count_at(0) - 1 : 0);

        int32_t totalCountSum = 0;
        Int32Vector expectedAlleleCountsByIndex(std::max(MAXIMUM_ALLELE_INDEX, next->maximum_allele_index()) + 1, _pool);
        for (int32_t i = 0, len = next->distinct_allele_count(); i < len; ++i) {
            int32_t count = next->allele_count_at(i);
            int32_t index = next->allele_index_at(i);
            expectedAlleleCountsByIndex[index] = count;

            // Check consistency of alleleCountAt(x) and alleleCountFor(alleleIndexAt(x))
            ASSERT_EQ(next->allele_count_for(index), count);
            totalCountSum += count;

            // check on counts of, in theory, unaffected allele counts.
            if (index > current->minimum_allele_index() + 1) {
                ASSERT_EQ(next->allele_count_for(index), current->allele_count_for(index));
            }
        }
        ASSERT_EQ(totalCountSum, ploidy);
        ASSERT_EQ(next->index(), current->index() + 1);
        ASSERT_EQ(next->ploidy(), ploidy);
    }
}