#ifndef ROVACA_HC_MANN_WHITNEY_U_H_
#define ROVACA_HC_MANN_WHITNEY_U_H_
#include "forward.h"
#include "genotype_macors.h"

namespace rovaca
{

typedef class MannWhitneyU MannWhitneyU, *pMannWhitneyU;

/*!
 * Mann-Whitney-U测试（又称为Wilcoxon秩和检验）是一种非参数统计方法，用于检验两个独立样本是否来自具有相同分布的总体。
 * 它不需要数据满足正态分布的假设，因此在数据不满足正态分布时，可以作为t检验的替代方法。
 * 以下是Mann-Whitney U测试的实现步骤：
 *   将两个独立样本合并为一个数据集。
 *   对合并后的数据集进行排序（从小到大），并为每个数据点分配一个秩次（排名）。
 *   将原始秩次分配回各自的样本中。
 *   分别计算两个样本的秩次和，记为R1和R2。
 *   计算Mann-Whitney U统计量。U1 = n1 * n2 + (n1 * (n1 + 1)) / 2 - R1，U2 = n1 * n2 + (n2 * (n2 + 1)) / 2 -
 * R2，其中n1和n2分别为两个样本的大小。 选择较小的U值作为最终的U统计量。 计算U统计量的标准误差（SE）和z得分。SE = sqrt((n1 * n2 * (n1 + n2 +
 * 1)) / 12)，z = (U - n1 * n2 / 2) / SE。
 *   根据z得分和给定的显著性水平（如0.05），查找对应的临界值，以判断是否拒绝原假设（两个样本来自具有相同分布的总体）。
 */
class MannWhitneyU
{
    /* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
private:
    typedef struct Rank
    {
        double value;
        float rank;
        int32_t series;
    } Rank, *pRank;

    /* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
public:
    enum TestType { FIRST_DOMINATES = 0, SECOND_DOMINATES, TWO_SIDED };

    typedef struct Result
    {
        double u, z, p, median_shift;
    } Result, *pResult;

    typedef struct TestStatistic
    {
        double u1, u2, true_u, num_of_ties_transformed;
    } TestStatistic, *pTestStatistic;

    struct RankedData
    {
        std::pmr::vector<pRank> rank;
        std::pmr::vector<size_t> num_of_ties;

        explicit RankedData(pMemoryPool pool)
            : rank(pool)
            , num_of_ties(pool)
        {}
    };

    typedef struct RankedData RankedData, *pRankedData;

    /* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
private:
    pMemoryPool _pool;

    /* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
public:
    static pMannWhitneyU create(pMemoryPool pool) { return new ALLOC_TYPE_IN_POOL(pool, MannWhitneyU) MannWhitneyU{pool}; }

    /*!
     * @brief Constructs a new rank sum test with the given data.
     * @param series1
     * @param series2
     * @param which_side
     * @return
     */
    pResult test(DoubleVector& series1, DoubleVector& series2, TestType which_side);

    pTestStatistic calculate_one_sided_u(DoubleVector& series1, DoubleVector& series2, TestType which_series_dominates);
    pTestStatistic calculate_two_sided_u(DoubleVector& series1, DoubleVector& series2);

    /**
     * Rank both groups together and return a TestStatistic object that includes U1, U2 and number of ties for sigma
     */
    pTestStatistic calculate_u1and_u2(DoubleVector& series1, DoubleVector& series2);

    pRankedData calculate_rank(DoubleVector& series1, DoubleVector& series2);

    double transform_ties(size_t num_of_ranks, const std::pmr::vector<size_t>& num_of_ties);

    /**
     * Calculates the Z score (i.e. standard deviations from the mean) of the rank sum
     * test statistics given input data of lengths n1 and n2 respectively, as well as the number of ties, for normal
     * approximation only.
     */
    static double calculate_z(double u, int32_t n1, int32_t n2, double nties, TestType which_side);

    /**
     * creates histogram of test statistics from a permutation test.
     * @param series1 data from group 1
     * @param series2 data from group 2
     * @param test_stat_u test statistic u from observed data
     * @return p-value based on histogram with u calculated for every possible permutation of group tag.
     */
    double permutation_test(DoubleVector& series1, DoubleVector& series2, double test_stat_u);

    /* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
private:
    explicit MannWhitneyU(pMemoryPool pool)
        : _pool(pool)
    {}

    /**
     * checks to see if the permutations have already been computed before creating them from scratch.
     * @param list_to_permute list of tags in numerical order to be permuted
     * @param num_of_permutations the number of permutations this list will have (n1+n2 choose n1)
     * @return set of all possible permutations for the given list.
     */
    std::pmr::vector<std::pmr::vector<int32_t>> get_permutations(std::pmr::vector<int32_t>& permute, int32_t num_of_permutations);
};

}  // namespace rovaca

#endif  // ROVACA_HC_MANN_WHITNEY_U_H_
