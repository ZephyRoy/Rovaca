#include "mann_whitney_u.h"

#include <algorithm>
#include <numeric>

#include "rovaca_logger.h"
#include "math_utils.h"

namespace rovaca
{

static constexpr size_t s_minimum_normal_n = 10;
static const double s_sqrt2 = std::sqrt(2.0);
static constexpr double s_double_precision = 1.0E-6;

static inline bool cigar_op_is_equal(double x, double y) { return std::abs(x - y) < s_double_precision; }

static inline double cumulative_probability(double x)
{
    if (std::abs(x) > 40.0) {
        return x < 0.0 ? 0.0 : 1.0;
    }
    return 0.5 * (1 + std::erf(x / s_sqrt2));
}

static inline double inverse_cumulative_probability(double p)
{
    CHECK_CONDITION_EXIT(p < 0.0 || p > 1.0, "p < 0.0 || p > 1.0");
    return s_sqrt2 * math_utils::erf_inv(2 * p - 1);
}

static inline double median(const DoubleVector& data)
{
    size_t len = data.size();
    size_t mid = len / 2;
    if ((len & 1) == 0) {
        return (data.at(mid) + data.at(mid - 1)) / 2.0;
    }
    else {
        return data.at(mid);
    }
}

MannWhitneyU::pResult MannWhitneyU::test(DoubleVector& series1, DoubleVector& series2, MannWhitneyU::TestType which_side)
{
    size_t n1 = series1.size();
    size_t n2 = series2.size();

    // If one of the groups is empty we return NaN
    if (n1 == 0 || n2 == 0) {
        return nullptr;
    }

    pTestStatistic result;
    if (ROVACA_UNLIKELY(which_side == TWO_SIDED)) {
        result = calculate_two_sided_u(series1, series2);
    }
    else {
        result = calculate_one_sided_u(series1, series2, which_side);
    }
    double u = result->true_u;
    double nties = result->num_of_ties_transformed;
    double z, p;
    if (n1 >= s_minimum_normal_n || n2 >= s_minimum_normal_n) {
        z = calculate_z(u, (int32_t)n1, (int32_t)n2, nties, which_side);
        p = 2.0 * cumulative_probability(z);
        if (which_side == TestType::TWO_SIDED) {
            p = p / 2;
        }
    }
    else {
        p = permutation_test(series1, series2, u);
        z = inverse_cumulative_probability(p);
    }

    double abs_median_diff = std::abs(median(series1) - median(series2));
    return new ALLOC_TYPE_IN_POOL(_pool, Result) Result{u, z, p, abs_median_diff};
}

MannWhitneyU::pTestStatistic MannWhitneyU::calculate_one_sided_u(DoubleVector& series1, DoubleVector& series2,
                                                                 MannWhitneyU::TestType which_series_dominates)
{
    pTestStatistic stat = calculate_u1and_u2(series1, series2);
    double u = which_series_dominates == FIRST_DOMINATES ? stat->u1 : stat->u2;
    return new ALLOC_TYPE_IN_POOL(_pool, TestStatistic) TestStatistic{NAN, NAN, u, stat->num_of_ties_transformed};
}

MannWhitneyU::pTestStatistic MannWhitneyU::calculate_two_sided_u(DoubleVector& series1, DoubleVector& series2)
{
    pTestStatistic stat = calculate_u1and_u2(series1, series2);
    double u = std::min(stat->u1, stat->u2);
    return new ALLOC_TYPE_IN_POOL(_pool, TestStatistic) TestStatistic{NAN, NAN, u, stat->num_of_ties_transformed};
}

MannWhitneyU::pTestStatistic MannWhitneyU::calculate_u1and_u2(DoubleVector& series1, DoubleVector& series2)
{
    pRankedData ranked = calculate_rank(series1, series2);
    const std::pmr::vector<pRank>& ranks = ranked->rank;
    const std::pmr::vector<size_t>& num_of_ties = ranked->num_of_ties;
    size_t length_of_ranks = ranks.size();

    double num_of_ties_for_sigma = transform_ties(length_of_ranks, num_of_ties);

    // Calculate R1 and R2 and U.
    float r1 = 0, r2 = 0;
    for (pRank r : ranks) {
        r1 += r->series == 1 ? r->rank : 0.0f;
        r2 += r->series == 1 ? 0.0f : r->rank;
    }

    double n1 = double(series1.size());
    double n2 = double(series2.size());
    double u1 = r1 - ((n1 * (n1 + 1)) / 2);
    double u2 = r2 - ((n2 * (n2 + 1)) / 2);

    return new ALLOC_TYPE_IN_POOL(_pool, TestStatistic) TestStatistic{u1, u2, NAN, num_of_ties_for_sigma};
}

MannWhitneyU::pRankedData MannWhitneyU::calculate_rank(DoubleVector& series1, DoubleVector& series2)
{
    std::sort(series1.begin(), series1.end());
    std::sort(series2.begin(), series2.end());

    pRankedData result = new ALLOC_TYPE_IN_POOL(_pool, RankedData) RankedData{_pool};

    size_t s1_size = series1.size(), s2_size = series2.size();
    size_t i = 0, j = 0, r = 0, len = s1_size + s2_size;
    result->rank.resize(len);
    while (r < len) {
        if (i >= s1_size) {
            result->rank[r++] = new ALLOC_TYPE_IN_POOL(_pool, Rank) Rank{series2[j++], (float)r + 1.0f, 2};
        }
        else if (j >= s2_size) {
            result->rank[r++] = new ALLOC_TYPE_IN_POOL(_pool, Rank) Rank{series1[i++], (float)r + 1.0f, 1};
        }
        else if (series1[i] <= series2[j]) {
            result->rank[r++] = new ALLOC_TYPE_IN_POOL(_pool, Rank) Rank{series1[i++], (float)r + 1.0f, 1};
        }
        else {
            result->rank[r++] = new ALLOC_TYPE_IN_POOL(_pool, Rank) Rank{series2[j++], (float)r + 1.0f, 2};
        }
    }

    float rank;
    size_t count;
    auto& ranks = result->rank;
    for (i = 0; i < len;) {
        rank = ranks.at(i)->rank;
        count = 1;
        for (j = i + 1; j < len && cigar_op_is_equal(ranks.at(i)->value, ranks.at(j)->value); ++j) {
            rank += ranks.at(j)->rank;
            ++count;
        }

        if (count > 1) {
            rank /= count;
            for (r = i; r < i + count; ++r) {
                ranks[r]->rank = rank;
            }
            result->num_of_ties.push_back(count);
        }

        // skip forward the right number of items
        i += count;
    }

    return result;
}

double MannWhitneyU::transform_ties(size_t num_of_ranks, const std::pmr::vector<size_t>& num_of_ties)
{
    DoubleVector transformed_ties(_pool);
    std::for_each(num_of_ties.begin(), num_of_ties.end(), [&](size_t count) {
        if (count != num_of_ranks) {
            transformed_ties.push_back(std::pow(count, 3) - (double)count);
        }
    });
    return std::accumulate(transformed_ties.begin(), transformed_ties.end(), 0.0);
}

double MannWhitneyU::calculate_z(double u, int32_t n1, int32_t n2, double nties, MannWhitneyU::TestType which_side)
{
    double m = (n1 * n2) / 2.0;

    // Adds a continuity correction
    double correction;
    if (which_side == TestType::TWO_SIDED) {
        correction = (u - m) >= 0 ? .5 : -.5;
    }
    else {
        correction = which_side == TestType::FIRST_DOMINATES ? -.5 : .5;
    }

    // If all the data is tied, the number of ties for sigma is set to 0. In order to get a p-value of .5 we need to remove the continuity
    // correction.
    if (nties == 0) {
        correction = 0;
    }
    double sigma = std::sqrt((n1 * n2 / 12.0) * ((n1 + n2 + 1) - nties / ((n1 + n2) * (n1 + n2 - 1))));
    return (u - m - correction) / sigma;
}

double MannWhitneyU::permutation_test(DoubleVector& series1, DoubleVector& series2, double test_stat_u)
{
    std::pmr::map<int32_t, int32_t> histo(_pool);
    size_t n1 = series1.size();
    size_t n2 = series2.size();

    pRankedData ranked_groups = calculate_rank(series1, series2);
    std::pmr::vector<pRank>& ranks = ranked_groups->rank;

    std::pmr::vector<int32_t> first_permutation(n1 + n2, 1, _pool);
    std::fill_n(first_permutation.begin(), n1, 0);

    int32_t num_of_perms = (int32_t)math_utils::binomial_coefficient(int32_t(n1 + n2), (int32_t)n2);
    std::pmr::vector<std::pmr::vector<int32_t>> all_permutations = get_permutations(first_permutation, num_of_perms);

    DoubleVector new_series1(n1, _pool);
    DoubleVector new_series2(n2, _pool);
    size_t series1end, series2end;
    for (const std::pmr::vector<int32_t>& curr_perm : all_permutations) {
        series1end = series2end = 0;
        for (size_t i = 0, len = curr_perm.size(); i < len; i++) {
            int32_t grouping = curr_perm.at(i);
            if (grouping == 0) {
                new_series1[series1end] = ranks[i]->rank;
                series1end++;
            }
            else {
                new_series2[series2end] = ranks[i]->rank;
                series2end++;
            }
        }
        assert(series1end == n1);
        assert(series2end == n2);

        double new_u = std::accumulate(new_series1.begin(), new_series1.end(), 0.0) - (((double)n1 * ((double)n1 + 1.0)) / 2.0);
        int32_t key = (int32_t)std::round(2 * new_u);
        if (!histo.count(key)) {
            histo.insert({key, 1});
        }
        else {
            ++histo.at(key);
        }
    }

    /**
     * in order to deal with edge cases where the observed value is also the most extreme value, we are taking half
     * of the count in the observed bin plus everything more extreme (in the first_dominates case the smaller bins)
     * and dividing by the total count of everything in the histogram. just using get_cumulative_distribution() gives
     * a p-value of 1 in the most extreme case which doesn't result in a usable z-score.
     */
    int32_t count = 0;
    int32_t rr = (int32_t)std::round(2 * test_stat_u);
    double sum_of_all_smaller_bins = histo.at(rr) / 2.0;
    std::for_each(histo.begin(), histo.end(), [&](const std::pair<int32_t, int32_t>& tup) {
        count += tup.second;
        if (tup.first < rr) {
            sum_of_all_smaller_bins += tup.second;
        }
    });

    return sum_of_all_smaller_bins / count;
}

std::pmr::vector<std::pmr::vector<int32_t>> MannWhitneyU::get_permutations(std::pmr::vector<int32_t>& temp, int32_t num_of_permutations)
{
    std::pmr::vector<std::pmr::vector<int32_t>> result(_pool);
    result.reserve(num_of_permutations);
    result.push_back(temp);
    int32_t len = (int32_t)temp.size();
    int32_t i, k, l, begin, end;
    while (true) {
        k = -1;
        for (i = len - 2; i >= 0; i--) {
            if (temp[i] < temp[i + 1]) {
                k = i;
                break;
            }
        }

        if (k == -1) {
            break;
        }

        l = -1;
        for (i = len - 1; i >= k + 1; i--) {
            if (temp[k] < temp[i]) {
                l = i;
                break;
            }
        }

        std::swap(temp.at(k), temp.at(l));

        end = len - 1;
        for (begin = k + 1; begin < end; begin++) {
            std::swap(temp.at(begin), temp.at(end));
            end--;
        }
        result.push_back(temp);
    }
    return result;
}

}  // namespace rovaca