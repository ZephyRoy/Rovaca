#include "math_utils.h"

#include <gtest/gtest.h>

#include <cstdint>
#include <memory_resource>
#include <vector>

using namespace rovaca;
using namespace math_utils;

static constexpr size_t s_max_bufer_size = 1024 * 1024 * 20;

typedef std::pmr::memory_resource MemoryPool, *pMemoryPool;

class MathUtilsUnitTest : public ::testing::Test
{
protected:
    void SetUp() override
    {
        _buffer = new uint8_t[s_max_bufer_size]{};
        _pool = new std::pmr::monotonic_buffer_resource(_buffer, s_max_bufer_size, std::pmr::null_memory_resource());
    }

    void TearDown() override
    {
        delete _pool;
        delete[] _buffer;
    }

    static constexpr double _precision = 1e-6;

    uint8_t* _buffer{};
    pMemoryPool _pool{};
};

TEST_F(MathUtilsUnitTest, testLogLog10conversions)
{
    std::pmr::vector<double> values({std::exp(1), 10.0, 0.5, 1, 1.5, 100.0}, _pool);
    double log_x, log10_x;
    for (double d : values) {
        log_x = std::log(d);
        log10_x = std::log10(d);
        ASSERT_NEAR(log10_to_log(log10_x), log_x, _precision) << "log10ToLog";
        ASSERT_NEAR(log_to_log10(log_x), log10_x, _precision) << "logToLog10";
        ASSERT_NEAR(log_to_log10(log10_to_log(log10_x)), log10_x, _precision) << "log10->log->log10";
        ASSERT_NEAR(log10_to_log(log_to_log10(log_x)), log_x, _precision) << "log->log10->log";
    }
}

TEST_F(MathUtilsUnitTest, testLog10Gamma)
{
    ASSERT_NEAR(log10gamma(4.0), 0.7781513, _precision);
    ASSERT_NEAR(log10gamma(10), 5.559763, _precision);
    ASSERT_NEAR(log10gamma(10654), 38280.532152137, _precision);
}

TEST_F(MathUtilsUnitTest, testLog1mexp)
{
    std::pmr::vector<double> keys({std::exp(1), 10, 1, 0, -1, -3, -10, -30, -100, -300, -1000, -3000}, _pool);
    double actual, expected;
    for (double x : keys) {
        actual = math_utils::log1mexp(x);
        expected = std::log(1 - std::exp(x));
        if (std::isnan(expected)) {
            ASSERT_TRUE(std::isnan(actual));
        }
        else if (std::isinf(actual) || std::isinf(expected)) {
            ASSERT_TRUE(std::isinf(actual) || std::isinf(expected));
        }
        else {
            ASSERT_NEAR(actual, expected, 1E-9);
        }
    }
}

TEST_F(MathUtilsUnitTest, testLog10OneMinusPow10)
{
    std::pmr::vector<double> keys({std::exp(1), 10, 1, 0, -1, -3, -10, -30, -100, -300, -1000, -3000}, _pool);
    double actual, expected;
    for (double x : keys) {
        actual = math_utils::log10one_minus_pow10(x);
        expected = std::log10(1 - std::pow(10.0, x));
        if (std::isnan(expected)) {
            ASSERT_TRUE(std::isnan(actual));
        }
        else if (std::isinf(actual) || std::isinf(expected)) {
            ASSERT_TRUE(std::isinf(actual) || std::isinf(expected));
        }
        else {
            ASSERT_NEAR(actual, expected, 1E-9);
        }
    }
}

TEST_F(MathUtilsUnitTest, testApplyToArray)
{
    std::pmr::vector<std::pmr::vector<double>> values(
        {{1, 2, 3, 4, 5}, {0.0}, {3, 2, 5, 6}, {19, -5, 22, 55, -1000, 2, 2, 2}, {-1, -1, -1, -1, -1}, {-1, -2, -3, -10, -1}}, _pool);
    for (const std::pmr::vector<double>& arr : values) {
        std::pmr::vector<double> copy(arr, _pool);

        std::pmr::vector<double> result1 = math_utils::apply_to_array(arr, [](double d) -> double { return std::exp(d); });

        std::pmr::vector<double> result2(_pool);
        result2.reserve(copy.size());
        std::for_each(copy.begin(), copy.end(), [&result2](double d) { result2.emplace_back(std::exp(d)); });

        for (size_t i = 0, len = arr.size(); i < len; ++i) {
            ASSERT_NEAR(result1.at(i), result2.at(i), _precision);
        }

        // make sure original array was not affected
        for (size_t i = 0, len = arr.size(); i < len; ++i) {
            ASSERT_NEAR(arr.at(i), copy.at(i), _precision);
        }
    }
}

TEST_F(MathUtilsUnitTest, testApplyToArrayInPlace)
{
    std::pmr::vector<std::pmr::vector<double>> values(
        {{1, 2, 3, 4, 5}, {0.0}, {3, 2, 5, 6}, {19, -5, 22, 55, -1000, 2, 2, 2}, {-1, -1, -1, -1, -1}, {-1, -2, -3, -10, -1}}, _pool);
    for (std::pmr::vector<double>& arr : values) {
        std::pmr::vector<double> copy(arr, _pool);

        std::pmr::vector<double>& result1 = math_utils::apply_to_array_in_place(arr, [](double d) -> double { return std::exp(d); });

        std::for_each(copy.begin(), copy.end(), [](double& d) { d = std::exp(d); });

        for (size_t i = 0, len = arr.size(); i < len; ++i) {
            ASSERT_NEAR(result1.at(i), copy.at(i), _precision);
        }

        // make sure original array was affected
        for (size_t i = 0, len = arr.size(); i < len; ++i) {
            ASSERT_NEAR(arr.at(i), copy.at(i), _precision);
        }
    }
}

TEST_F(MathUtilsUnitTest, testNormalizeFromReal)
{
    std::pmr::vector<double> values{{1.0, 2.0, 3.0}, _pool};
    std::pmr::vector<double> actual = math_utils::normalize_sum_to_one(values);
    std::pmr::vector<double> expected{{1.0 / 6.0, 2.0 / 6.0, 3.0 / 6.0}, _pool};

    for (size_t i = 0, len = actual.size(); i < len; ++i) {
        ASSERT_NEAR(expected.at(i), actual.at(i), _precision);
    }
}

TEST_F(MathUtilsUnitTest, testNormalize)
{
    std::pmr::vector<double> values1{{log10(3.0), log10(2.0), log10(1.0)}, _pool};
    std::pmr::vector<double> values2{{log10(3.0), log10(2.0), log10(1.0)}, _pool};
    std::pmr::vector<double> values3{{3.0, 2.0, 1.0}, _pool};
    std::pmr::vector<double> normalized = math_utils::normalize_from_log10to_linear_space(values1);
    std::pmr::vector<double> normalizedLog10 = math_utils::normalize_log10(values2);
    std::pmr::vector<double> normalizedLogInLog10 = math_utils::scale_log_space_array_for_numerical_stability(values3);

    std::pmr::vector<double> normalizedExpected{{3.0 / 6.0, 2.0 / 6.0, 1.0 / 6.0}, _pool};
    std::pmr::vector<double> normalizedLog10Expected{{log10(3.0 / 6.0), log10(2.0 / 6.0), log10(1.0 / 6.0)}, _pool};
    std::pmr::vector<double> normalizedLogInLogExpected{{0.0, -1.0, -2.0}, _pool};

    for (size_t i = 0, len = values1.size(); i < len; ++i) {
        ASSERT_NEAR(normalizedExpected.at(i), normalized.at(i), _precision);
        ASSERT_NEAR(normalizedLog10Expected.at(i), normalizedLog10.at(i), _precision);
        ASSERT_NEAR(normalizedLogInLogExpected.at(i), normalizedLogInLog10.at(i), _precision);
    }
}

TEST_F(MathUtilsUnitTest, testLog10Factorial)
{
    ASSERT_NEAR(MathUtils::log10factorial(4), 1.3802112, 1e-6);
    ASSERT_NEAR(MathUtils::log10factorial(10), 6.559763, 1e-6);
    ASSERT_NEAR(MathUtils::log10factorial(200), 374.896888, 1e-3);
    ASSERT_NEAR(MathUtils::log10factorial(12342), 45138.2626503, 1e-1);

    int small_start = 1;
    int med_start = 200;
    int large_start = 12342;
    double log10Factorial_small = 0;
    double log10Factorial_middle = MathUtils::log10factorial(med_start);
    double log10Factorial_large = MathUtils::log10factorial(large_start);
    for (int i = 1; i < 1000; i++) {
        log10Factorial_small += log10(i + small_start);
        log10Factorial_middle += log10(i + med_start);
        log10Factorial_large += log10(i + large_start);

        ASSERT_NEAR(MathUtils::log10factorial(small_start + i), log10Factorial_small, 1e-6);
        ASSERT_NEAR(MathUtils::log10factorial(med_start + i), log10Factorial_middle, 1e-3);
        ASSERT_NEAR(MathUtils::log10factorial(large_start + i), log10Factorial_large, 1e-1);
    }
}

TEST_F(MathUtilsUnitTest, testLog10BinomialCoefficient)
{
    std::pmr::vector<double> z_vals{{0.999, 0.9, 0.8, 0.5, 0.2, 0.01, 0.0001}, _pool};
    std::pmr::vector<int32_t> exponent{{5, 15, 25, 50, 100}, _pool};

    double logz, expected_log;
    for (double z : z_vals) {
        logz = std::log10(z);
        for (int32_t exp : exponent) {
            expected_log = exp * std::log10(1 + z);
            std::pmr::vector<double> newtonArray_log(1 + exp, 0, _pool);
            for (int32_t k = 0; k <= exp; ++k) {
                newtonArray_log[k] = math_utils::log10binomial_coefficient(exp, k) + k * logz;
            }
            ASSERT_NEAR(math_utils::log10_sum_log10_1(newtonArray_log), expected_log, _precision);
        }
    }

    // results from Wolfram Alpha
    ASSERT_NEAR(math_utils::log10binomial_coefficient(4, 2), 0.7781513, 1e-6);
    ASSERT_NEAR(math_utils::log10binomial_coefficient(10, 3), 2.079181, 1e-6);
    ASSERT_NEAR(math_utils::log10binomial_coefficient(103928, 119), 400.2156, 1e-4);
}

TEST_F(MathUtilsUnitTest, testApproximateLogSumLog)
{
    double requiredPrecision = 1E-4;
    std::pmr::vector<double> values{{0.0, 0.0, 0.0}, _pool};
    ASSERT_NEAR(MathUtils::approximate_log10sum_log10(values), log10(3), requiredPrecision);
    ASSERT_NEAR(MathUtils::approximate_log10sum_log10(values, 0), 0.0, requiredPrecision);
    ASSERT_NEAR(MathUtils::approximate_log10sum_log10(values, 3), log10(3), requiredPrecision);
    ASSERT_NEAR(MathUtils::approximate_log10sum_log10(values, 2), log10(2), requiredPrecision);
    ASSERT_NEAR(MathUtils::approximate_log10sum_log10(values, 1), 0.0, requiredPrecision);

    std::pmr::vector<double> values2{{NEGATIVE_INFINITY, NEGATIVE_INFINITY, NEGATIVE_INFINITY}, _pool};
    ASSERT_TRUE(std::isinf(MathUtils::approximate_log10sum_log10(values2)));
}