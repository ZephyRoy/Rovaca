#ifndef ROVACA_HC_MATH_UTILS_H_
#define ROVACA_HC_MATH_UTILS_H_
#include <array>
#include <boost/math/special_functions/digamma.hpp>
#include <boost/math/special_functions/erf.hpp>
#include <cmath>
#include <cstdint>
#include <functional>
#include <memory>
#include <memory_resource>
#include <numeric>
#include <vector>

#include "rovaca_logger.h"

#define POSITIVE_INFINITY (INFINITY)
#define NEGATIVE_INFINITY (-INFINITY)

namespace rovaca
{

/*!
 * @brief 与类无关的函数放在这里 😄
 */
namespace math_utils
{

static const double s_log_10 = std::log(10.0);
static const double s_log10_e = std::log10(std::exp(1.0));
static const double s_inv_log_2 = 1.0 / std::log(2.0);
static const double s_inv_log_10 = 1.0 / std::log(10.0);
static const double s_log10_p_of_zero = -1000000.0;
static const double s_log_one_third = -std::log(3.0);
static const double s_log10_one_half = std::log10(0.5);
static const double s_log10_one_third = -std::log10(3.0);
static const double s_log1mexp_threshold = log(0.5);

double log_to_log10(double ln);
double log10_to_log(double log10);
double erf_inv(double x);
double log10gamma(double x);
int32_t fast_round(double d);
double log1mexp(double a);
double log10one_minus_pow10(double a);
double sum_log10(const std::pmr::vector<double>& log10_values);

std::pmr::vector<double> ebe_add(const std::pmr::vector<double>& a, const std::pmr::vector<double>& b);
std::pmr::vector<double> ebe_subtract(const std::pmr::vector<double>& a, const std::pmr::vector<double>& b);

template <typename T>
size_t max_element_index(const std::pmr::vector<T>& arr, size_t begin, size_t end)
{
    CHECK_CONDITION_EXIT(begin >= end, "invalid range specified");
    size_t max_index = begin;
    for (size_t i = begin + 1; i < end; i++) {
        if (arr[i] > arr[max_index]) {
            max_index = i;
        }
    }
    return max_index;
}

template <typename T>
size_t max_element_index(const std::pmr::vector<T>& arr)
{
    return max_element_index(arr, 0, arr.size());
}

std::pmr::vector<double> apply_to_array(const std::pmr::vector<double>& arr, const std::function<double(double)>& func);
std::pmr::vector<double> apply_to_array(const std::pmr::vector<int32_t>& arr, const std::function<double(int32_t)>& func);
std::pmr::vector<double>& apply_to_array_in_place(std::pmr::vector<double>& arr, const std::function<double(double)>& func);

void add_to_array_in_place(std::pmr::vector<double>& array, const std::pmr::vector<double>& summand);
void add_to_array_in_place(std::pmr::vector<int32_t>& array, const std::pmr::vector<int32_t>& summand);

/*! @brief normalizes the real-space probability array */
std::pmr::vector<double> normalize_sum_to_one(const std::pmr::vector<double>& arr);

/*!
 * @brief Given an array of log space (log or log10) values, subtract all values by the array maximum so that the max element in log space
 * is zero.
 * This is equivalent to dividing by the maximum element in real space and is useful for avoiding underflow/overflow when the array's values
 * matter only up to an arbitrary normalizing factor, for example, an array of likelihoods.
 */
std::pmr::vector<double>& scale_log_space_array_for_numerical_stability(std::pmr::vector<double>& arr);

// log10SumLog10
double log10_sum_log10_1(double a, double b);
double log10_sum_log10_1(const std::pmr::vector<double>& log10values, size_t begin, size_t end);
double log10_sum_log10_1(const std::pmr::vector<double>& log10values, size_t begin);
double log10_sum_log10_1(const std::pmr::vector<double>& log10values);

// log10sumLog10
double log10_sum_log10_2(const std::pmr::vector<double>& log10values, size_t begin, size_t end);
double log10_sum_log10_2(const std::pmr::vector<double>& log10values, size_t begin);
double log10_sum_log10_2(const std::pmr::vector<double>& log10values);

std::pmr::vector<double> normalize_log10(std::pmr::vector<double>& arr);
std::pmr::vector<double> normalize_log10(std::pmr::vector<double>& arr, bool take_log10of_output, bool in_place);

std::pmr::vector<double> normalize_from_log10to_linear_space(std::pmr::vector<double>& arr);

std::pmr::vector<double> normalize_from_log10(const std::pmr::vector<double>& arr);
std::pmr::vector<double> normalize_from_log10(const std::pmr::vector<double>& arr, bool take_log10of_output);
std::pmr::vector<double> normalize_from_log10(const std::pmr::vector<double>& arr, bool take_log10of_output, bool keep_in_log_space);

double binomial_coefficient(int32_t n, int32_t k);
double log10binomial_coefficient(int32_t n, int32_t k);

}  // namespace math_utils

class MathUtils
{
    /* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
private:
    class Log10Cache
    {
        /* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
    private:
        static constexpr int32_t s_max_size = 10000000;
        std::array<double, s_max_size> _cache{};

        /* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
    public:
        Log10Cache(const Log10Cache&) = delete;
        void operator=(const Log10Cache&) = delete;
        ~Log10Cache() = default;

        static Log10Cache* singleton()
        {
            static std::unique_ptr<Log10Cache> s(new Log10Cache{});
            return s.get();
        }

        double get(int32_t i) { return i < s_max_size ? _cache.at(i) : compute(i); }

        /* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
    private:
        static double compute(int32_t i) { return std::log10(i); }

        Log10Cache()
        {
            for (int32_t i = 0; i < s_max_size; ++i) {
                _cache[i] = compute(i);
            }
        }
    };

    class Log10FactorialCache
    {
        /* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
    private:
        static constexpr int32_t s_max_size = 10000;
        std::array<double, s_max_size> _cache{};

        /* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
    public:
        Log10FactorialCache(const Log10FactorialCache&) = delete;
        void operator=(const Log10FactorialCache&) = delete;
        ~Log10FactorialCache() = default;

        static Log10FactorialCache* singleton()
        {
            static std::unique_ptr<Log10FactorialCache> s{new Log10FactorialCache{}};
            return s.get();
        }

        double get(int32_t i) { return i < s_max_size ? _cache.at(i) : compute(i); }

        /* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
    private:
        static double compute(int32_t i) { return math_utils::log10gamma(i + 1); }

        Log10FactorialCache()
        {
            for (int32_t i = 0; i < s_max_size; ++i) {
                _cache[i] = compute(i);
            }
        }
    };

#if 0
    class DigammaCache
    {
        /* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
    private:
        static constexpr int32_t s_max_size = 100000;
        std::array<double, s_max_size> _cache{};

        /* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
    public:
        DigammaCache(const DigammaCache&) = delete;
        void operator=(const DigammaCache&) = delete;
        DigammaCache()
        {
            for (int32_t i = 0; i < s_max_size; ++i) {
                _cache[i] = compute(i);
            }
        }
        ~DigammaCache() = default;
        double get(int32_t i) { return i < s_max_size ? _cache.at(i) : compute(i); }

        /* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
    private:
        static double compute(int32_t i) { return boost::math::digamma(i); }
    };
#endif

    class JacobianLogTable
    {
        /* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
    public:
        static constexpr double s_max_tolerance = 8.0L;

        /* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
    private:
        static constexpr double s_table_step = 0.0001L;
        static constexpr double s_inv_step = 1.0 / s_table_step;

        std::array<double, static_cast<size_t>(s_max_tolerance / s_table_step) + 1> _cache{};

        /* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
    public:
        JacobianLogTable(const JacobianLogTable&) = delete;
        void operator=(const JacobianLogTable&) = delete;
        static JacobianLogTable* singleton()
        {
            static std::unique_ptr<JacobianLogTable> s{new JacobianLogTable{}};
            return s.get();
        }

        double get(double difference) { return _cache.at((size_t)std::round(difference * s_inv_step)); }

        /* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
    private:
        JacobianLogTable()
        {
            for (size_t i = 0, size = _cache.size(); i < size; ++i) {
                _cache[i] = std::log10(1.0 + std::pow(10.0, -s_table_step * (double)i));
            }
        }
    };

    /* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
public:
    static double log10factorial(int32_t i) { return Log10FactorialCache::singleton()->get(i); }
    static double log10(int32_t i) { return Log10Cache::singleton()->get(i); }

    static double approximate_log10sum_log10(const std::pmr::vector<double>& values, size_t begin, size_t end);
    static double approximate_log10sum_log10(const std::pmr::vector<double>& values, size_t end);
    static double approximate_log10sum_log10(const std::pmr::vector<double>& values)
    {
        return approximate_log10sum_log10(values, values.size());
    }
    static double approximate_log10sum_log10(double a, double b);
    static double approximate_log10sum_log10(double a, double b, double c)
    {
        return approximate_log10sum_log10(a, approximate_log10sum_log10(b, c));
    }
};

}  // namespace rovaca

#endif  // ROVACA_HC_MATH_UTILS_H_
