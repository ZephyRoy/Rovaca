#include "math_utils.h"

namespace rovaca
{

namespace math_utils
{

double log_to_log10(double ln) { return ln * s_log10_e; }

double log10_to_log(double log10) { return log10 * s_log_10; }

double log10gamma(double x) { return log_to_log10(std::lgamma(x)); }

double erf_inv(double x) { return boost::math::erf_inv(x); }

int32_t fast_round(double d) { return (d > 0.0) ? (int32_t)(d + 0.5l) : (int32_t)(d - 0.5l); }

double log1mexp(double a)
{
    if (a > 0) return NAN;
    if (a == 0) return NEGATIVE_INFINITY;
    return (a < s_log1mexp_threshold) ? std::log1p(-std::exp(a)) : std::log(-std::expm1(a));
}

double log10one_minus_pow10(double a)
{
    if (a > 0) return NAN;
    if (a == 0) return NEGATIVE_INFINITY;
    return log1mexp(a * s_log_10) * s_inv_log_10;
}

double sum_log10(const std::pmr::vector<double>& log10_values) { return std::pow(10.0, log10_sum_log10_1(log10_values)); }

std::pmr::vector<double> ebe_add(const std::pmr::vector<double>& a, const std::pmr::vector<double>& b)
{
    CHECK_CONDITION_EXIT(a.size() != b.size(), "a.size() != b.size()");
    std::pmr::vector<double> result(a, a.get_allocator());
    for (size_t i = 0, len = a.size(); i < len; ++i) {
        result[i] += b.at(i);
    }
    return result;
}

std::pmr::vector<double> ebe_subtract(const std::pmr::vector<double>& a, const std::pmr::vector<double>& b)
{
    CHECK_CONDITION_EXIT(a.size() != b.size(), "a.size() != b.size()");
    std::pmr::vector<double> result(a, a.get_allocator());
    for (size_t i = 0, len = a.size(); i < len; ++i) {
        result[i] -= b.at(i);
    }
    return result;
}

std::pmr::vector<double> apply_to_array(const std::pmr::vector<double>& arr, const std::function<double(double)>& func)
{
    std::pmr::vector<double> result(arr.get_allocator());
    result.reserve(arr.size());
    for (const double& d : arr) {
        result.emplace_back(func(d));
    }
    return result;
}

std::pmr::vector<double> apply_to_array(const std::pmr::vector<int32_t>& arr, const std::function<double(int32_t)>& func)
{
    std::pmr::vector<double> result(arr.get_allocator());
    result.reserve(arr.size());
    for (const int32_t& d : arr) {
        result.emplace_back(func(d));
    }
    return result;
}

std::pmr::vector<double>& apply_to_array_in_place(std::pmr::vector<double>& arr, const std::function<double(double)>& func)
{
    for (double& i : arr) {
        i = func(i);
    }
    return arr;
}

std::pmr::vector<double> normalize_sum_to_one(const std::pmr::vector<double>& arr)
{
    std::pmr::vector<double> result(arr.get_allocator());
    if (arr.empty()) {
        result = arr;
    }
    else {
        double sum = std::accumulate(arr.begin(), arr.end(), 0.0);
        CHECK_CONDITION_EXIT(sum < 0, "values in probability array sum to a negative number");
        result = apply_to_array(arr, [&sum](double x) -> double { return x / sum; });
    }
    return result;
}

std::pmr::vector<double>& scale_log_space_array_for_numerical_stability(std::pmr::vector<double>& arr)
{
    double max_value = std::max_element(arr.begin(), arr.end()).operator*();
    return apply_to_array_in_place(arr, [&max_value](double x) -> double { return x - max_value; });
}

double log10_sum_log10_1(double a, double b)
{
    return a > b ? a + std::log10(1 + std::pow(10.0, b - a)) : b + std::log10(1 + std::pow(10.0, a - b));
}

double log10_sum_log10_1(const std::pmr::vector<double>& log10values, size_t begin, size_t end)
{
    if (begin >= end) {
        return NEGATIVE_INFINITY;
    }

    size_t max_value_index = max_element_index(log10values, begin, end);
    double max_value = log10values.at(max_value_index);
    if (max_value == NEGATIVE_INFINITY) {
        return max_value;
    }

    double cur_val, sum = 1.0l;
    for (size_t i = begin; i < end; ++i) {
        cur_val = log10values.at(i);
        if (i == max_value_index || cur_val == NEGATIVE_INFINITY) {
            continue;
        }
        else {
            sum += std::pow(10, cur_val - max_value);
        }
    }
    CHECK_CONDITION_EXIT(std::isnan(sum) || sum == POSITIVE_INFINITY, "log10 p: values must be non-infinite and non-nan");
    return max_value + (sum != 1.0 ? std::log10(sum) : 0.0);
}

double log10_sum_log10_1(const std::pmr::vector<double>& log10values, size_t begin)
{
    return log10_sum_log10_1(log10values, begin, log10values.size());
}

double log10_sum_log10_1(const std::pmr::vector<double>& log10values) { return log10_sum_log10_1(log10values, 0); }

double log10_sum_log10_2(const std::pmr::vector<double>& log10values, size_t begin, size_t end)
{
    if (end < begin + 2) {
        return begin == end ? NEGATIVE_INFINITY : log10values.at(begin);
    }
    size_t max_value_index = max_element_index(log10values, begin, end);
    double max_value = log10values.at(max_value_index);
    if (max_value == NEGATIVE_INFINITY) {
        return max_value;
    }
    double sum = 1.0l;
    for (size_t i = begin; i < end; ++i) {
        sum += i == max_value_index ? 0.0 : std::pow(10.0, log10values.at(i) - max_value);
    }
    CHECK_CONDITION_EXIT(std::isnan(sum) || sum == POSITIVE_INFINITY, "log10 p: values must be non-infinite and non-nan");
    return max_value + std::log10(sum);
}

double log10_sum_log10_2(const std::pmr::vector<double>& log10values, size_t begin)
{
    return log10_sum_log10_2(log10values, begin, log10values.size());
}

double log10_sum_log10_2(const std::pmr::vector<double>& log10values) { return log10_sum_log10_2(log10values, 0); }

std::pmr::vector<double> normalize_log10(std::pmr::vector<double>& arr, bool take_log10of_output, bool in_place)
{
    double log10sum = log10_sum_log10_1(arr);
    std::pmr::vector<double> result(arr.get_allocator());
    if (in_place) {
        result = apply_to_array_in_place(arr, [&](double x) { return x - log10sum; });
    }
    else {
        result = apply_to_array(arr, [&](double x) { return x - log10sum; });
    }
    if (!take_log10of_output) {
        result = apply_to_array_in_place(arr, [](double x) { return std::pow(10.0, x); });
    }
    return result;
}

std::pmr::vector<double> normalize_log10(std::pmr::vector<double>& arr) { return normalize_log10(arr, true, true); }

std::pmr::vector<double> normalize_from_log10to_linear_space(std::pmr::vector<double>& arr) { return normalize_log10(arr, false, true); }

double binomial_coefficient(int32_t n, int32_t k) { return std::pow(10, log10binomial_coefficient(n, k)); }

double log10binomial_coefficient(int32_t n, int32_t k)
{
    return MathUtils::log10factorial(n) - MathUtils::log10factorial(k) - MathUtils::log10factorial(n - k);
}

std::pmr::vector<double> normalize_from_log10(const std::pmr::vector<double>& arr) { return normalize_from_log10(arr, false); }

std::pmr::vector<double> normalize_from_log10(const std::pmr::vector<double>& arr, bool take_log10of_output)
{
    return normalize_from_log10(arr, take_log10of_output, false);
}

std::pmr::vector<double> normalize_from_log10(const std::pmr::vector<double>& arr, bool take_log10of_output, bool keep_in_log_space)
{
    double max_value = std::max_element(arr.begin(), arr.end()).operator*();

    std::pmr::vector<double> normalized(arr.get_allocator());
    normalized.reserve(arr.size());
    if (keep_in_log_space) {
        std::for_each(arr.begin(), arr.end(), [&](double x) { normalized.push_back(x - max_value); });
        return normalized;
    }

    std::for_each(arr.begin(), arr.end(), [&](double x) { normalized.push_back(std::pow(10.0, x - max_value)); });

    double x, sum = std::accumulate(normalized.begin(), normalized.end(), 0.0);
    for (int32_t i = 0, len = (int32_t)arr.size(); i < len; ++i) {
        x = normalized.at(i) / sum;
        if (take_log10of_output) {
            x = std::log10(x);
            if (x < s_log10_p_of_zero || std::isinf(x)) {
                x = arr.at(i) - max_value;
            }
        }
        normalized.at(i) = x;
    }
    return normalized;
}

void add_to_array_in_place(std::pmr::vector<double>& array, const std::pmr::vector<double>& summand)
{
    CHECK_CONDITION_EXIT(array.size() != summand.size(), "arrays must have same length");
    for (size_t i = 0, len = array.size(); i < len; ++i) {
        array[i] += summand.at(i);
    }
}

void add_to_array_in_place(std::pmr::vector<int32_t>& array, const std::pmr::vector<int32_t>& summand)
{
    CHECK_CONDITION_EXIT(array.size() != summand.size(), "arrays must have same length");
    for (size_t i = 0, len = array.size(); i < len; ++i) {
        array[i] += summand.at(i);
    }
}

}  // namespace math_utils

double MathUtils::approximate_log10sum_log10(const std::pmr::vector<double>& values, size_t begin, size_t end)
{
    CHECK_CONDITION_EXIT(values.empty(), "values is empty");
    if (begin == end) {
        return NEGATIVE_INFINITY;
    }

    size_t max_ele_index = std::distance(std::begin(values), std::max_element(values.begin(), values.end()));
    double val, diff, approx_sum = values.at(max_ele_index);
    for (size_t i = begin; i < end; ++i) {
        val = values.at(i);
        if (i == max_ele_index || val == NEGATIVE_INFINITY) {
            continue;
        }
        diff = approx_sum - val;
        approx_sum += diff < JacobianLogTable::s_max_tolerance ? JacobianLogTable::singleton()->get(diff) : 0.0;
    }
    return approx_sum;
}

double MathUtils::approximate_log10sum_log10(const std::pmr::vector<double>& values, size_t end)
{
    size_t max_ele_index = std::distance(std::begin(values), std::max_element(values.begin(), values.end()));
    double val, diff, approx_sum = values.at(max_ele_index);
    for (size_t i = 0; i < end; ++i) {
        val = values.at(i);
        if (i == max_ele_index || val == NEGATIVE_INFINITY) {
            continue;
        }
        diff = approx_sum - val;
        approx_sum += diff < JacobianLogTable::s_max_tolerance ? JacobianLogTable::singleton()->get(diff) : 0.0;
    }
    return approx_sum;
}

double MathUtils::approximate_log10sum_log10(double a, double b)
{
    if (a > b) {
        return approximate_log10sum_log10(b, a);
    }
    else if (a == NEGATIVE_INFINITY) {
        return b;
    }
    double diff = b - a;
    return b + (diff < JacobianLogTable::s_max_tolerance ? JacobianLogTable::singleton()->get(diff) : 0.0);
}

}  // namespace rovaca