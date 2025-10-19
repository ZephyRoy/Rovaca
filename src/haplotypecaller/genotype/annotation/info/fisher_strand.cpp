#include "fisher_strand.h"

#include "index_range.hpp"
#include "info_data.hpp"
#include "math_utils.h"
#include "quality_utils.h"
#include "utils/hypergeometric_distribution.h"

namespace rovaca
{

static constexpr double MIN_PVALUE = 1E-320;
static constexpr double REL_ERR = 1 - 10e-7;
static constexpr int32_t s_min_count = 2;
static constexpr int32_t s_target_table_size = 200;

void FisherStrand::calculate_annotation_from_gt_field([[maybe_unused]] pGenotypesContext gc, [[maybe_unused]] pInfoData target,
                                                      [[maybe_unused]] pMemoryPool pool)
{}

void FisherStrand::calculate_annotation_from_likelihoods(pVariant vc, pRALikelihoods likelihoods, pInfoData target, pMemoryPool pool)
{
    Int32Vector2D table = get_contingency_table(vc, likelihoods, s_min_count, pool);
    double p_value = p_value_for_contingency_table(table);
    double fs = annotation_for_one_table(p_value);
    target->set_fs(fs);
}

double FisherStrand::p_value_for_contingency_table(const Int32Vector2D& original_table)
{
    Int32Vector2D normalized_table = normalize_contingency_table(original_table);
    return two_sided_pvalue(normalized_table);
}

double FisherStrand::two_sided_pvalue(const Int32Vector2D& normalized_table)
{
    const Int32Vector2D& x = normalized_table;
    int32_t m = x[0][0] + x[0][1];
    int32_t n = x[1][0] + x[1][1];
    int32_t k = x[0][0] + x[1][0];
    int32_t lo = std::max(0, k - n);
    int32_t hi = std::min(k, m);

    IndexRange support(lo, hi + 1);
    // special case, support has only one value
    if (support.size() <= 1) {
        return 1.0;
    }

    pMemoryPool pool = normalized_table.get_allocator().resource();
    HypergeometricDistribution dist(m + n, m, k);
    DoubleVector logds = support.map_to_double([&](int32_t i) { return dist.log_probability(i); }, pool);
    double threshold = logds[x[0][0] - lo] * REL_ERR;

    DoubleVector log10ds(pool);
    std::for_each(logds.begin(), logds.end(), [&](double d) {
        if (d <= threshold) {
            log10ds.push_back(math_utils::log_to_log10(d));
        }
    });
    double p_value = math_utils::sum_log10(log10ds);
    return std::min(1.0, p_value);
}

double FisherStrand::make_value_object_for_annotation(double p_value)
{
    return QualityUtils::phred_scale_error_rate(std::max(p_value, MIN_PVALUE));
}

Int32Vector2D FisherStrand::normalize_contingency_table(const Int32Vector2D& table)
{
    int32_t sum = table[0][0] + table[0][1] + table[1][0] + table[1][1];
    if (sum <= s_target_table_size * 2) {
        return {table, table.get_allocator()};
    }

    double norm_factor = (double)sum / (double)s_target_table_size;

    return {{
                {(int32_t)(table[0][0] / norm_factor), (int32_t)(table[0][1] / norm_factor)},
                {(int32_t)(table[1][0] / norm_factor), (int32_t)(table[1][1] / norm_factor)},
            },
            table.get_allocator()};
}

}  // namespace rovaca