#include "hypergeometric_distribution.h"

#include "utils/saddle_point_expansion.h"

namespace rovaca
{

double HypergeometricDistribution::log_probability(int32_t x) const
{
    double ret;
    std::pair<int32_t, int32_t> domain = get_domain(_population_size, _number_of_successes, _sample_size);

    if (x < domain.first || x > domain.second) {
        ret = NEGATIVE_INFINITY;
    }
    else {
        double p = (double)_sample_size / (double)_population_size;
        double q = (double)(_population_size - _sample_size) / (double)_population_size;
        double p1 = SaddlePointExpansion::log_binomial_probability(x, _number_of_successes, p, q);
        double p2 = SaddlePointExpansion::log_binomial_probability(_sample_size - x, _population_size - _number_of_successes, p, q);
        double p3 = SaddlePointExpansion::log_binomial_probability(_sample_size, _population_size, p, q);
        ret = p1 + p2 - p3;
    }

    return ret;
}

}  // namespace rovaca