#ifndef ROVACA_HC_HYPERGEOMETRIC_DISTRIBUTION_H_
#define ROVACA_HC_HYPERGEOMETRIC_DISTRIBUTION_H_
#include <cmath>
#include <cstdint>
#include <utility>

#include "genotype_macors.h"

namespace rovaca
{

class HypergeometricDistribution
{
    /* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
private:
    int32_t _population_size;
    int32_t _number_of_successes;
    int32_t _sample_size;

    /* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
public:
    DISALLOW_COPY_AND_ASSIGN(HypergeometricDistribution);
    HypergeometricDistribution(int32_t population_size, int32_t number_of_successes, int32_t sample_size)
        : _population_size(population_size)
        , _number_of_successes(number_of_successes)
        , _sample_size(sample_size)
    {}

    double log_probability(int32_t x) const;

    /* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
private:
    static std::pair<int32_t, int32_t> get_domain(int32_t n, int32_t m, int32_t k)
    {
        return {get_lower_domain(n, m, k), get_upper_domain(m, k)};
    }

    static int32_t get_upper_domain(int32_t m, int32_t k) { return std::min(k, m); }
    static int32_t get_lower_domain(int32_t n, int32_t m, int32_t k) { return std::max(0, m - n + k); }
};

}  // namespace rovaca

#endif  // ROVACA_HC_HYPERGEOMETRIC_DISTRIBUTION_H_
