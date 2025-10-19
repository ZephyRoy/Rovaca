#ifndef ROVACA_HC_GENOTYPE_NUM_LIKELIHOODS_CACHE_H_
#define ROVACA_HC_GENOTYPE_NUM_LIKELIHOODS_CACHE_H_
#include <boost/math/special_functions/binomial.hpp>
#include <cstdint>
#include <memory>

#include "rovaca_logger.h"
#include "genotype_macors.h"

namespace rovaca
{

typedef class GenotypeNumLikelihoodsCache GenotypeNumLikelihoodsCache, *pGenotypeNumLikelihoodsCache;

class GenotypeNumLikelihoodsCache
{
    static constexpr int32_t s_default_ploidy = 3;
    static constexpr int32_t s_default_n_alleles = 8;

    /* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
private:
    int32_t _static_cache[s_default_n_alleles][s_default_ploidy]{};

    /* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
public:
    static pGenotypeNumLikelihoodsCache get_instance()
    {
        static std::unique_ptr<GenotypeNumLikelihoodsCache> cache(new GenotypeNumLikelihoodsCache{});
        return cache.get();
    }

    int32_t get(int32_t num_alleles, int32_t ploidy)
    {
        CHECK_CONDITION_EXIT(num_alleles <= 0 || ploidy <= 0,
                             "num_alleles and ploidy must both exceed 0, but they are num_alleles: {}, ploidy: {}", num_alleles, ploidy);
        if (ROVACA_LIKELY(num_alleles < s_default_n_alleles && ploidy < s_default_ploidy)) {
            return _static_cache[num_alleles - 1][ploidy - 1];
        }

        return (int32_t)boost::math::binomial_coefficient<double>(num_alleles + ploidy - 1, ploidy);
    }

    /* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
private:
    GenotypeNumLikelihoodsCache()
    {
        for (int32_t a = 0; a < s_default_n_alleles; ++a) {
            for (int32_t p = 0; p < s_default_ploidy; ++p) {
                _static_cache[a][p] = (int32_t)boost::math::binomial_coefficient<double>(a + p + 1, p + 1);
            }
        }
    }
};

}  // namespace rovaca

#endif  // ROVACA_HC_GENOTYPE_NUM_LIKELIHOODS_CACHE_H_
