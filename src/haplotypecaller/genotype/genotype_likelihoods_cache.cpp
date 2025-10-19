#include "genotype_likelihoods_cache.h"

#include <algorithm>

#include "genotype_likelihoods.h"
#include "rovaca_logger.h"
#include "math_utils.h"
#include "quality_utils.h"

namespace rovaca
{

#define DEFAULT_PLOIDY                (2)
#define MAX_N_INDEL_INFORMATIVE_READS (40)

pGenotypeLikelihoods GenotypeLikelihoodsCache::get_indel_pls(int32_t ploidy, int32_t n_informative_reads)
{
    int32_t num_reads = n_informative_reads > MAX_N_INDEL_INFORMATIVE_READS ? MAX_N_INDEL_INFORMATIVE_READS : n_informative_reads;
    return _cache.at(ploidy).at(num_reads);
}

GenotypeLikelihoodsCache::GenotypeLikelihoodsCache()
    : _data()
    , _pool(_data, s_buffer_size, std::pmr::null_memory_resource())
    , _cache(&_pool)
{
    init_cache();
}

void GenotypeLikelihoodsCache::init_cache()
{
    QualityUtils::init();

    static const double s_indel_error_rate = -4.5;
    static const uint8_t s_indel_qual = (uint8_t)std::round(s_indel_error_rate * -10.0);
    static const double s_no_indel_likelihood = QualityUtils::qual_to_prob_log10(s_indel_qual);
    static const double s_indel_likelihood = QualityUtils::qual_to_error_prob_log10(s_indel_qual);

    _cache.resize(DEFAULT_PLOIDY + 1);
    std::for_each(_cache.begin(), _cache.end(),
                  [&](std::pmr::vector<pGenotypeLikelihoods>& c) -> void { c.resize(MAX_N_INDEL_INFORMATIVE_READS + 1); });
    std::pmr::vector<pGenotypeLikelihoods>& result = _cache.at(DEFAULT_PLOIDY);

    double denominator = -MathUtils::log10(DEFAULT_PLOIDY);
    DoubleVector first_likelihoods{DEFAULT_PLOIDY + 1, &_pool};
    result[0] = GenotypeLikelihoods::create(std::move(first_likelihoods), &_pool);
    for (int32_t n_informative_reads = 1; n_informative_reads <= MAX_N_INDEL_INFORMATIVE_READS; ++n_informative_reads) {
        DoubleVector pls{DEFAULT_PLOIDY + 1, &_pool};
        pls[0] = n_informative_reads * s_no_indel_likelihood;
        for (int32_t alt_count = 1; alt_count <= DEFAULT_PLOIDY; ++alt_count) {
            double ref_likelihood_accum = s_no_indel_likelihood + MathUtils::log10(DEFAULT_PLOIDY - alt_count);
            double alt_likelihood_accum = s_indel_likelihood + MathUtils::log10(alt_count);
            double sum = MathUtils::approximate_log10sum_log10(ref_likelihood_accum, alt_likelihood_accum);
            pls[alt_count] = n_informative_reads * (sum + denominator);
        }
        result[n_informative_reads] = GenotypeLikelihoods::create(std::move(pls), &_pool);
    }
}

}  // namespace rovaca