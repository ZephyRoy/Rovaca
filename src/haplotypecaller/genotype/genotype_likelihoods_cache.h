#ifndef ROVACA_HC_GENOTYPE_LIKELIHOODS_CACHE_H_
#define ROVACA_HC_GENOTYPE_LIKELIHOODS_CACHE_H_
#include <memory>

#include "forward.h"

namespace rovaca
{

/*!
 * @brief 此cache仅用于GVCF模式中，不允许使用破坏cache对象的函数，如std::move
 */
class GenotypeLikelihoodsCache
{
    static constexpr std::size_t s_buffer_size = 10240;

    /* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
private:
    uint8_t _data[s_buffer_size];
    std::pmr::monotonic_buffer_resource _pool;

    std::pmr::vector<std::pmr::vector<pGenotypeLikelihoods>> _cache;

    /* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
public:
    static pGenotypeLikelihoodsCache get_instance()
    {
        static std::unique_ptr<GenotypeLikelihoodsCache> cache(new GenotypeLikelihoodsCache{});
        return cache.get();
    }

    pGenotypeLikelihoods get_indel_pls(int32_t ploidy, int32_t n_informative_reads);

    /* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
private:
    GenotypeLikelihoodsCache();

    void init_cache();
};

}  // namespace rovaca

#endif  // ROVACA_HC_GENOTYPE_LIKELIHOODS_CACHE_H_
