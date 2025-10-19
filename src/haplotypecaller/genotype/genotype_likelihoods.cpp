#include "genotype_likelihoods.h"

#include <algorithm>
#include <array>

#include "genotype_macors.h"
#include "math_utils.h"

namespace rovaca
{

static constexpr double s_max_pl = (double)INT32_MAX;
static constexpr size_t s_arr_size = 1326;

std::array<GenotypeLikelihoodsAllelePair, s_arr_size> GenotypeLikelihoods::s_diploid_plindex_to_allele_index = calculate_diploid_plcache();

pGenotypeNumLikelihoodsCache GenotypeLikelihoods::s_num_likelihood_cache = GenotypeNumLikelihoodsCache::get_instance();

pGenotypeLikelihoods GenotypeLikelihoods::create(DoubleVector&& log10likelihoods, pMemoryPool pool)
{
    return new ALLOC_TYPE_IN_POOL(pool, GenotypeLikelihoods) GenotypeLikelihoods{std::forward<DoubleVector&&>(log10likelihoods)};
}

pGenotypeLikelihoods GenotypeLikelihoods::create(Int32Vector&& likelihoods, pMemoryPool pool)
{
    return new ALLOC_TYPE_IN_POOL(pool, GenotypeLikelihoods) GenotypeLikelihoods{std::forward<Int32Vector>(likelihoods)};
}

pGenotypeLikelihoods GenotypeLikelihoods::create(const Int32Vector& likelihoods, pMemoryPool pool)
{
    Int32Vector pls(likelihoods, pool);
    return new ALLOC_TYPE_IN_POOL(pool, GenotypeLikelihoods) GenotypeLikelihoods{std::move(pls)};
}

const GenotypeLikelihoodsAllelePair& GenotypeLikelihoods::get_allele_pair(int32_t pl_index)
{
    CHECK_CONDITION_EXIT(pl_index >= (int32_t)s_arr_size, "out of range");
    return s_diploid_plindex_to_allele_index.at(pl_index);
}

double GenotypeLikelihoods::get_gq_log10from_likelihoods(int32_t i_of_choosen_genotype, const DoubleVector& likelihoods, pMemoryPool pool)
{
    if (likelihoods.empty()) {
        return NEGATIVE_INFINITY;
    }

    double qual = NEGATIVE_INFINITY;
    for (int32_t i = 0, len = (int32_t)likelihoods.size(); i < len; ++i) {
        if (i == i_of_choosen_genotype) {
            continue;
        }
        if (likelihoods.at(i) > qual) {
            qual = likelihoods.at(i);
        }
    }

    qual = likelihoods.at(i_of_choosen_genotype) - qual;
    if (qual < 0.0) {
        // QUAL can be negative if the chosen genotype is not the most likely one individually.
        // In this case, we compute the actual genotype probability and QUAL is the likelihood of it not being the chosen one
        DoubleVector copy_value{likelihoods, pool};
        DoubleVector normalized = math_utils::normalize_from_log10(copy_value);
        return std::log10(1.0 - normalized.at(i_of_choosen_genotype));
    }
    else {
        // invert the size, as this is the probability of making an error
        return -1 * qual;
    }
}

Int32Vector GenotypeLikelihoods::gls_to_pls(const DoubleVector& gls)
{
    Int32Vector pls(gls.get_allocator());
    pls.reserve(gls.size());
    double adjust = std::max_element(gls.begin(), gls.end()).operator*();
    std::for_each(gls.begin(), gls.end(),
                  [&](double d) { pls.emplace_back((int32_t)std::round(std::min(-10.0 * (d - adjust), s_max_pl))); });
    return pls;
}

DoubleVector GenotypeLikelihoods::pls_to_gls(const Int32Vector& pls)
{
    DoubleVector gls(pls.get_allocator());
    gls.reserve(pls.size());
    std::for_each(pls.begin(), pls.end(), [&](int32_t i) { gls.emplace_back(i / -10.0); });
    return gls;
}

Int32Vector GenotypeLikelihoods::get_plindices_of_alleles(int32_t allele1index, int32_t allele2index, pMemoryPool pool)
{
    Int32Vector indexes{3, pool};
    indexes[0] = calculate_plindex(allele1index, allele1index);
    indexes[1] = calculate_plindex(allele1index, allele2index);
    indexes[2] = calculate_plindex(allele2index, allele2index);
    return indexes;
}

std::array<GenotypeLikelihoodsAllelePair, s_arr_size> GenotypeLikelihoods::calculate_diploid_plcache()
{
    std::array<GenotypeLikelihoodsAllelePair, s_arr_size> result{};

    int32_t alt_alleles = s_max_diploid_alt_alleles_that_can_be_genotyped;

    for (int allele1 = 0; allele1 <= alt_alleles; allele1++) {
        for (int allele2 = allele1; allele2 <= alt_alleles; allele2++) {
            result[calculate_plindex(allele1, allele2)] = {allele1, allele2};
        }
    }

    // a bit of sanity checking
    for (int32_t i = 0; i < (int32_t)s_arr_size; i++) {
        CHECK_CONDITION_EXIT(INVALID_INT == result.at(i).allele_index1 || INVALID_INT == result.at(i).allele_index2,
                             "bug: cache entry {} is unexpected null", i);
    }

    return result;
}

}  // namespace rovaca