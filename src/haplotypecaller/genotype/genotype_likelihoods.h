#ifndef ROVACA_HC_GENOTYPE_LIKELIHOODS_H_
#define ROVACA_HC_GENOTYPE_LIKELIHOODS_H_
#include "forward.h"
#include "genotype_num_likelihoods_cache.hpp"

namespace rovaca
{

struct GenotypeLikelihoodsAllelePair
{
    int32_t allele_index1{INVALID_INT};
    int32_t allele_index2{INVALID_INT};
};

class GenotypeLikelihoods
{
    /* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
private:
    static constexpr size_t s_arr_size = 1326;

    static pGenotypeNumLikelihoodsCache s_num_likelihood_cache;
    static std::array<GenotypeLikelihoodsAllelePair, s_arr_size> s_diploid_plindex_to_allele_index;

    /* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
public:
    static constexpr size_t s_max_diploid_alt_alleles_that_can_be_genotyped = 50;

    DoubleVector _log10likelihoods;
    Int32Vector _pls;

    static pGenotypeLikelihoods create(DoubleVector&& log10likelihoods, pMemoryPool pool);
    static pGenotypeLikelihoods create(Int32Vector&& likelihoods, pMemoryPool pool);
    static pGenotypeLikelihoods create(const Int32Vector& likelihoods, pMemoryPool pool);

    /*!
     * @brief as per the vcf spec: "the ordering of genotypes for the likelihoods is given by: f(j/k) = (k*(k+1)/2)+j.
     * in other words, for biallelic sites the ordering is: aa,ab,bb; for triallelic sites the ordering is: aa,ab,bb,ac,bc,cc, etc."
     * assumes that allele1index < allele2index
     * @param allele1index
     * @param allele2index
     * @return
     */
    static int32_t calculate_plindex(int32_t allele1index, int32_t allele2index)
    {
        return (allele2index * (allele2index + 1) / 2) + allele1index;
    }

    /*!
     * @brief get the diploid allele index pair for the given pl index
     * @param pl_index
     * @return
     */
    static const GenotypeLikelihoodsAllelePair& get_allele_pair(int32_t pl_index);

    /**
     * compute how many likelihood elements are associated with the given number of alleles
     * equivalent to asking in how many ways n non-negative integers can add up to p is s(n,p)
     * where p = ploidy (number of chromosomes) and n = total # of alleles.
     * each chromosome can be in one single state (0,...,n-1) and there are p of them.
     * naive solution would be to store n*p likelihoods, but this is not necessary because we can't distinguish chromosome states, but
     * rather only total number of alt allele counts in all chromosomes.
     *
     * for example, s(3,2) = 6: for alleles a,b,c, on a diploid organism we have six possible genotypes:
     * aa,ab,bb,ac,bc,cc.
     * another way of expressing is with vector (#of a alleles, # of b alleles, # of c alleles)
     * which is then, for ordering above, (2,0,0), (1,1,0), (0,2,0), (1,1,0), (0,1,1), (0,0,2)
     * in general, for p=2 (regular biallelic), then s(n,2) = n*(n+1)/2
     *
     * note this method caches the value for most common num allele / ploidy combinations for efficiency
     *
     * for non-cached values, the result is calculated via a call to calc_num_likelihoods,
     * which uses the apache commons combinatorics_utils class
     * using the formula (num_alleles + ploidy - 1) choose ploidy
     *
     *   @param  num_alleles      number of alleles (including ref)
     *   @param  ploidy          ploidy, or number of chromosomes in set
     *   @return    number of likelihood elements we need to hold.
     */
    static int32_t num_likelihoods(int32_t num_alleles, int32_t ploidy) { return s_num_likelihood_cache->get(num_alleles, ploidy); }

    static double get_gq_log10from_likelihoods(int32_t i_of_choosen_genotype, const DoubleVector& likelihoods, pMemoryPool pool);

    static Int32Vector gls_to_pls(const DoubleVector& gls);
    static DoubleVector pls_to_gls(const Int32Vector& pls);

    static Int32Vector get_plindices_of_alleles(int32_t allele1index, int32_t allele2index, pMemoryPool pool);

    /* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
private:
    explicit GenotypeLikelihoods(DoubleVector&& log10likelihoods)
        : _log10likelihoods(std::move(log10likelihoods))
        , _pls(gls_to_pls(_log10likelihoods))
    {}

    explicit GenotypeLikelihoods(Int32Vector&& likelihoods)
        : _log10likelihoods(pls_to_gls(likelihoods))
        , _pls(std::move(likelihoods))
    {}

    static std::array<GenotypeLikelihoodsAllelePair, s_arr_size> calculate_diploid_plcache();
};

}  // namespace rovaca

#endif  // ROVACA_HC_GENOTYPE_LIKELIHOODS_H_
