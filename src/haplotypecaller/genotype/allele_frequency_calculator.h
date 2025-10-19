#ifndef ROVACA_HC_ALLELE_FREQUENCY_CALCULATOR_H_
#define ROVACA_HC_ALLELE_FREQUENCY_CALCULATOR_H_
#include "forward.h"
#include "genotype_macors.h"

namespace rovaca
{

class AlleleFrequencyCalculator
{
    /* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
private:
    double _ref_pseudocount;
    double _snp_pseudocount;
    double _indel_pseudocount;
    int32_t _default_ploidy;

    /* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
public:
    ~AlleleFrequencyCalculator() = default;
    DISALLOW_COPY_AND_ASSIGN(AlleleFrequencyCalculator);

    static pAlleleFrequencyCalculator create(double ref_pseudocount, double snp_pseudocount, double indel_pseudocount,
                                             int32_t default_ploidy);

    static pAlleleFrequencyCalculator make_calculator(int32_t ploidy, double snp_heterozygosity, double indel_heterozygosity,
                                                      double heterozygosity_standard_deviation);

    /*! @brief compute the probability of the alleles segregating given the genotype likelihoods of the samples in vc */
    pAFCalculationResult calculate(pVariant vc) const { return calculate(vc, _default_ploidy); }
    pAFCalculationResult calculate(pVariant vc, int32_t default_ploidy) const;

    /* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
private:
    AlleleFrequencyCalculator(double ref_pseudocount, double snp_pseudocount, double indel_pseudocount, int32_t default_ploidy);

    /*! @brief private function that actually calculates allele frequencies etc */
    pAFCalculationResult calculate(int32_t num_allele, const AlleleVector& alleles, pGenotypesContext genotypes, int32_t default_ploidy,
                                   uint32_t ref_length) const;

    /*!
     * @brief effective_allele_counts[allele a] = sum_{genotypes g} (posterior_probability(g) * num_copies of a in g), which we denote as
     * sum [n_g p_g] for numerical stability we will do this in log space: count = sum 10^(log (n_g p_g)) = sum 10^(log n_g + log p_g)
     * thanks to the log-sum-exp trick this lets us work with log posteriors alone
     * @param genotypes
     * @param log10allele_frequencies
     * @return
     */
    static DoubleVector effective_allele_counts(pGenotypesContext genotypes, const DoubleVector& log10allele_frequencies);

    static DoubleVector log10normalized_genotype_posteriors(pGenotype g, pGenotypeLikelihoodCalculator calculator,
                                                            const DoubleVector& log10allele_frequencies);
};

}  // namespace rovaca

#endif  // ROVACA_HC_ALLELE_FREQUENCY_CALCULATOR_H_
