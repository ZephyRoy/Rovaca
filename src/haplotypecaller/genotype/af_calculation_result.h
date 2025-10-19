#ifndef ROVACA_HC_AF_CALCULATION_RESULT_H_
#define ROVACA_HC_AF_CALCULATION_RESULT_H_
#include "allele.h"
#include "forward.h"
#include "genotype_macors.h"

namespace rovaca
{

class AFCalculationResult
{
    /* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
private:
    double _log10posterior_of_no_variant;
    AlleleToDoubleMap _log10p_ref_by_allele;
    Int32Vector _allele_counts_of_mle;
    AlleleVector _alleles_used_in_genotyping;

    /* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
public:
    ~AFCalculationResult() = default;
    DISALLOW_COPY_AND_ASSIGN(AFCalculationResult);
    static pAFCalculationResult create(double log10posterior_of_no_variant, AlleleToDoubleMap&& log10p_ref_by_allele,
                                       Int32Vector&& allele_counts_of_mle, const AlleleVector& alleles_used_in_genotyping,
                                       pMemoryPool pool);

    const AlleleVector& get_alleles_used_in_genotyping() const { return _alleles_used_in_genotyping; }
    int32_t get_allele_count_at_mle(pAllele allele) const;
    const Int32Vector& get_allele_counts_of_mle() const { return _allele_counts_of_mle; }
    double log10prob_only_ref_allele_exists() const { return _log10posterior_of_no_variant; }
    double log10prob_variant_present() const;

    double get_log10posterior_of_allele_absent(pAllele allele);
    bool passes_threshold(pAllele allele, double phred_scale_qual_threshold);

    /* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
private:
    AFCalculationResult(double log10posterior_of_no_variant, AlleleToDoubleMap&& log10p_ref_by_allele, Int32Vector&& allele_counts_of_mle,
                        AlleleVector&& alleles_used_in_genotyping);
};

}  // namespace rovaca

#endif  // ROVACA_HC_AF_CALCULATION_RESULT_H_
