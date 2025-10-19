#include "genotype_prior_calculator.h"

#include "rovaca_logger.h"

namespace rovaca
{

#define LOG10_SNP_NORMALIZATION_CONSTANT (0.4771212547)

enum AlleleType { REF = 0, SNP, INDEL, OTHER };

pGenotypePriorCalculator GenotypePriorCalculator::assuming_hw(double snp_het, double indel_het)
{
    CHECK_CONDITION_EXIT(snp_het < 0, "snp-het in log10 scale must be 0 or a negative");
    CHECK_CONDITION_EXIT(indel_het < 0, "indel-het in log10 scale must be 0 or a negative");
    double other_het = std::max(indel_het, snp_het);
    return new GenotypePriorCalculator(snp_het, snp_het * 2, indel_het, indel_het * 2, other_het, other_het * 2);
}

GenotypePriorCalculator::GenotypePriorCalculator(double snp_het, double snp_hom, double indel_het, double indel_hom, double other_het,
                                                 double other_hom)
{
    // REFs: by convention ref log10 likelihoods are set to 0.
    _het_values[REF] = 0.00;
    _hom_values[REF] = 0.00;

    // SNPs: normalized for all possible mutations (number of nucs (4) - 1)
    _het_values[SNP] = snp_het - LOG10_SNP_NORMALIZATION_CONSTANT;
    _hom_values[SNP] = snp_hom - LOG10_SNP_NORMALIZATION_CONSTANT;
    // INDELs:
    _het_values[INDEL] = indel_het;
    _hom_values[INDEL] = indel_hom;
    // Others:
    _het_values[OTHER] = other_het;
    _hom_values[OTHER] = other_hom;

    for (int32_t i = 0; i < NUMBER_OF_ALLELE_TYPES; ++i) {
        _diff_values[i] = _hom_values[i] - _het_values[i];
    }
}

}  // namespace rovaca