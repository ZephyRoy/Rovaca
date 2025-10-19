#ifndef ROVACA_HC_GENOTYPE_PRIOR_CALCULATOR_H_
#define ROVACA_HC_GENOTYPE_PRIOR_CALCULATOR_H_
#include <vector>

#include "forward.h"
#include "genotype_macors.h"

namespace rovaca
{

class GenotypePriorCalculator
{
#define NUMBER_OF_ALLELE_TYPES 4

    /* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
private:
    double _het_values[NUMBER_OF_ALLELE_TYPES]{};
    double _hom_values[NUMBER_OF_ALLELE_TYPES]{};
    double _diff_values[NUMBER_OF_ALLELE_TYPES]{};

    /* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
public:
    DISALLOW_COPY_AND_ASSIGN(GenotypePriorCalculator);
    static pGenotypePriorCalculator assuming_hw(double snp_het, double indel_het);

    /* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
private:
    GenotypePriorCalculator(double snp_het, double snp_hom, double indel_het, double indel_hom, double other_het, double other_hom);
};

}  // namespace rovaca

#endif  // ROVACA_HC_GENOTYPE_PRIOR_CALCULATOR_H_
