#ifndef ROVACA_HC_INTERFACE_GENOTYPE_ANNOTATION_H_
#define ROVACA_HC_INTERFACE_GENOTYPE_ANNOTATION_H_
#include "forward.h"

namespace rovaca
{

class InterfaceGenotypeAnnotation
{
    /* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
public:
    /*!
     * @brief computes the annotation for the given genotype and the likelihoods per read. expected to modified the passed genotype builder
     * @param ref
     * @param vc
     * @param g
     * @param likelihoods
     */
    virtual void annotate(pRefFragment ref, pVariant vc, pGenotype g, pRALikelihoods likelihoods, pMemoryPool pool) = 0;

    virtual ~InterfaceGenotypeAnnotation() = default;
};

}  // namespace rovaca

#endif  // ROVACA_HC_INTERFACE_GENOTYPE_ANNOTATION_H_
