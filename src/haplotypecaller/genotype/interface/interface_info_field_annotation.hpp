#ifndef ROVACA_HC_INTERFACE_INFO_FIELD_ANNOTATION_H_
#define ROVACA_HC_INTERFACE_INFO_FIELD_ANNOTATION_H_
#include "forward.h"

namespace rovaca
{

class InterfaceInfoFieldAnnotation
{
    /* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
public:
    /*!
     * @brief computes the annotation for the given variant and the likelihoods per read
     * @param ref in|
     * @param vc in|
     * @param likelihoods in|
     * @param info out|
     */
    virtual void annotate(pRefFragment ref, pVariant vc, pRALikelihoods likelihoods, pInfoData target, pMemoryPool pool) = 0;

    virtual ~InterfaceInfoFieldAnnotation() = default;
};

}  // namespace rovaca

#endif  // ROVACA_HC_INTERFACE_INFO_FIELD_ANNOTATION_H_
