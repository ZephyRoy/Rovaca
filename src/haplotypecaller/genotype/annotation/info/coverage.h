#ifndef ROVACA_HC_COVERAGE_H_
#define ROVACA_HC_COVERAGE_H_
#include "interface/interface_info_field_annotation.hpp"

namespace rovaca
{

/*!
 * @brief 计算Info列：DP
 * 描述每个位点的测序深度的属性
 * Coverage属性可以帮助我们了解每个位点的测序深度，从而更好地评估每个位点的可靠性和质量。
 * 通常情况下，高测序深度可以提高变异检测的灵敏度和准确性，但也可能会增加假阳性的风险。
 * 因此，在进行变异鉴定时，需要综合考虑多个属性和其他信息，以确定每个位点的可靠性和质量
 */
class Coverage : public InterfaceInfoFieldAnnotation
{
    /* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
public:
    void annotate(pRefFragment ref, pVariant vc, pRALikelihoods likelihoods, pInfoData target, pMemoryPool pool) override;

    ~Coverage() override = default;
};

}  // namespace rovaca

#endif  // ROVACA_HC_COVERAGE_H_
