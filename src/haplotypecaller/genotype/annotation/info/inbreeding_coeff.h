#ifndef ROVACA_HC_INBREEDING_COEFF_H_
#define ROVACA_HC_INBREEDING_COEFF_H_
#include "interface/interface_info_field_annotation.hpp"

namespace rovaca
{

/*!
 * @brief 计算Info列：InbreedingCoeff
 * 描述每个样本的近亲交配系数的属性。具体来说，它是一个浮点数，表示该样本的近亲交配系数。
 * InbreedingCoeff的值越高，表示该样本的近亲交配程度越高，可能存在遗传缺陷或者疾病等问题。
 * 相反，InbreedingCoeff的值越低，表示该样本的近亲交配程度越低，可能存在更多的遗传多样性和健康基因。
 */
class InbreedingCoeff : public InterfaceInfoFieldAnnotation
{
    /* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
public:
    void annotate(pRefFragment ref, pVariant vc, pRALikelihoods likelihoods, pInfoData target, pMemoryPool pool) override;

    ~InbreedingCoeff() override = default;
};

}  // namespace rovaca

#endif  // ROVACA_HC_INBREEDING_COEFF_H_
