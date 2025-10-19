#ifndef ROVACA_HC_CHROMOSOME_COUNTS_H_
#define ROVACA_HC_CHROMOSOME_COUNTS_H_
#include "interface/interface_info_field_annotation.hpp"

namespace rovaca
{

/*!
 * @brief 计算Info列：AN AC AF
 * 描述每个等位基因在每个染色体上的计数的属性，含多个整数的列表，每个整数表示该等位基因在对应染色体上的计数
 * ChromosomeCounts属性可以帮助我们了解每个等位基因在不同染色体上的分布情况，从而更好地评估每个位点的可靠性和质量。
 */
class ChromosomeCounts : public InterfaceInfoFieldAnnotation
{
    /* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
public:
    void annotate(pRefFragment ref, pVariant vc, pRALikelihoods likelihoods, pInfoData target, pMemoryPool pool) override;

    ~ChromosomeCounts() override = default;
};

}  // namespace rovaca

#endif  // ROVACA_HC_CHROMOSOME_COUNTS_H_
