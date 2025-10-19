#ifndef ROVACA_HC_QUAL_BY_DEPTH_H_
#define ROVACA_HC_QUAL_BY_DEPTH_H_
#include "interface/interface_info_field_annotation.hpp"

namespace rovaca
{

/*!
 * @brief 计算Info列：QD
 * 描述每个位点的质量深度比的属性。具体来说，它是一个浮点数，表示该位点上的质量值和深度值之比。
 * QualByDepth的值越高，表示该位点上的质量值和深度值之比越高，可能存在更高的变异质量和可靠性。
 * 相反，QualByDepth的值越低，表示该位点上的质量值和深度值之比越低，可能存在更多的测序偏差和错误。
 */
class QualByDepth : public InterfaceInfoFieldAnnotation
{
    /* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
public:
    void annotate(pRefFragment ref, pVariant vc, pRALikelihoods likelihoods, pInfoData target, pMemoryPool pool) override;

    static int32_t get_depth(pGenotypesContext gc, pRALikelihoods likelihoods);

    ~QualByDepth() override = default;
};

}  // namespace rovaca

#endif  // ROVACA_HC_QUAL_BY_DEPTH_H_
