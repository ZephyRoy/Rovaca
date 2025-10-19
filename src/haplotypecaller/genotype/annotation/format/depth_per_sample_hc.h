#ifndef ROVACA_HC_DEPTH_PER_SAMPLE_HC_H_
#define ROVACA_HC_DEPTH_PER_SAMPLE_HC_H_
#include "interface/interface_genotype_annotation.hpp"

namespace rovaca
{

/*!
 * @brief 计算Format列：DP
 * 在gVCF文件中，DP列是一个必需的FORMAT列，表示每个位点的总深度。DP的全称是Depth，表示该位点的测序深度，即该位点被测序的总次数。
 *
 * DP列的值是一个整数，表示该位点的总深度。例如，一个DP值为10的位点表示该位点被测序了10次。
 *
 * DP列的值可以用于检测测序深度是否足够，以及检测样本污染和杂合性等问题。
 * 如果一个位点的DP值过低，那么可能需要增加测序深度，以提高分型质量。如果一个位点的DP值过高，那么可能存在样本污染或杂合性等问题，需要进一步检查。
 *
 * 需要注意的是，DP列是一个必需的FORMAT列，所有的gVCF文件都包含该列。在使用gVCF文件进行变异检测和分型时，DP列的值是非常重要的，可以影响变异检测的灵敏度和特异性。
 */
class DepthPerSampleHC : public InterfaceGenotypeAnnotation
{
    /* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
public:
    void annotate(pRefFragment ref, pVariant vc, pGenotype g, pRALikelihoods likelihoods, pMemoryPool pool) override;

    ~DepthPerSampleHC() override = default;
};

}  // namespace rovaca

#endif  // ROVACA_HC_DEPTH_PER_SAMPLE_HC_H_
