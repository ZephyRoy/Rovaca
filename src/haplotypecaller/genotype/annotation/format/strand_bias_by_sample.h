#ifndef ROVACA_HC_STRAND_BIAS_BY_SAMPLE_H_
#define ROVACA_HC_STRAND_BIAS_BY_SAMPLE_H_
#include "interface/interface_genotype_annotation.hpp"

namespace rovaca
{

/*!
 * @brief 计算Format列：SB
 * 在gVCF文件中，SB列是一个可选的FORMAT列，表示单个碱基的分型质量分布。SB的全称是Strand Bias，表示该碱基在正链和负链上的分型质量分布。
 *
 * SB列的值是一个四元组，分别表示在正链和负链上，参考碱基和非参考碱基的分型质量。
 * 例如，一个SB值为0,0,0,0的碱基表示在正链和负链上，参考碱基和非参考碱基的分型质量都为0，即该碱基没有分型质量偏差。
 *
 * SB列的值可以帮助检测分型质量偏差，例如由于PCR扩增偏差、测序偏差或样本污染等原因导致的偏差。
 * 如果一个碱基的SB值明显偏离了期望的分布，那么可能需要进一步检查该碱基的分型质量，以确定是否存在偏差。
 *
 * 需要注意的是，SB列是一个可选的FORMAT列，不是所有的gVCF文件都包含该列。如果你的gVCF文件中没有SB列，那么就无法使用该列来检测分型质量偏差。
 */
class StrandBiasBySample : public InterfaceGenotypeAnnotation
{
    /* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
public:
    void annotate(pRefFragment ref, pVariant vc, pGenotype g, pRALikelihoods likelihoods, pMemoryPool pool) override;

    static Int32Vector get_contingency_array(const Int32Vector2D& table, pMemoryPool pool);

    ~StrandBiasBySample() override = default;
};

}  // namespace rovaca

#endif  // ROVACA_HC_STRAND_BIAS_BY_SAMPLE_H_
