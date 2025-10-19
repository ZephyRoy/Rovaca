#ifndef ROVACA_HC_STRAND_ODDS_RATIO_H_
#define ROVACA_HC_STRAND_ODDS_RATIO_H_
#include "strand_bias_test.h"

namespace rovaca
{

/*!
 * @brief 计算Info列：SOR
 * StrandOddsRatio是指所有reads的正负链比值的对数比（log odds ratio）。该值用于评估SNP变异是否偏向于某一条链。
 * 具体来说，该值是指在变异位点上，正链上的reads数与负链上的reads数的比值，除以参考基因组上正链上的reads数与负链上的reads数的比值，然后取对数。
 * 如果该值越大，说明变异位点上的reads更可能来自于正链，如果该值越小，则说明变异位点上的reads更可能来自于负链。
 * 该值的绝对值越大，说明偏向性越明显。StrandOddsRatio的值可以用于评估SNP变异的可信度。
 */
class StrandOddsRatio : public StrandBiasTest
{
    /* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
public:
    static double calculate_sor(const Int32Vector2D& table);

    ~StrandOddsRatio() override = default;

    /* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
protected:
    void calculate_annotation_from_gt_field(pGenotypesContext gc, pInfoData target, pMemoryPool pool) override;

    void calculate_annotation_from_likelihoods(pVariant vc, pRALikelihoods likelihoods, pInfoData target, pMemoryPool pool) override;
};

}  // namespace rovaca

#endif  // ROVACA_HC_STRAND_ODDS_RATIO_H_
