#ifndef ROVACA_HC_FISHER_STRAND_H_
#define ROVACA_HC_FISHER_STRAND_H_
#include "strand_bias_test.h"

namespace rovaca
{

/*!
 * @brief 计算Info列：FS
 * 描述每个位点的链特异性偏差的属性。具体来说，它是一个基于Fisher精确检验的统计量，用于比较正链和负链上的参考和变异碱基的分布情况
 * FisherStrand的值越低，表示正链和负链上的参考和变异碱基的分布情况越均衡，不存在链特异性偏差。
 * 相反，FisherStrand的值越高，表示正链和负链上的参考和变异碱基的分布情况越不均衡，可能存在链特异性偏差。
 */
class FisherStrand : public StrandBiasTest
{
    /* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
public:
    static double p_value_for_contingency_table(const Int32Vector2D& original_table);

    /*!
     * @brief computes the 2-sided pvalue of the fisher's exact test on a normalized table that ensures that the sum of all four entries is
     * less than 2 * 200.
     * @param normalized_table
     * @return
     */
    static double two_sided_pvalue(const Int32Vector2D& normalized_table);

    static double make_value_object_for_annotation(double p_value);
    static double annotation_for_one_table(double p_value) { return make_value_object_for_annotation(p_value); }

    ~FisherStrand() override = default;

    /* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
protected:
    void calculate_annotation_from_gt_field(pGenotypesContext gc, pInfoData target, pMemoryPool pool) override;
    void calculate_annotation_from_likelihoods(pVariant vc, pRALikelihoods likelihoods, pInfoData target, pMemoryPool pool) override;

    /* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
private:
    /**
     * normalize the table so that the entries are not too large.
     * note that this method does not necessarily make a copy of the table being passed in!
     * @param table  the original table
     * @return a normalized version of the table or the original table if it is already normalized
     */
    static Int32Vector2D normalize_contingency_table(const Int32Vector2D& table);
};

}  // namespace rovaca

#endif  // ROVACA_HC_FISHER_STRAND_H_
