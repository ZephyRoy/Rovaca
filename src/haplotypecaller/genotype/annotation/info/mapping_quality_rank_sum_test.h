#ifndef ROVACA_HC_MAPPING_QUALITY_RANK_SUM_TEST_H_
#define ROVACA_HC_MAPPING_QUALITY_RANK_SUM_TEST_H_
#include "rank_sum_test.h"

namespace rovaca
{

/*!
 * @brief 计算Info列：MQRankSum
 * 描述每个位点的比对质量的属性。具体来说，它是一个基于Wilcoxon秩和检验的统计量，用于比较参考和变异碱基的比对质量分布情况。
 * MappingQualityRankSumTest的值越高，表示参考和变异碱基的比对质量分布越不均衡，可能存在比对错误或者测序偏差等问题。
 * 相反，MappingQualityRankSumTest的值越低，表示参考和变异碱基的比对质量分布越均衡，可能存在更高的比对准确性和测序质量。
 */
class MappingQualityRankSumTest : public RankSumTest
{
    /* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
public:
    void set_values_to_info(pInfoData info, void *data) override;

    OptionalDouble get_element_for_read(pReadRecord read, pVariant vc) override;

    ~MappingQualityRankSumTest() override = default;
};

}  // namespace rovaca

#endif  // ROVACA_HC_MAPPING_QUALITY_RANK_SUM_TEST_H_
