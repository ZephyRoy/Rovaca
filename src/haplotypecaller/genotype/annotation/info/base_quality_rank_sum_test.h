#ifndef ROVACA_HC_BASE_QUALITY_RANK_SUM_TEST_H_
#define ROVACA_HC_BASE_QUALITY_RANK_SUM_TEST_H_
#include "rank_sum_test.h"

namespace rovaca
{

/*!
 * @brief 计算Info列：BaseQRankSum
 * 基于Wilcoxon秩和检验的统计量，用于比较参考和变异碱基的质量分数分布。
 * BaseQRankSum的值越高，表示变异碱基的质量分数相对于参考碱基更高，因此更可能是真实的变异。
 * 相反，BaseQRankSum的值越低，表示变异碱基的质量分数相对于参考碱基更低，因此更可能是假阳性。
 */
class BaseQualityRankSumTest : public RankSumTest
{
    /* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
public:
    void set_values_to_info(pInfoData info, void *data) override;

    OptionalDouble get_element_for_read(pReadRecord read, pVariant vc) override { return get_read_base_quality(read, vc); }

    static OptionalDouble get_read_base_quality(pReadRecord read, pVariant vc);

    ~BaseQualityRankSumTest() override = default;
};

}  // namespace rovaca

#endif  // ROVACA_HC_BASE_QUALITY_RANK_SUM_TEST_H_
