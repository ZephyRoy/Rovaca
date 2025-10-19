#ifndef ROVACA_HC_READ_POS_RANK_SUM_TEST_H_
#define ROVACA_HC_READ_POS_RANK_SUM_TEST_H_
#include "rank_sum_test.h"

namespace rovaca
{

/*!
 * @brief 计算Info列：ReadPosRankSum
 * ReadPosRankSumTest是指所有reads的位置排名和测试（Rank Sum Test）的结果。
 * 该测试是一种用于检测SNP变异是否偏向于read的某一端的统计方法。
 * 具体来说，该测试会将所有reads按照其在变异位点上的位置排名，然后将排名在变异位点上的reads与排名在变异位点之外的reads的比较得到一个统计值。
 * 如果该值越大，说明变异位点上的reads更可能来自于某一端，如果该值越小，则说明变异位点上的reads更可能来自于另一端。该测试的结果可以用于评估SNP变异的可信度。
 */
class ReadPosRankSumTest : public RankSumTest
{
    /* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
public:
    void set_values_to_info(pInfoData info, void *data) override;

    OptionalDouble get_element_for_read(pReadRecord read, pVariant vc) override { return get_read_position(read, vc); }

    static OptionalDouble get_read_position(pReadRecord read, pVariant vc);

    ~ReadPosRankSumTest() override = default;
};

}  // namespace rovaca

#endif  // ROVACA_HC_READ_POS_RANK_SUM_TEST_H_
