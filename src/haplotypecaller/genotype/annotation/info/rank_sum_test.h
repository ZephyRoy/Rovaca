#ifndef ROVACA_HC_RANK_SUM_TEST_H_
#define ROVACA_HC_RANK_SUM_TEST_H_
#include "forward.h"
#include "genotype_macors.h"
#include "interface/interface_info_field_annotation.hpp"

namespace rovaca
{

class RankSumTest : public InterfaceInfoFieldAnnotation
{
    /* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
public:
    /*!
     * @brief 此类有多个派生配，计算函数主体annotate在此类实现，为了方便每个派生类设置自己的值，使用此函数
     * @param info
     * @param data
     */
    virtual void set_values_to_info(pInfoData info, void* data) = 0;

    virtual OptionalDouble get_element_for_read(pReadRecord read, pVariant vc) = 0;

    void annotate(pRefFragment ref, pVariant vc, pRALikelihoods likelihoods, pInfoData target, pMemoryPool pool) override;

    void fill_quals_from_likelihood(pVariant vc, pRALikelihoods likelihoods, DoubleVector& ref_quals, DoubleVector& alt_quals);

    static bool is_usable_read(pReadRecord read, pVariant vc);

    OptionalDouble get_element_for_read(pReadRecord read, pVariant vc, void* best_allele);

    ~RankSumTest() override = default;
};

}  // namespace rovaca

#endif  // ROVACA_HC_RANK_SUM_TEST_H_
