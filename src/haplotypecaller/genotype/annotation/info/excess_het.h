#ifndef ROVACA_HC_EXCESS_HET_H_
#define ROVACA_HC_EXCESS_HET_H_
#include "interface/interface_info_field_annotation.hpp"

namespace rovaca
{

/*!
 * @brief 计算Info列：ExcessHet
 * 描述每个位点的杂合度的属性。具体来说，它是一个浮点数，表示在该位点上观察到的杂合度与期望杂合度之间的差异。
 * ExcessHet的值越高，表示该位点上观察到的杂合度与期望杂合度之间的差异越大，可能存在杂合性较高的个体或者样本污染等问题。
 * 相反，ExcessHet的值越低，表示该位点上观察到的杂合度与期望杂合度之间的差异越小，可能存在杂合性较低的个体或者样本缺失等问题。
 */
class ExcessHet : public InterfaceInfoFieldAnnotation
{
    /* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
public:
    void annotate(pRefFragment ref, pVariant vc, pRALikelihoods likelihoods, pInfoData target, pMemoryPool pool) override;

    static std::pair<int32_t, double> calculate_eh(pGenotypeCounts t, int32_t sample_count, pMemoryPool pool);

    ~ExcessHet() override = default;

    /* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
private:
    static std::pair<int32_t, double> calculate_eh(pVariant vc, pGenotypesContext genotypes, pMemoryPool pool);

    /**
     * note that this method is not accurate for very small p-values. beyond 1.0e-16 there is no guarantee that the p-value is accurate,
     * just that it is in fact smaller than 1.0e-16 (and therefore we should filter it). it would be more computationally expensive to
     * calculate accuracy beyond a given threshold. here we have enough accuracy to filter anything below a p-value of 10e-6.
     * @param het_count number of observed hets (n_ab)
     * @param ref_count number of observed hom_refs (n_aa)
     * @param hom_count number of observed hom_vars (n_bb)
     * @return right sided p-value or the probability of getting the observed or higher number of hets given the sample size (n) and the
     * observed number of allele a (rare_copies)
     */
    static double exact_test(int32_t het_count, int32_t ref_count, int32_t hom_count, pMemoryPool pool);
};

}  // namespace rovaca

#endif  // ROVACA_HC_EXCESS_HET_H_
