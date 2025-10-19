#ifndef ROVACA_HC_RMS_MAPPING_QUALITY_H_
#define ROVACA_HC_RMS_MAPPING_QUALITY_H_
#include "interface/interface_info_field_annotation.hpp"

namespace rovaca
{

/*!
 * @brief 计算Info列：MQ(vcf) 或 RAW_MQandDP(gvcf)
 * RMSMappingQuality是指所有reads的比对质量（mapping quality）的均方根值（root mean square）。
 * 比对质量是指一个read与参考基因组的比对得分，通常用Phred质量值表示，范围是0-60，数值越高表示比对质量越好。
 * RMSMappingQuality的值越高，说明该样本的比对质量越高，可信度也就越高。
 */
class RMSMappingQuality : public InterfaceInfoFieldAnnotation
{
    /* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
private:
    bool _use_raw;

    /* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
public:
    explicit RMSMappingQuality(bool use_raw)
        : InterfaceInfoFieldAnnotation()
        , _use_raw(use_raw)
    {}

    void annotate(pRefFragment ref, pVariant vc, pRALikelihoods likelihoods, pInfoData target, pMemoryPool pool) override;

    ~RMSMappingQuality() override = default;
};

}  // namespace rovaca

#endif  // ROVACA_HC_RMS_MAPPING_QUALITY_H_
