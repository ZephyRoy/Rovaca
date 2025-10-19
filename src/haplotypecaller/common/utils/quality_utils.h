#ifndef ROVACA_HC_QUALITY_UTILS_H_
#define ROVACA_HC_QUALITY_UTILS_H_
#include <cstdint>

namespace rovaca
{

/*!
 * @brief 其中uint8_t版本使用cache
 */
namespace QualityUtils
{

static constexpr uint8_t s_mapping_quality_unavailable = 255;

/*! @brief convert a phred-scaled quality score to its probability of being wrong (q30 => 0.001) */
double qual_to_error_prob(double qual);
double qual_to_error_prob(uint8_t qual);

/*! @brief onvert a phred-scaled quality score to its log10 probability of being wrong (q30 => log10(0.001)) */
double qual_to_error_prob_log10(double qual);
double qual_to_error_prob_log10(uint8_t qual);

/*! @brief convert a phred-scaled quality score to its probability of being true (q30 => 0.999) */
double qual_to_prob(double qual);
double qual_to_prob(uint8_t qual);

/*! @brief convert a phred-scaled quality score to its log10 probability of being true (q30 => log10(0.999)) */
double qual_to_prob_log10(uint8_t qual);

/*! @brief convert a log-probability to a phred-scaled value  ( log(0.001) => 30 ) */
double log_prob_to_phred(double prob);

/*! @brief convert a probability of being wrong to a phred-scaled quality score as a double */
double phred_scale_error_rate(double error_rate);

/*! @brief convert a log10 probability of being wrong to a phred-scaled quality score as a double */
double phred_scale_log10error_rate(double error_rate_log10);

// 添加 init 方法，避免静态对象初始化顺序导致的错误
void init();

}  // namespace QualityUtils

}  // namespace rovaca

#endif  // ROVACA_HC_QUALITY_UTILS_H_