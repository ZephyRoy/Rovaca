#include "quality_utils.h"

#include <array>
#include <cmath>

namespace rovaca
{

const double s_phred_to_log_prob_multiplier = -std::log10(10) / 10.0;
const double s_log_prob_to_phred_multiplier = 1.0 / s_phred_to_log_prob_multiplier;
const double s_min_log10_scaled_qual = std::log10(4.9E-324);

/*! @brief maximum sense quality value */
static constexpr int32_t s_max_qual = 254;

std::array<double, s_max_qual + 1> s_qual_to_error_prob_cache{};
std::array<double, s_max_qual + 1> s_qual_to_prob_log10cache{};

double QualityUtils::qual_to_error_prob(double qual) { return std::pow(10.0, qual / -10.0); }

double QualityUtils::qual_to_error_prob(uint8_t qual) { return s_qual_to_error_prob_cache.at(qual & 0xff); }

double QualityUtils::qual_to_error_prob_log10(double qual) { return qual * -0.1; }

double QualityUtils::qual_to_error_prob_log10(uint8_t qual) { return qual_to_error_prob_log10(double(qual & 0xff)); }

double QualityUtils::qual_to_prob(double qual) { return 1.0 - qual_to_error_prob(qual); }

double QualityUtils::qual_to_prob(uint8_t qual) { return 1.0 - qual_to_error_prob(qual); }

double QualityUtils::qual_to_prob_log10(uint8_t qual) { return s_qual_to_prob_log10cache.at(qual & 0xff); }

double QualityUtils::log_prob_to_phred(double prob) { return prob * s_log_prob_to_phred_multiplier; }

double QualityUtils::phred_scale_error_rate(double error_rate) { return phred_scale_log10error_rate(std::log10(error_rate)); }

double QualityUtils::phred_scale_log10error_rate(double error_rate_log10)
{
    return std::abs(-10.0 * std::max(error_rate_log10, s_min_log10_scaled_qual));
}

void QualityUtils::init()
{
    static bool init = false;

    if (!init) {
        for (int32_t i = 0; i <= s_max_qual; ++i) {
            s_qual_to_error_prob_cache[i] = std::pow(10.0, (double)i / -10.0);
        }

        for (int32_t i = 0; i <= s_max_qual; ++i) {
            s_qual_to_prob_log10cache[i] = std::log10(1.0 - s_qual_to_error_prob_cache.at(i));
        }

        init = true;
    }
}

}  // namespace rovaca