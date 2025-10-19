#include "quality_utils.h"

#include <gtest/gtest.h>

#include <cmath>

using namespace rovaca;
using namespace QualityUtils;
static constexpr double precision = 1e-6;

TEST(QualityUtilsUnitTest, testQualCaches)
{
    ASSERT_NEAR(qual_to_error_prob((uint8_t)20), 0.01, precision);
    ASSERT_NEAR(qual_to_error_prob_log10((uint8_t)20), std::log10(0.01), precision);
    ASSERT_NEAR(qual_to_prob((uint8_t)20), 0.99, precision);
    ASSERT_NEAR(qual_to_prob_log10((uint8_t)20), std::log10(0.99), precision);

    ASSERT_NEAR(qual_to_error_prob((uint8_t)30), 0.001, precision);
    ASSERT_NEAR(qual_to_error_prob_log10((uint8_t)30), std::log10(0.001), precision);
    ASSERT_NEAR(qual_to_prob((uint8_t)30), 0.999, precision);
    ASSERT_NEAR(qual_to_prob_log10((uint8_t)30), std::log10(0.999), precision);

    ASSERT_NEAR(qual_to_error_prob((uint8_t)40), 0.0001, precision);
    ASSERT_NEAR(qual_to_error_prob_log10((uint8_t)40), std::log10(0.0001), precision);
    ASSERT_NEAR(qual_to_prob((uint8_t)40), 0.9999, precision);
    ASSERT_NEAR(qual_to_prob_log10((uint8_t)40), std::log10(0.9999), precision);
}