#include "saddle_point_expansion.h"

#include "math_utils.h"

namespace rovaca
{

namespace SaddlePointExpansion
{

static const double PI = std::acos(-1.0);
static const double TWO_PI = PI * 2;
static const double s_half_log_2_pi = 0.5 * std::log(TWO_PI);
static const double s_exact_stirling_errors[31] = {
    0.0,                           /* 0.0 */
    0.1534264097200273452913848,   /* 0.5 */
    0.0810614667953272582196702,   /* 1.0 */
    0.0548141210519176538961390,   /* 1.5 */
    0.0413406959554092940938221,   /* 2.0 */
    0.03316287351993628748511048,  /* 2.5 */
    0.02767792568499833914878929,  /* 3.0 */
    0.02374616365629749597132920,  /* 3.5 */
    0.02079067210376509311152277,  /* 4.0 */
    0.01848845053267318523077934,  /* 4.5 */
    0.01664469118982119216319487,  /* 5.0 */
    0.01513497322191737887351255,  /* 5.5 */
    0.01387612882307074799874573,  /* 6.0 */
    0.01281046524292022692424986,  /* 6.5 */
    0.01189670994589177009505572,  /* 7.0 */
    0.01110455975820691732662991,  /* 7.5 */
    0.010411265261972096497478567, /* 8.0 */
    0.009799416126158803298389475, /* 8.5 */
    0.009255462182712732917728637, /* 9.0 */
    0.008768700134139385462952823, /* 9.5 */
    0.008330563433362871256469318, /* 10.0 */
    0.007934114564314020547248100, /* 10.5 */
    0.007573675487951840794972024, /* 11.0 */
    0.007244554301320383179543912, /* 11.5 */
    0.006942840107209529865664152, /* 12.0 */
    0.006665247032707682442354394, /* 12.5 */
    0.006408994188004207068439631, /* 13.0 */
    0.006171712263039457647532867, /* 13.5 */
    0.005951370112758847735624416, /* 14.0 */
    0.005746216513010115682023589, /* 14.5 */
    0.005554733551962801371038690  /* 15.0 */
};

double get_stirling_error(double z)
{
    double ret;
    if (z < 15.0) {
        double z2 = 2.0 * z;
        if (std::floor(z2) == z2) {
            ret = s_exact_stirling_errors[(int32_t)z2];
        }
        else {
            ret = std::lgamma(z + 1.0) - (z + 0.5) * std::log(z) + z - s_half_log_2_pi;
        }
    }
    else {
        double z2 = z * z;
        ret = (0.083333333333333333333 -
               (0.00277777777777777777778 -
                (0.00079365079365079365079365 - (0.000595238095238095238095238 - 0.0008417508417508417508417508 / z2) / z2) / z2) /
                   z2) /
              z;
    }
    return ret;
}

double get_deviance_part(double x, double mu)
{
    double ret;
    if (std::abs(x - mu) < 0.1 * (x + mu)) {
        double d = x - mu;
        double v = d / (x + mu);
        double s1 = v * d;
        double s = NAN;
        double ej = 2.0 * x * v;
        v *= v;
        int j = 1;
        while (s1 != s) {
            s = s1;
            ej *= v;
            s1 = s + ej / ((j * 2) + 1);
            ++j;
        }
        ret = s1;
    }
    else {
        ret = x * std::log(x / mu) + mu - x;
    }
    return ret;
}

double log_binomial_probability(int x, int n, double p, double q)
{
    double ret;
    if (x == 0) {
        if (p < 0.1) {
            ret = -get_deviance_part(n, n * q) - n * p;
        }
        else {
            ret = n * std::log(q);
        }
    }
    else if (x == n) {
        if (q < 0.1) {
            ret = -get_deviance_part(n, n * p) - n * q;
        }
        else {
            ret = n * std::log(p);
        }
    }
    else {
        ret = get_stirling_error(n) - get_stirling_error(x) - get_stirling_error(n - x) - get_deviance_part(x, n * p) -
              get_deviance_part(n - x, n * q);
        double f = (TWO_PI * x * (n - x)) / n;
        ret = -0.5 * std::log(f) + ret;
    }
    return ret;
}

}  // namespace SaddlePointExpansion

}  // namespace rovaca