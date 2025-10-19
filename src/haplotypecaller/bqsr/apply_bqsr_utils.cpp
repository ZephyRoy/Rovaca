#include "apply_bqsr_utils.h"

#include <float.h>

#include <algorithm>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <limits>
#include <string>

static const double TWO_PI = 2 * 3.141592653;
static const double HALF_LOG_2_PI = 0.5 * std::log(TWO_PI);

static double g_log10_qemp_prior[MAX_USABLE_Q_SCORE + 1] = {
    -0.22579135264472738, -2.2257913526447273, -8.225791352644727,  -18.225791352644727, -32.22579135264473,  -50.22579135264473,
    -72.22579135264473,   -98.22579135264473,  -128.22579135264473, -162.22579135264473, -200.22579135264473, -242.22579135264473,
    -288.2257913526447,   -338.2257913526447,  -392.2257913526447,  -450.2257913526447,  -512.2257913526447,  -578.2257913526447,
    -648.2257913526447,   -722.2257913526447,  -800.2257913526447,  -882.2257913526447,  -968.2257913526447,  -1058.2257913526448,
    -1152.2257913526448,  -1250.2257913526448, -1352.2257913526448, -1458.2257913526448, -1568.2257913526448, -1682.2257913526448,
    -1800.2257913526448,  -1922.2257913526448, -2048.2257913526446, -2178.2257913526446, -2312.2257913526446, -2450.2257913526446,
    -2592.2257913526446,  -2738.2257913526446, -2888.2257913526446, -3042.2257913526446, -3200.2257913526446};

static double EXACT_STIRLING_ERRORS[128] = {
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

static double getDeviancePart(double x, double mu)
{
    double ret;
    if (std::abs(x - mu) < 0.1 * (x + mu)) {
        double d = x - mu;
        double v = d / (x + mu);
        double s1 = v * d;
        double s = std::numeric_limits<double>::quiet_NaN();
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

static double getStirlingError(double z)
{
    double ret;
    if (z < 15.0) {
        double z2 = 2.0 * z;
        if (std::floor(z2) == z2) {
            ret = EXACT_STIRLING_ERRORS[(int)z2];
        }
        else {
            ret = std::lgamma(z + 1.0) - (z + 0.5) * std::log(z) + z - HALF_LOG_2_PI;
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

static double logBinomialProbability(int x, int n, double p, double q)
{
    double ret;
    if (x == 0) {
        if (p < 0.1) {
            ret = -getDeviancePart(n, n * q) - n * p;
        }
        else {
            ret = n * std::log(q);
        }
    }
    else if (x == n) {
        if (q < 0.1) {
            ret = -getDeviancePart(n, n * p) - n * q;
        }
        else {
            ret = n * std::log(p);
        }
    }
    else {
        ret =
            getStirlingError(n) - getStirlingError(x) - getStirlingError(n - x) - getDeviancePart(x, n * p) - getDeviancePart(n - x, n * q);
        double f = (TWO_PI * x * (n - x)) / n;
        ret = -0.5 * std::log(f) + ret;
    }
    return ret;
}

double log10_binomial_probability(int n, int k, double p)
{
    if (n == 0) {
        return (k == 0) ? 0.0 : -DBL_MAX;
    }
    double ret = 0.0;
    if (k < 0 || k > n) {
        ret = -DBL_MAX;
    }
    else {
        ret = logBinomialProbability(k, n, p, 1.0 - p);
    }
    return ret;
}

double log_binomial_coeffcient(int n, int k)
{
    if (n < 0) {}
    if (k > n || k < 0) {}
    return log10_factorial(n) - log10_factorial(k) - log10_factorial(n - k);
}

double log10_factorial(int n)
{
    if (n < MAX_LGGAMMA_CACHE) return g_log10_factorial_cache[n];
    return std::lgamma(n + 1) * LOG10_E;
}

void clip_low_qual_ends_in_method_WRITE_NS(bam1_t* read, uint8_t* clipped_bases, int start, int stop)
{
    int i = 0;
    int end = std::min(read->core.l_qseq, stop + 1);
    for (i = start; i < end; ++i) {
        clipped_bases[i] = 'N';
    }
}

// Initialize complement array
uint8_t init_complement[256] = {0};
void initialize_complement()
{
    init_complement['A'] = 'T';
    init_complement['C'] = 'G';
    init_complement['G'] = 'C';
    init_complement['T'] = 'A';
    init_complement['a'] = 'T';
    init_complement['c'] = 'G';
    init_complement['g'] = 'C';
    init_complement['t'] = 'A';
}

void reverse_bam_seq(bam1_t* read, uint8_t* clipped_bases)
{
    uint8_t* seq = bam_get_seq(read);
    int32_t read_len = read->core.l_qseq;

    static bool is_initialized = false;
    if (!is_initialized) {
        initialize_complement();
        is_initialized = true;
    }

    for (int i = 0; i < read_len; i++) {
        uint8_t base = seq_nt16_str[bam_seqi(seq, read_len - 1 - i)];
        clipped_bases[i] = init_complement[base] ? init_complement[base] : base;  // Use complement if available, else use base
    }
}

#include <stdexcept>
#include <string>

int parseInt(const std::string& s, int radix)
{
    if (s.empty()) {
        throw std::invalid_argument("empty string");
    }
    else if (radix < 2) {
        throw std::invalid_argument("radix less than 2");
    }
    else if (radix > 36) {
        throw std::invalid_argument("radix greater than 36");
    }
    else {
        bool negative = false;
        int i = 0;
        int len = s.length();
        int limit = -2147483647;

        if (len <= 0) {
            throw std::invalid_argument("invalid input string");
        }
        else {
            char firstChar = s[0];

            if (firstChar < '0') {
                if (firstChar == '-') {
                    negative = true;
                    limit = INT_MIN;
                }
                else if (firstChar != '+') {
                    throw std::invalid_argument("invalid input string");
                }

                if (len == 1) {
                    throw std::invalid_argument("invalid input string");
                }

                ++i;
            }

            int multmin = limit / radix;
            int result = 0;
            int digit;

            while (i < len) {
                digit = s[i++];
                if (isdigit(digit)) {
                    digit -= '0';
                }
                else if (isalpha(digit)) {
                    digit = tolower(digit) - 'a' + 10;
                }
                else {
                    throw std::invalid_argument("invalid input string");
                }

                if (digit < 0 || result < multmin) {
                    throw std::invalid_argument("invalid input string");
                }

                result *= radix;

                if (result < limit + digit) {
                    throw std::invalid_argument("invalid input string");
                }

                result -= digit;
            }

            return negative ? result : -result;
        }
    }
}

double convert_to_double(const std::string& str, int n)
{
    double value = std::stod(str);
    (void)n;
    return value;
}

int bqsr_covariate_key_from_context(uint8_t* dna, int start, int end)
{
    int key = end - start;
    int bit_offset = LENGTH_BITS;
    int base_index, i;

    for (i = start; i < end; i++) {
        base_index = g_bqsr_seq_table[dna[i]];
        if (base_index == -1) {  // ignore non-ACGT bases
            return -1;
        }
        key |= (base_index << bit_offset);
        bit_offset += 2;
    }
    return key;
}

double empirical_quality_bayesian_estimate(long n_observations, long n_errors, double Q_reported)
{
    int num_bins = (MAX_REASONABLE_Q_SCORE + 1) * RESOLUTION_BINS_PER_QUAL;
    double QEmpOfBin, log10_posteriors, log10_posteriors_max = -DBL_MAX;
    int max_index = 0;
    for (int i = 0; i < num_bins; ++i) {
        QEmpOfBin = (double)i;
        log10_posteriors = log10_qual_emp_prior(QEmpOfBin, Q_reported) + log10_qual_emp_likelihood(QEmpOfBin, n_observations, n_errors);
        if (log10_posteriors > log10_posteriors_max) {
            log10_posteriors_max = log10_posteriors;
            max_index = i;
        }
    }
    return (double)max_index;
}

double log10_qual_emp_prior(double Qempirical, double Qreported)
{
    const int difference = std::min(std::abs(static_cast<int>(Qempirical - Qreported)), MAX_USABLE_Q_SCORE);
    return g_log10_qemp_prior[difference];
}

double log10_qual_emp_likelihood(double Qempirical, long nObservations, long nErrors)
{
    if (nObservations == 0) return 0.0;
    if (nObservations > MAX_NUMBER_OF_OBSERVATIONS) {
        // we need to decrease nErrors by the same fraction that we are decreasing nObservations
        double fraction = (double)MAX_NUMBER_OF_OBSERVATIONS / (double)nObservations;
        nErrors = std::round((double)nErrors * fraction);
        nObservations = MAX_NUMBER_OF_OBSERVATIONS;
    }
    // this is just a straight binomial PDF
    double log10Prob = log10_binomial_probability((int)nObservations, (int)nErrors, std::pow(10.0, Qempirical * (-0.1)));
    if (std::isinf(log10Prob) || std::isnan(log10Prob)) log10Prob = -DBL_MAX;
    return log10Prob;
}