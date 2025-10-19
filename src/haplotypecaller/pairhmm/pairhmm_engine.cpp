#include "pairhmm_engine.h"

#include "core/avx.h"
#include "core/avx_impl.h"
#include "rovaca_logger.h"
#include "pairhmm_internal.h"

float (*s_compute_full_prob_float_old_pairhmm)(const TestCase& tc) = nullptr;
double (*s_compute_full_prob_double_old_pairhmm)(const TestCase& tc) = nullptr;

namespace rovaca
{

void init_pairhmm_ptr(bool use_old);

DoubleVector2D (*call_pairhmm)(const HaplotypeVector& hs, ReadVector& rs, int32_t min_quality_threshold, PcrIndelModel pcr_option,
                               pMemoryPool pool);

bool avx512_supported() { return is_avx512_supported(); }

void init_pairhmm_ptr(bool use_old)
{
    bool avx512_supported = is_avx512_supported();
    bool avx2_supported = is_avx2_supported();

    if (avx512_supported) {
        RovacaLogger::info("pairhmm using avx512 instruction");
    }
    else if (avx2_supported) {
        RovacaLogger::info("pairhmm using avx2 instruction");
    }
    else {
        RovacaLogger::error("the machine must support avx2 or avx512 instruction");
        exit(EXIT_FAILURE);
    }

    if (use_old || !avx512_supported) {
        call_pairhmm = call_old_pairhmm;
    }
    else {
        call_pairhmm = call_new_pairhmm002;
    }

    ConvertChar::init();
    if (avx512_supported) {
        s_compute_full_prob_float_old_pairhmm = compute_fp_avx512_s;
        s_compute_full_prob_double_old_pairhmm = compute_fp_avx512_d;
    }
    else {
        s_compute_full_prob_float_old_pairhmm = compute_fp_avx2_s;
        s_compute_full_prob_double_old_pairhmm = compute_fp_avx2_d;
    }
}

}  // namespace rovaca