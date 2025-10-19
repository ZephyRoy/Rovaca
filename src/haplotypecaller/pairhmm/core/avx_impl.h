#ifndef AVX_IMPL_H_
#define AVX_IMPL_H_
#include "common.h"

float compute_fp_avx2_s(const TestCase &tc);

double compute_fp_avx2_d(const TestCase &tc);

float compute_fp_avx512_s(const TestCase &tc);

double compute_fp_avx512_d(const TestCase &tc);

#endif  // AVX_IMPL_H_
