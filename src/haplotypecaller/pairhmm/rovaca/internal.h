#ifndef INTERNAL_H
#define INTERNAL_H
#include <cstdint>

namespace rovaca
{

constexpr float k_min_accepted = 1e-28f;
constexpr float l_log10_initial_constant_f = 36.1236000061;         // log10f(ldexpf(1.f, 120.f))
constexpr double l_log10_initial_constant_d = 307.050595577260822;  // log10(ldexp(1.0, 1020.0))

constexpr uint32_t k_avx256_bit_count = 256;
constexpr uint32_t k_avx512_bit_count = 512;

}  // namespace rovaca

#endif  // INTERNAL_H
