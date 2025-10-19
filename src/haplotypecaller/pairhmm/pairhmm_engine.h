#ifndef __PAIRHMM_ENGINE_H__
#define __PAIRHMM_ENGINE_H__

#include "../common/enum.h"
#include "../genotype/forward.h"

namespace rovaca
{

extern bool avx512_supported();

extern void init_pairhmm_ptr(bool use_old);

extern DoubleVector2D (*call_pairhmm)(const HaplotypeVector& hs, ReadVector& rs, int32_t min_quality_threshold, PcrIndelModel pcr_option,
                                      pMemoryPool pool);

}  // namespace rovaca

#endif  // __PAIRHMM_ENGINE_H__