#ifndef ROVACA_HC_HAPLOTYPE_UTILS_H_
#define ROVACA_HC_HAPLOTYPE_UTILS_H_
#include "forward.h"

namespace rovaca
{

namespace HaplotypeUtils
{

/*! @brief 创建一个新的单倍型，从这个单倍型派生，正好跨越所提供的位置 */
pHaplotype trim(pHaplotype h, pSimpleInterval interval, bool ignore_ref_state, pMemoryPool pool);
pHaplotype trim(pHaplotype h, pSimpleInterval interval, pMemoryPool pool);

}  // namespace HaplotypeUtils

}  // namespace rovaca

#endif  // ROVACA_HC_HAPLOTYPE_UTILS_H_
