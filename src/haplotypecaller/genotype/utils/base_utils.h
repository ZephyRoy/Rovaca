#ifndef ROVACA_HC_BASE_UTILS_H_
#define ROVACA_HC_BASE_UTILS_H_
#include <random>

#include "forward.h"

namespace rovaca
{

namespace BaseUtils
{

bool is_regular_base(uint8_t byte);
bool is_all_regular_bases(pBases bases);

int32_t simple_base_to_base_index(uint8_t byte);

/*! @brief 使用整个bases_str构建pBases */
pBases bases_create(const char* bases_str, pMemoryPool pool);
/*! @brief 使用bases_str开始的num个字符构建pBases */
pBases bases_create(const char* bases_str, uint32_t num, pMemoryPool pool);

void fill_with_random_bases(pBases dest, int32_t start, int32_t end, std::mt19937& gen);

}  // namespace BaseUtils

}  // namespace rovaca

#endif  // ROVACA_HC_BASE_UTILS_H_
