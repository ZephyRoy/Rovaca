#ifndef AVX_512_FLOAT_H
#define AVX_512_FLOAT_H
#include <vector>

#include "common.h"

namespace rovaca
{

/**
 * \brief 计算一个tc中的结果: 1*n 或 n*n
 * \param tc 包含reads、单倍体，以及对应于result中的坐标
 * \param result 存储计算结果，由外部开辟空间
 */
void compute_full_prob_avx512_float_1xn(const TestCase& tc, float* result);

}  // namespace rovaca

#endif  // AVX_512_FLOAT_H
