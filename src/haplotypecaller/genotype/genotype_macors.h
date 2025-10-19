#ifndef ROVACA_HC_GENOTYPE_MACORS_H_
#define ROVACA_HC_GENOTYPE_MACORS_H_
#include <cmath>

namespace rovaca
{

#define un_used(x)                 ((void)x)

#define ROVACA_DELETE_FUNCTION(decl) decl = delete
#define DISALLOW_COPY(TypeName)    ROVACA_DELETE_FUNCTION(TypeName(const TypeName&))
#define DISALLOW_ASSIGN(TypeName)  ROVACA_DELETE_FUNCTION(void operator=(const TypeName&))
#define DISALLOW_COPY_AND_ASSIGN(TypeName)           \
    ROVACA_DELETE_FUNCTION(TypeName(const TypeName&)); \
    ROVACA_DELETE_FUNCTION(void operator=(const TypeName&))

#define INVALID_INT       (-1)
#define POSITIVE_INFINITY (INFINITY)
#define NEGATIVE_INFINITY (-INFINITY)

#define ROVACA_LIKELY(x)    __builtin_expect(!!(x), 1)
#define ROVACA_UNLIKELY(x)  __builtin_expect(!!(x), 0)

/*! @brief 状态值，表示未被初始化的ciagr_op */
#define BAM_CUNINITIALIZE (99)

#define MAX_GENOTYPE_QUAL (99)

#define HOM_GT_CUNT       (2)

/*! @brief 分配固定数量内存 */
#define ALLOC_MEM_IN_POOL(pool, size) (pool->allocate(size))
/*! @brief 分配type类型内存 */
#define ALLOC_TYPE_IN_POOL(pool, type) (pool->allocate(sizeof(type)))
/*! @brief 分配柔性数组内存 */
#define ALLOC_FLEXIBLE_IN_POOL(pool, type1, num, type2) (pool->allocate(sizeof(type1) + sizeof(type2) * (num + 1)))

}  // namespace rovaca

#endif  // ROVACA_HC_GENOTYPE_MACORS_H_
