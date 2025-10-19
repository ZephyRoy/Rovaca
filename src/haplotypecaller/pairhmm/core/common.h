#ifndef COMMON_H_
#define COMMON_H_
#include <x86intrin.h>

#include <cassert>
#include <cstdint>
#include <cstring>

#define ALIGNED32   __attribute__((aligned(32)))
#define likely(x)   __builtin_expect(!!(x), 1)
#define unlikely(x) __builtin_expect(!!(x), 0)

typedef struct
{
    uint32_t rslen;
    uint32_t haplen;
    const ALIGNED32 uint8_t* ALIGNED32 hap;
    const ALIGNED32 uint8_t* ALIGNED32 rs;
    const ALIGNED32 uint8_t* ALIGNED32 q;
    const ALIGNED32 uint8_t* ALIGNED32 i;
    const ALIGNED32 uint8_t* ALIGNED32 d;
    const ALIGNED32 uint8_t* ALIGNED32 c;
    ALIGNED32 void* ALIGNED32 debug_data;
} TestCase;

constexpr uint32_t cache_len = 'T' - 'A' + 1;

struct ConvertChar
{
public:
    static inline uint8_t k_conversion_table[cache_len];

    static void init()
    {
        // memset(k_conversion_table, '0', 255);

        k_conversion_table['A' - 'A'] = 0;
        k_conversion_table['C' - 'A'] = 1;
        k_conversion_table['T' - 'A'] = 2;
        k_conversion_table['G' - 'A'] = 3;
        k_conversion_table['N' - 'A'] = 4;
    }

    static inline uint8_t get(uint8_t input) { return k_conversion_table[input - 'A']; }
};

// avx2_s
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
typedef struct
{
    ALIGNED32 __m256* ALIGNED32 p_mm_arr;
    ALIGNED32 __m256* ALIGNED32 p_gapm_arr;
    ALIGNED32 __m256* ALIGNED32 p_mx_arr;
    ALIGNED32 __m256* ALIGNED32 p_xx_arr;
    ALIGNED32 __m256* ALIGNED32 p_my_arr;
    ALIGNED32 __m256* ALIGNED32 p_yy_arr;
} ProbalityPtrDataAvx2S;

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

// avx2_d
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
typedef struct
{
    ALIGNED32 __m256d* ALIGNED32 p_mm_arr;
    ALIGNED32 __m256d* ALIGNED32 p_gapm_arr;
    ALIGNED32 __m256d* ALIGNED32 p_mx_arr;
    ALIGNED32 __m256d* ALIGNED32 p_xx_arr;
    ALIGNED32 __m256d* ALIGNED32 p_my_arr;
    ALIGNED32 __m256d* ALIGNED32 p_yy_arr;
} ProbalityPtrDataAvx2D;

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

// avx512_s
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
typedef struct
{
    ALIGNED32 __m512* ALIGNED32 p_mm_arr;
    ALIGNED32 __m512* ALIGNED32 p_gapm_arr;
    ALIGNED32 __m512* ALIGNED32 p_mx_arr;
    ALIGNED32 __m512* ALIGNED32 p_xx_arr;
    ALIGNED32 __m512* ALIGNED32 p_my_arr;
    ALIGNED32 __m512* ALIGNED32 p_yy_arr;
} ProbalityPtrDataAvx512S;

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

// avx512_d
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
typedef struct
{
    ALIGNED32 __m512d* ALIGNED32 p_mm_arr;
    ALIGNED32 __m512d* ALIGNED32 p_gapm_arr;
    ALIGNED32 __m512d* ALIGNED32 p_mx_arr;
    ALIGNED32 __m512d* ALIGNED32 p_xx_arr;
    ALIGNED32 __m512d* ALIGNED32 p_my_arr;
    ALIGNED32 __m512d* ALIGNED32 p_yy_arr;
} ProbalityPtrDataAvx512D;

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#endif  // COMMON_H_
