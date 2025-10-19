#ifndef COMMON_H
#define COMMON_H
#include <cstdint>
#include <vector>

namespace rovaca
{

#define MIN_ACCEPTED 1e-28f

#define ALIGNED      __attribute__((aligned(64)))

constexpr uint8_t k_dummy = 'D';

struct ALIGNED TestCase
{
    uint64_t read_num : 32;
    uint64_t hap_num : 32;
    uint64_t min_read_len : 32;
    uint64_t min_hap_len : 32;
    uint64_t max_read_len : 32;
    uint64_t max_hap_len : 32;
    ALIGNED uint32_t *read_len_arr;
    ALIGNED uint32_t *hap_len_arr;
    ALIGNED int32_t *reads;
    ALIGNED int32_t *haps;
    ALIGNED float *distm;
    ALIGNED float *_1_distm;
    ALIGNED float *gapm;
    ALIGNED float *mm;
    ALIGNED float *mi;
    ALIGNED float *ii;
    ALIGNED float *md;
    ALIGNED float *dd;
};

struct DebugData
{
    uint32_t read_idx, hap_idx;
    float dist;
    float m[11];
    float i[7];
    float d[7];
};

}  // namespace rovaca

#endif  // COMMON_H
