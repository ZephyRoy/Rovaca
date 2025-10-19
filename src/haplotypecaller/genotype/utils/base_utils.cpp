#include "base_utils.h"

#include <algorithm>
#include <cstring>

#include "genotype_macors.h"
#include "genotype_struct.h"

namespace rovaca
{

static constexpr int32_t base_index_map[256] = {
    -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
    -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, 0,  -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
    -1, 0,  -1, 1,  -1, -1, -1, 2,  -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, 3,  -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
    -1, 0,  -1, 1,  -1, -1, -1, 2,  -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, 3,  -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
    -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
    -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
    -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
    -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1};

bool BaseUtils::is_regular_base(uint8_t byte) { return simple_base_to_base_index(byte) != -1; }

bool BaseUtils::is_all_regular_bases(pBases bases)
{
    return std::all_of(bases->data, bases->data + bases->num, [](uint8_t b) { return is_regular_base(b); });
}

int32_t BaseUtils::simple_base_to_base_index(uint8_t byte) { return base_index_map[byte & 0xff]; }

pBases BaseUtils::bases_create(const char *bases_str, pMemoryPool pool)
{
    uint32_t len = strlen(bases_str);
    auto bases = new ALLOC_FLEXIBLE_IN_POOL(pool, Bases, len, uint8_t) Bases{len};
    memcpy(bases->data, bases_str, len * sizeof(uint8_t));
    return bases;
}

pBases BaseUtils::bases_create(const char *bases_str, uint32_t num, pMemoryPool pool)
{
    pBases bases = new ALLOC_FLEXIBLE_IN_POOL(pool, Bases, num, uint8_t) Bases{num};
    memcpy(bases->data, bases_str, num * sizeof(uint8_t));
    return bases;
}

void BaseUtils::fill_with_random_bases(pBases dest, int32_t start, int32_t end, std::mt19937 &gen)
{
    static const char *bases_str = "ACGT";
    std::uniform_int_distribution<int32_t> dis(0, 3);
    for (int32_t i = start; i < end; ++i) {
        dest->data[i] = bases_str[dis(gen)];
    }
}

}  // namespace rovaca