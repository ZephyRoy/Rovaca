#include "haplotype.h"

#include <cstring>

#include "genotype_macors.h"

namespace rovaca
{

pHaplotype Haplotype::create(pMemoryPool pool) { return new ALLOC_TYPE_IN_POOL(pool, Haplotype) Haplotype{}; }

void Haplotype::init_haplotype(const char *bases_str, uint8_t is_ref, pMemoryPool pool)
{
    uint32_t num = strlen(bases_str) + 1;
    auto bases = new ALLOC_FLEXIBLE_IN_POOL(pool, Bases, num, uint8_t) Bases{num};
    memcpy(bases->data, bases_str, num * sizeof(uint8_t));
    return init_haplotype(bases_str, num, is_ref, pool);
}

void Haplotype::init_haplotype(const char *bases_str, uint32_t num, uint8_t is_ref, pMemoryPool pool)
{
    auto bases = new ALLOC_FLEXIBLE_IN_POOL(pool, Bases, num, uint8_t) Bases{num};
    memcpy(bases->data, bases_str, num * sizeof(uint8_t));
    init_allele(bases, is_ref);
}

}  // namespace rovaca