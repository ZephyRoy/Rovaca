#ifndef ROVACA_HC_GENOTYPE_STRUCT_H_
#define ROVACA_HC_GENOTYPE_STRUCT_H_
#include <cstdint>

namespace rovaca
{

static constexpr std::size_t s_info_buffer_size = 512;  // 具体大小待议

/*! @brief Allele和Haplotype中字符串的存储结构，整个Bases对象源于pMemoryPool */
typedef struct Bases
{
    uint32_t num;
    uint8_t data[];
} Bases, *pBases;

/*! @brief Fasta片段 */
typedef struct RefFragment
{
    uint32_t len;
    uint8_t* data;
} RefFragment, *pRefFragment;

/*! @brief Haplotype中的cigar，整个Cigar对象源于pMemoryPool */
typedef struct Cigar
{
    uint32_t num;
    uint32_t data[];
} Cigar, *pCigar;

}  // namespace rovaca

#endif  // ROVACA_HC_GENOTYPE_STRUCT_H_
