#ifndef ROVACA_HC_ALIGNMENT_UTILS_H_
#define ROVACA_HC_ALIGNMENT_UTILS_H_
#include "cigar_builder.h"
#include "forward.h"
#include "index_range.hpp"

namespace rovaca
{

namespace AlignmentUtils
{

static constexpr uint8_t s_gap_qual_character{0};
static constexpr uint8_t s_gap_base_character{45};

typedef class CigarPairTransform CigarPairTransform, *pCigarPairTransform;
typedef class CigarPairTransformManager CigarPairTransformManager, *pCigarPairTransformManager;

uint32_t length_on_read(uint32_t cigar_elenment);
uint32_t length_on_reference(uint32_t cigar_elenment);

int64_t read_start_on_reference_haplotype(pCigar haplotype_vs_ref_cigar, int64_t read_start_on_haplotype);

pCigar trim_cigar_by_bases(pCigar cigar, int64_t start, int64_t end, pMemoryPool pool);

pCigar trim_cigar_by_reference(pCigar cigar, int64_t start, int64_t end, pMemoryPool pool);

pCigar trim_cigar(pCigar cigar, int64_t start, int64_t end, bool by_reference, pMemoryPool pool);

pCigar apply_cigar_to_cigar(pCigar first_to_second, pCigar second_to_third, pMemoryPool pool);

pCigarPairTransform get_transformer(uint32_t op12, uint32_t op23);

bool last_base_on_right_is_same(const std::pmr::vector<uint8_t*>& sequences, std::pmr::vector<IndexRange>& bounds);

bool first_base_on_left_is_same(const std::pmr::vector<uint8_t*>& sequences, std::pmr::vector<IndexRange>& bounds);

bool next_base_on_left_is_same(const std::pmr::vector<uint8_t*>& sequences, std::pmr::vector<IndexRange>& bounds);

std::pair<int32_t, int32_t> normalize_alleles(const std::pmr::vector<uint8_t*>& sequences, std::pmr::vector<IndexRange>& bounds,
                                              int32_t max_shift, bool trim);

CigarBuilder::Result left_align_indels(pCigar cigar, pBases ref_bases, pBases read_bases, int64_t read_start, pMemoryPool pool);

pCigar append_clipped_elements_from_cigar_to_cigar(pCigar cigar_to_have_clipped_elements_added, const uint32_t* original_clipped_cigar,
                                                   uint32_t cigar_num, pMemoryPool pool);

void create_read_aligned_to_ref(pReadRecord original_read, pHaplotype haplotype, pHaplotype ref_haplotype, int64_t reference_start,
                                bool is_informative, p_lib_sw_avx sw, pMemoryPool pool, pBamDataPool bam_pool);

/*!
 * @brief 从给定的碱基序列和CIGAR字符串中获取覆盖参考区间的碱基序列
 * @note 参考区间起始和结束位置是0-based的
 */
pBases get_bases_covering_ref_interval(int64_t ref_start, int64_t ref_end, pBases bases, int64_t bases_start_on_ref,
                                       pCigar bases_to_ref_cigar, pMemoryPool pool);

/*!
 * @brief 将 Read 序列按照其 cigar 字符串对齐到 Ref 序列上，并返回对齐后的碱基和质量信息。
 *        Read 序列中存在但 Ref 序列中不存在的碱基，则会被丢弃。
 *        Ref 序列中存在但 Read 序列中不存在的碱基，则会用GAP_CHARACTER值填充 Read 序列中的碱基，并用0填充质量信息。
 */
std::pair<pBases, pBases> get_bases_and_base_qualities_aligned_one_to_one(pReadRecord read, pMemoryPool pool);
std::pair<pBases, pBases> get_bases_and_base_qualities_aligned_one_to_one(pReadRecord read, uint8_t base_pad, uint8_t qual_pad,
                                                                          pMemoryPool pool);

}  // namespace AlignmentUtils

}  // namespace rovaca

#endif  // ROVACA_HC_ALIGNMENT_UTILS_H_