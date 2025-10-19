#ifndef ROVACA_HC_ADAPTER_UTILS_H_
#define ROVACA_HC_ADAPTER_UTILS_H_
#include "../forward.h"
#include "htslib/vcf.h"

namespace rovaca
{

namespace AdapterUtils
{

/*!
 * @brief 有几个 interval:
 *      original: 原始的interval
 *      original_padded: original 左右扩 100
 *
 *      variant: 由 original 根据变异位点trim出的包含所有变异位点的最小interval
 *      padded_variant: variant 左右扩 20
 *
 *      left_flank: gvcf 模式下左侧无变异区间 [original->get_start(), variant->get_start()]
 *      left_flank_padded: left_flank 左右扩 100
 *
 *      right_flank: gvcf 模式下右侧无变异区间 [variant->get_stop() + 1, original->get_stop()]
 *      right_flank_padded: right_flank 左右扩 100
 *
 *      new_span: original->intersect(left_flank 或 right_flank)
 *      new_span_padded: padded->intersect(left_flank_padded 或 right_flank_padded)
 * @param original
 * @param variant
 * @param padding
 * @param chr_len
 * @param pool
 * @return
 */
IntervalPair non_variant_left_flank_region(pSimpleInterval original, pSimpleInterval original_padded, pSimpleInterval variant,
                                           int64_t padding, int64_t chr_len, pMemoryPool pool);
IntervalPair non_variant_right_flank_region(pSimpleInterval original, pSimpleInterval original_padded, pSimpleInterval variant,
                                            int64_t padding, int64_t chr_len, pMemoryPool pool);

/*!
 * @brief 根据组装出来的单倍体，计算 all_variation_events(包含了 get_variation_events).
 * 筛选出与 original 重叠的变异体，并计算它们的起始位置和结束位置的最小值和最大值创建一个新的SimpleInterval对象，得到 variant
 * 根据 all_variation_events 是否存在 indel 计算 padding 大小，得到 padded_variant
 * @param haplotypes
 * @param ref
 * @param ref_loc
 * @param original
 * @param original_padded
 * @param args
 * @param pool
 * @return
 */
IntervalPair trim_region(HaplotypeVector& haplotypes, pRefFragment ref, pSimpleInterval ref_loc, pSimpleInterval original,
                         pSimpleInterval original_padded, pHCArgs args, pMemoryPool pool);

/*!
 * @brief 根据 trim_region 得到的 padded_variant 裁剪 reads
 * @param reads
 * @param padded
 * @param pool
 * @param bpool
 * @return
 */
ReadHashSet trim_reads_by_region(const ReadHashSet& reads, pSimpleInterval padded, pMemoryPool pool, pBamDataPool bpool,
                                 std::pmr::list<bam1_t*>& extra_memory_reads);

/*!
 * @brief 根据 trim_region 得到的 padded_variant 裁剪 haplotype
 * @param haplotypes
 * @param padded
 * @param pool
 * @return
 */
HaplotypeVector trim_haplotype_by_region(const HaplotypeVector& haplotypes, pSimpleInterval padded, pMemoryPool pool);

/*!
 * @brief 针对 trim 后的 read进行的第一次过滤，此次过滤出来的reads不再需要
 * @param reads trimed reads
 * @return
 */
ReadList filter_non_passing_reads1(ReadHashSet& reads, int32_t minimum_read_length_after_trimming, pMemoryPool pool);

/*!
 * @brief 二次过滤，保留过滤的 read 后续使用
 * @param reads
 * @param pool
 * @return
 */
ReadList filter_non_passing_reads2(ReadHashSet& reads, uint8_t mapping_quality_threshold, pMemoryPool pool);

void vcf_info2bcf(pInfoData info, bcf1_t* rec);
void gvcf_info2bcf(pInfoData info, bcf1_t* rec);

/*!
 * @brief 变异点转bcf1_t
 */
bcf1_t* variant2bcf(bam_hdr_t* bam_header, bcf_hdr_t* vcf_header, pVariant vc, bcf1_t* b);

/*!
 * @brief 非变异点转bcf1_t
 */
bcf1_t* variant2bcf(bam_hdr_t* bam_header, bcf_hdr_t* vcf_header, int32_t tid, int64_t pos, pAllele ref, int32_t block_end, int32_t dp,
                    int32_t gq, int32_t min_dp, const std::vector<int32_t>& pls, bcf1_t* b);

}  // namespace AdapterUtils

}  // namespace rovaca

#endif  // ROVACA_HC_ADAPTER_UTILS_H_
