#ifndef ROVACA_HC_BLOCK_COMBINER_H_
#define ROVACA_HC_BLOCK_COMBINER_H_
#include <algorithm>
#include <cstddef>
#include <cstdint>

#include "forward.h"
#include "htslib/vcf.h"

namespace rovaca
{

class BlockCombiner
{
    /* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
private:
    bool floor_blocks_;
    bool is_gvcf_model_;
    bcf_hdr_t* vcf_herder_;
    bam_hdr_t* bam_herder_;
    std::map<int32_t, std::pair<int32_t, int32_t>> gq_partitions_;

    bcf1_t* output_;
    int32_t sample_id_;
    int32_t contig_of_next_available_start_;
    pHomRefBlock current_block_;
    int64_t next_available_start_;

    size_t del_line_offset_;
    int64_t next_variant_pos_;

    /* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
public:
    BlockCombiner(const std::vector<int32_t>& gq_partitions, bcf_hdr_t* vcf_herder, bam_hdr_t* bam_herder, bool floor_blocks, bool is_gvcf);
    ~BlockCombiner();

    /*!
     * @brief 添加 vc
     * @param vc
     */
    void submit(pVariant vc, kstring_t* s);
    void submit_vector(const VariantVector& vcs, kstring_t* s)
    {
        std::for_each(vcs.begin(), vcs.end(), [&](pVariant vc) { submit(vc, s); });
    }

    /*!
     * @brief 一个 task 完成后, 调用 force_output()
     */
    void force_output(kstring_t* s);

    size_t del_line_offset() const { return del_line_offset_; }
    int64_t next_variant_pos() const { return next_variant_pos_; }

    /* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
private:
    void clear_combiner();

    static std::map<int32_t, std::pair<int32_t, int32_t>> parse_partitions(const std::vector<int32_t>& gq_partitions);

    void emit_current_block(kstring_t* s);

    bcf1_t* add_hom_ref_site(pVariant vc, pGenotype g);

    bool genotype_can_be_merged_in_current_block(pGenotype g);

    pHomRefBlock create_new_block(pVariant vc, pGenotype g);
};

}  // namespace rovaca

#endif  // ROVACA_HC_BLOCK_COMBINER_H_
