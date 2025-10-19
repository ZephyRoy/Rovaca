#include "block_combiner.h"

#include <algorithm>
#include <cstdint>

#include "forward.h"
#include "genotype.h"
#include "genotype_macors.h"
#include "genotypes_context.hpp"
#include "hom_ref_block.h"
#include "rovaca_logger.h"
#include "utils/adapter_utils.h"
#include "variant.h"

namespace rovaca
{

BlockCombiner::BlockCombiner(const std::vector<int32_t> &gq_partitions, bcf_hdr_t *vcf_herder, bam_hdr_t *bam_herder, bool floor_blocks,
                             bool is_gvcf)
    : floor_blocks_(floor_blocks)
    , is_gvcf_model_(is_gvcf)
    , vcf_herder_(vcf_herder)
    , bam_herder_(bam_herder)
    , gq_partitions_(parse_partitions(gq_partitions))
    , output_{bcf_init()}
    , sample_id_(INVALID_INT)
    , contig_of_next_available_start_(INVALID_INT)
    , current_block_(nullptr)
    , next_available_start_(INVALID_INT)
    , del_line_offset_(0)
    , next_variant_pos_(INVALID_INT)
{}

BlockCombiner::~BlockCombiner() { bcf_destroy1(output_); }

void BlockCombiner::submit(pVariant vc, kstring_t *s)
{
    int32_t ret;
    uint32_t ref_length = vc->ref_allele()->length();
    bool is_del = ref_length > 1;
    bool is_variant = vc->has_log10_error();
    if (is_variant && is_del) {
        del_line_offset_ = 1;
        next_variant_pos_ = vc->get_start() + ref_length;
    }
    else if (del_line_offset_ && ((is_variant && !is_del) || vc->get_start() >= next_variant_pos_)) {
        del_line_offset_ = 0;
        next_variant_pos_ = INVALID_INT;
    }

    // 此处考虑添加小于 next_variant_pos_ 直接返回的逻辑

    if (!is_gvcf_model_) {
        AdapterUtils::variant2bcf(bam_herder_, vcf_herder_, vc, output_);
        ret = vcf_format1(vcf_herder_, output_, s);
        CHECK_CONDITION_EXIT(ret != 0, "error: vcf_format1");
    }
    else {
        if (sample_id_ == INVALID_INT) {
            sample_id_ = vc->genotype()->at(0)->sample_id();
        }

        if (current_block_ != nullptr && !current_block_->is_contiguous(vc)) {
            emit_current_block(s);
        }

        pAllele non_ref_allele = StaticAllele::get_instance()->_non_ref_allele.get();
        pGenotype g = vc->genotype()->at(0);
        if ((g->is_hom_ref() || (g->is_no_call() && g->has_likelihoods() && g->pl()[0] == 0)) &&
            vc->has_allele(non_ref_allele, false, false) && vc->is_biallelic()) {
            bcf1_t *maybe_completed_band = add_hom_ref_site(vc, g);
            if (maybe_completed_band) {
                ret = vcf_format1(vcf_herder_, output_, s);
                CHECK_CONDITION_EXIT(ret != 0, "error: vcf_format1");
            }
        }
        else {
            emit_current_block(s);
            next_available_start_ = vc->get_stop();
            contig_of_next_available_start_ = vc->get_tid();
            AdapterUtils::variant2bcf(bam_herder_, vcf_herder_, vc, output_);
            ret = vcf_format1(vcf_herder_, output_, s);
            CHECK_CONDITION_EXIT(ret != 0, "error: vcf_format1");
        }
    }
}

void BlockCombiner::force_output(kstring_t *s)
{
    emit_current_block(s);  // 将所有位点转换为 bcf1_t
    clear_combiner();
}

void BlockCombiner::clear_combiner()
{
    del_line_offset_ = 0;
    bcf_clear1(output_);
    sample_id_ = INVALID_INT;
    contig_of_next_available_start_ = INVALID_INT;
    next_available_start_ = INVALID_INT;
}

std::map<int32_t, std::pair<int32_t, int32_t>> BlockCombiner::parse_partitions(const std::vector<int32_t> &gq_partitions)
{
    std::map<int32_t, std::pair<int32_t, int32_t>> result;
    int last_threshold{0};

    for (int32_t num : gq_partitions) {
        if (last_threshold + 1 != num) {
            for (int32_t i = last_threshold + 1; i < num; ++i) {
                result.insert({i, {last_threshold, num}});
            }
        }
        CHECK_CONDITION_EXIT(num < 0, "num < 0");
        CHECK_CONDITION_EXIT(num > MAX_GENOTYPE_QUAL + 1, "num > MAX_GENOTYPE_QUAL + 1");
        CHECK_CONDITION_EXIT(num < last_threshold, " num < last_threshold");
        CHECK_CONDITION_EXIT(num == last_threshold, " num == last_threshold");
        result.insert({last_threshold, {last_threshold, num}});
        last_threshold = num;
    }

    if (last_threshold < MAX_GENOTYPE_QUAL) {
        for (int32_t i = last_threshold; i <= MAX_GENOTYPE_QUAL; ++i) {
            result.insert({i, {last_threshold, MAX_GENOTYPE_QUAL + 1}});
        }
    }

    if (MAX_GENOTYPE_QUAL == last_threshold) {
        result.insert({last_threshold, {last_threshold, MAX_GENOTYPE_QUAL + 1}});
    }

    return result;
}

void BlockCombiner::emit_current_block(kstring_t *s)
{
    if (nullptr != current_block_) {
        current_block_->to_variant(vcf_herder_, bam_herder_, sample_id_, floor_blocks_, output_);
        int32_t ret = vcf_format1(vcf_herder_, output_, s);
        CHECK_CONDITION_EXIT(ret != 0, "error: vcf_format1");
        delete current_block_;
        current_block_ = nullptr;
    }
}

bcf1_t *BlockCombiner::add_hom_ref_site(pVariant vc, pGenotype g)
{
    if (next_available_start_ != INVALID_INT) {
        if (vc->get_start() <= next_available_start_ && vc->get_tid() == contig_of_next_available_start_) {
            if (vc->get_stop() <= next_available_start_) {
                return nullptr;
            }
        }
        next_available_start_ = INVALID_INT;
        contig_of_next_available_start_ = INVALID_INT;
    }

    bcf1_t *result;
    if (genotype_can_be_merged_in_current_block(g)) {
        current_block_->add(vc->get_start(), vc->has_end_key() ? vc->end_key() : vc->get_start(), g);
        result = nullptr;
    }
    else {
        result = current_block_ ? current_block_->to_variant(vcf_herder_, bam_herder_, sample_id_, floor_blocks_, output_) : nullptr;
        delete current_block_;
        current_block_ = create_new_block(vc, g);
    }

    return result;
}

bool BlockCombiner::genotype_can_be_merged_in_current_block(pGenotype g)
{
    return current_block_ && current_block_->within_bounds(std::min(g->get_gq(), MAX_GENOTYPE_QUAL)) &&
           current_block_->get_ploidy() == g->get_ploidy() &&
           (current_block_->get_min_pls().size() == g->pl().size() || current_block_->get_min_pls().empty() || !g->has_likelihoods());
}

pHomRefBlock BlockCombiner::create_new_block(pVariant vc, pGenotype g)
{
    int32_t gq = g->has_gq() ? std::min(g->get_gq(), MAX_GENOTYPE_QUAL) : 0;
    CHECK_CONDITION_EXIT(!gq_partitions_.count(gq), "gq {} didn't fit into any partition", gq);
    const auto &tup = gq_partitions_.at(gq);
    return new HomRefBlock{vc, tup.first, tup.second, g->get_ploidy()};
}

}  // namespace rovaca