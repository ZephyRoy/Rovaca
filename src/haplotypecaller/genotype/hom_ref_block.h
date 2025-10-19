#ifndef ROVACA_HC_HOM_REF_BLOCK_H_
#define ROVACA_HC_HOM_REF_BLOCK_H_
#include <unordered_set>

#include "forward.h"
#include "htslib/vcf.h"
#include "interface/interface_locatable.hpp"

namespace rovaca
{

class HomRefBlock : public InterfaceLocatable
{
    /* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
private:
    int32_t tid_;
    int64_t start_;
    int64_t end_;

    pAllele ref_allele_;
    pAllele alt_allele_;

    int32_t lower_gq_bound_;
    int32_t upper_gq_bound_;

    int32_t ploidy_;
    int32_t min_gq_;
    int32_t min_dp_;
    std::multiset<int32_t> dps_;
    std::vector<int32_t> min_pls_;

    /* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
public:
    ~HomRefBlock() override = default;
    int32_t get_tid() const override { return tid_; }
    int64_t get_start() const override { return start_; }
    int64_t get_stop() const override { return end_; }
    void set_tid(int32_t tid) override { tid_ = tid; }
    void set_start(int64_t start) override { start_ = start; }
    void set_stop(int64_t stop) override { end_ = stop; }

    pAllele ref() const { return ref_allele_; }
    pAllele alt() const { return alt_allele_; }
    int64_t size() const { return end_ - start_ + 1; }

    HomRefBlock(pVariant start_vc, int32_t lower_gq_bound, int32_t upper_gq_bound, int32_t default_ploidy);

    bool is_contiguous(pVariant vc) const;
    int32_t get_ploidy() const { return ploidy_; }
    const std::vector<int32_t>& get_min_pls() const { return min_pls_; }
    bool within_bounds(int32_t gq) const { return gq >= lower_gq_bound_ && gq < upper_gq_bound_; }

    void add(int64_t pos, int64_t end, pGenotype g);
    bcf1_t* to_variant(bcf_hdr_t* vcf_herder, bam_hdr_t* bam_herder, int32_t sample_id, bool floor_blocks, bcf1_t* b);

    /* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
private:
    int32_t median_dp() const;

    int32_t min_gq() const { return min_gq_; }
    int32_t gq_lowerbound() const { return lower_gq_bound_; }
};

}  // namespace rovaca

#endif  // ROVACA_HC_HOM_REF_BLOCK_H_
