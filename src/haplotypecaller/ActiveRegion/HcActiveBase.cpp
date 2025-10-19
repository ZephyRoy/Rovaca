#include "HcActiveBase.h"

#include <algorithm>
#include <iostream>

const int HC_REGION_REF_MODEL_DELETION_QUAL = 30;
const uint8_t HC_REGION_MIN_BASE_QUAL = 10;
const int HC_REGION_HQ_BASE_QUALITY_SOFTCLIP_THRESHOLD = 28;
const double HC_REGION_BASES_ACTIVE_HQ_BASES_THRESHOLD = 6.0;
const double HC_REGION_MAX_PROB_DISTANCE = 100;
const double HC_REGION_RC_SNP_P = 0.010000;
const double HC_REGION_RC_REF_P = 10.000000;

unsigned int HcActiveBase::get_high_quality_soft_clips(bam1_t* bam)
{
    unsigned int ret = 0;
    uint8_t* qual = bam_get_qual(bam);
    uint32_t* cigar = bam_get_cigar(bam);
    int qpos = 0;
    for (uint32_t i = 0; i < bam->core.n_cigar; i++) {
        uint32_t cigar_len = bam_cigar_oplen(cigar[i]);
        uint32_t op = bam_cigar_op(cigar[i]);

        if (op == BAM_CSOFT_CLIP) {
            for (uint32_t j = 0; j < cigar_len; j++) {
                if (qual[j + qpos] > HC_REGION_HQ_BASE_QUALITY_SOFTCLIP_THRESHOLD) {
                    ret++;
                }
            }
        }
        if (bam_cigar_type(bam_cigar_op(op)) & 1) qpos += cigar_len;
    }
    return ret;
}

static inline bool is_base_inside_adaptor(bam1_t* bam, hts_pos_t base_pos, int adaptor_boundary)
{
    if (adaptor_boundary == INT_MIN || bam->core.isize > 100) {
        return false;
    }
    return ((bam->core.flag & BAM_FREVERSE) != 0) ? base_pos <= adaptor_boundary : base_pos >= adaptor_boundary;
}

static bool hasWellDefFra(bam1_t* read)
{
    if (read->core.isize == 0 || !(read->core.flag & BAM_FPAIRED) || (read->core.flag & (BAM_FMUNMAP | BAM_FUNMAP)) ||
        ((read->core.flag & BAM_FREVERSE) != 0) == ((read->core.flag & BAM_FMREVERSE) != 0)) {
        return false;
    }
    if ((read->core.flag & BAM_FREVERSE) != 0) {
        return (read->core.pos + bam_cigar2rlen(read->core.n_cigar, bam_get_cigar(read))) > read->core.mpos;
    }
    else {
        return read->core.pos <= read->core.mpos + read->core.isize;
    }
}

int HcActiveBase::get_adaptor_boundary(bam1_t* bam)
{
    if (!hasWellDefFra(bam)) {
        return INT_MIN;
    }
    else if (bam->core.flag & BAM_FREVERSE) {
        return bam->core.mpos + 1 - 1;
    }
    else {
        return bam->core.pos + 1 + abs(bam->core.isize);
    }
}

static inline bool bases_is_snp(uint8_t base, uint8_t ref) { return !(ref == base); }

static inline int bases_is_before_delete_start(uint32_t this_cigar, uint32_t next_cigar, uint32_t bases_pos, uint32_t cigar_len)
{
    return (bases_pos == cigar_len - 1) && this_cigar != BAM_CDEL && next_cigar == BAM_CDEL;
}

static inline bool bases_is_after_delete_end(uint32_t this_cigar, uint32_t prev_cigar, uint32_t bases_pos)
{
    return (bases_pos == 0) && this_cigar != BAM_CDEL && prev_cigar == BAM_CDEL;
}

static inline bool bases_is_before_insert_start(uint32_t next_cigar, uint32_t bases_pos, uint32_t cigar_len)
{
    return (bases_pos == cigar_len - 1) && next_cigar == BAM_CINS;
}

static inline bool bases_is_after_insert_end(uint32_t prev_cigar, uint32_t bases_pos) { return (bases_pos == 0) && prev_cigar == BAM_CINS; }

static inline bool bases_is_next_soft_clip(uint32_t next_cigar, uint32_t prev_cigar, uint32_t bases_pos, uint32_t cigar_len)
{
    return ((bases_pos == cigar_len - 1) && next_cigar == BAM_CSOFT_CLIP) || ((bases_pos == 0) && prev_cigar == BAM_CSOFT_CLIP);
}
static inline bool bases_is_delete(uint32_t this_cigar) { return this_cigar == BAM_CDEL; }

int HcActiveBase::process_bam_to_slot(bam1_t* bam, __attribute__((unused)) int idx, __attribute__((unused)) int min_quality)
{
    uint8_t* qual = bam_get_qual(bam);
    uint8_t* seq = bam_get_seq(bam);
    uint32_t* cigar = bam_get_cigar(bam);
    int adaptorBoundary = get_adaptor_boundary(bam);
    unsigned int hq_soft = get_high_quality_soft_clips(bam);
    hts_pos_t pos = bam->core.pos;
    int qpos = 0;
    for (hts_pos_t iter = m_ringbuffer_hist.m_start; iter < pos; iter++) {
        compute_slot_likelihood(iter - get_acutal_start());
    }
    if (pos > m_ringbuffer_hist.m_start) m_ringbuffer_hist.reset(pos);
    for (uint32_t i = 0; i < bam->core.n_cigar; i++) {
        uint32_t this_op = bam_cigar_op(cigar[i]);
        int cigar_len = bam_cigar_oplen(cigar[i]);
        if (bam_cigar_type(this_op) & 2) {
            uint32_t prev_op = i > 0 ? bam_cigar_op(cigar[i - 1]) : 0xffff;
            uint32_t next_op = i < bam->core.n_cigar - 1 ? bam_cigar_op(cigar[i + 1]) : 0xffff;

            if (bases_is_delete(this_op)) {
                for (int k = 0; k < cigar_len; k++) {
                    int offset = get_actual_offset(bam->core.tid, pos + k);
                    if (offset < 0 || is_base_inside_adaptor(bam, pos + k + 1, adaptorBoundary)) {
                        continue;
                    }
                    m_ringbuffer_hist.increase_hist_storage(pos + k, HCActiveBaseStatus::VariantIDX, HC_REGION_REF_MODEL_DELETION_QUAL);
                }
            }
            else {
                {
                    uint8_t base_qual = qual[qpos];
                    char base_seq = seq_nt16_str[bam_seqi(seq, qpos)];
                    int offset = get_actual_offset(bam->core.tid, pos);
                    if (offset >= 0 && !is_base_inside_adaptor(bam, pos + 1, adaptorBoundary)) {
                        if (prev_op == BAM_CSOFT_CLIP) {
                            if (base_qual > HC_REGION_MIN_BASE_QUAL) m_ringbuffer_hist.add_hq_soft_clips(pos, hq_soft);
                            m_ringbuffer_hist.increase_hist_storage(pos, HCActiveBaseStatus::VariantIDX, base_qual);
                        }
                        else if (bases_is_snp((uint8_t)base_seq, (uint8_t)m_ref[pos]) || prev_op == BAM_CINS || prev_op == BAM_CDEL ||
                                 (cigar_len == 1 && (next_op == BAM_CINS || next_op == BAM_CDEL))) {
                            m_ringbuffer_hist.increase_hist_storage(pos, HCActiveBaseStatus::VariantIDX, base_qual);
                        }

                        else {
                            m_ringbuffer_hist.increase_hist_storage(pos, HCActiveBaseStatus::NONREF, base_qual);
                        }
                    }
                }

                for (int k = 1; k < cigar_len - 1; k++) {
                    int offset = get_actual_offset(bam->core.tid, pos + k);
                    if (offset < 0 || is_base_inside_adaptor(bam, pos + k + 1, adaptorBoundary)) {
                        continue;
                    }
                    uint8_t base_qual = qual[qpos + k];
                    // High Qual data structure
                    char base_seq = seq_nt16_str[bam_seqi(seq, qpos + k)];

                    if (__glibc_unlikely(bases_is_snp((uint8_t)base_seq, (uint8_t)m_ref[pos + k]))) {
                        m_ringbuffer_hist.increase_hist_storage(pos + k, HCActiveBaseStatus::VariantIDX, base_qual);
                    }
                    else {
                        m_ringbuffer_hist.increase_hist_storage(pos + k, HCActiveBaseStatus::NONREF, base_qual);
                    }
                }
                if (cigar_len > 1) {
                    uint8_t base_qual = qual[qpos + cigar_len - 1];
                    char base_seq = seq_nt16_str[bam_seqi(seq, qpos + cigar_len - 1)];
                    int offset = get_actual_offset(bam->core.tid, pos + cigar_len - 1);
                    if (offset >= 0 && !is_base_inside_adaptor(bam, pos + cigar_len, adaptorBoundary)) {
                        if (next_op == BAM_CSOFT_CLIP) {
                            if (base_qual > HC_REGION_MIN_BASE_QUAL) m_ringbuffer_hist.add_hq_soft_clips(pos + cigar_len - 1, hq_soft);
                            m_ringbuffer_hist.increase_hist_storage(pos + cigar_len - 1, HCActiveBaseStatus::VariantIDX, base_qual);
                        }
                        else if (bases_is_snp((uint8_t)base_seq, (uint8_t)m_ref[pos + cigar_len - 1]) || next_op == BAM_CINS ||
                                 next_op == BAM_CDEL) {
                            m_ringbuffer_hist.increase_hist_storage(pos + cigar_len - 1, HCActiveBaseStatus::VariantIDX, base_qual);
                        }

                        else {
                            m_ringbuffer_hist.increase_hist_storage(pos + cigar_len - 1, HCActiveBaseStatus::NONREF, base_qual);
                        }
                    }
                }
            }
            pos += cigar_len;
        }
        if (bam_cigar_type(bam_cigar_op(this_op)) & 1) {
            qpos += cigar_len;
        }
    }
    return 0;
}
void HcActiveBase::compute_slot_likelihood(int offset)
{
    if (get_acutal_start() + offset < m_ringbuffer_hist.m_start) return;
    compute_genotype_PL(offset);
    double active_value = compute_biallelic_non_ref_posterior();
    compute_extension_length(offset, active_value);
}

void HcActiveBase::likelihood_and_count(int qual, HCActiveBaseStatus status, int count)
{
    int offset = status * (m_ploidy + 1) * HCBASEMAXQUAL + qual * (m_ploidy + 1);
    for (int i = 0; i < m_ploidy + 1; i++) {
        m_likehood[i] += count * cache[offset + i];
    }
}
void HcActiveBase::likelihood_and_count(int qual, bool is_alt, int count)
{
    double referenceLikelihood(0.0);
    double nonRefLikelihood(0.0);

    if (is_alt) {
        nonRefLikelihood = qual_to_prob_log10((uint8_t)qual);
        referenceLikelihood = qual_to_error_prob_log10((uint8_t)qual) + s_log10_one_third;
    }
    else {
        referenceLikelihood = qual_to_prob_log10((uint8_t)qual);
        nonRefLikelihood = qual_to_error_prob_log10((uint8_t)qual) + s_log10_one_third;
    }
    m_likehood[0] += count * (referenceLikelihood + m_log10_ploidy);
    m_likehood[m_ploidy] += count * (nonRefLikelihood + m_log10_ploidy);
    for (int i = 1, j = m_ploidy - 1; i < m_ploidy; i++, j--) {
        m_likehood[i] += count * MathUtils::approximate_log10sum_log10(referenceLikelihood + MathUtils::log10(j),
                                                                       nonRefLikelihood + MathUtils::log10(i));
    }
}

void HcActiveBase::compute_extension_length(int offset, double acitve_value)
{
    ActiveResult& result = get_active_actual_result(offset);
    HighCount& hq_result = m_ringbuffer_hist.get_high_count(offset + get_acutal_start());
    //[offset];
    result.acitve_value = acitve_value;
    if (acitve_value == 0) {
        result.state = BaseState::HC_REGION_BASES_ACTIVE_NONE;
        result.extend_length = 0;
    }
    else {
        if (hq_result.mean > HC_REGION_BASES_ACTIVE_HQ_BASES_THRESHOLD) {
            result.state = BaseState::HC_REGION_BASES_ACTIVE_HQ_SOFT_CLIPS;
            result.extend_length = std::min(HC_REGION_MAX_PROB_DISTANCE, hq_result.mean);
        }
        else {
            result.state = BaseState::HC_REGION_BASES_ACTIVE_ONE_SLOT;
            result.extend_length = 0;
        }
    }
    hq_result.clear();
}

void HcActiveBase::compute_genotype_PL(int offset)
{
    std::fill(m_likehood.begin(), m_likehood.end(), 0);
    HCHistStorage& hist_inst = m_ringbuffer_hist.get_hist_inst(get_acutal_start() + offset);
    int read_counts = 0;
    for (int i = 1; i >= 0; i--) {
        uint8_t min_qual = hist_inst.min_idx[i] > HC_REGION_MIN_BASE_QUAL + 1 ? hist_inst.min_idx[i] : HC_REGION_MIN_BASE_QUAL + 1;
        uint8_t max_qual = hist_inst.max_idx[i];
        // bool is_alt = i == HCActiveBaseStatus::VariantIDX;
        HCActiveBaseStatus status = static_cast<HCActiveBaseStatus>(i);
        for (uint8_t qual = min_qual; qual <= max_qual; qual++) {
            const int count = hist_inst.hist_count[i][qual];
            if (count == 0) continue;
            likelihood_and_count(qual, status, count);
            read_counts += count;
            hist_inst.hist_count[i][qual] = 0;
        }
    }
    double denominator = read_counts * m_log10_ploidy;

    for (int i = 0; i < m_ploidy + 1; i++) {
        m_likehood[i] -= denominator;
    }
}

double HcActiveBase::compute_biallelic_non_ref_posterior()
{
    auto max_iter = std::max_element(m_likehood.begin(), m_likehood.end());
    double adjust = *max_iter;
    for (int i = 0; i < m_ploidy + 1; i++) {
        int pls = std::min((int)round(-10 * (m_likehood[i] - adjust)), INT32_MAX);
        m_likehood[i] = pls / -10.0;
    }
    if (max_element_index(m_likehood, 0, m_ploidy + 1) == 0) return 0.0;

    for (int i = 0; i < m_ploidy + 1; i++) {
        m_likehood[i] = m_likehood[i] + log10binomial_coefficient(m_ploidy, i) +
                        log_to_log10(lgamma(i + HC_REGION_RC_SNP_P) + lgamma(m_ploidy - i + HC_REGION_RC_REF_P));
        // std::cerr << "i=" << i << "," << log10binomial_coefficient(m_ploidy, i) << ","
        //           << log_to_log10(lgamma(i + HC_REGION_RC_SNP_P) + lgamma(m_ploidy - i + HC_REGION_RC_REF_P));
    }
    if (max_element_index(m_likehood, 0, m_ploidy + 1) == 0) return 0.0;

    normalize_log10(m_likehood, false, true);
    return 1 - m_likehood[0];
}

void HcActiveBase::increase_hist(int offset, HCActiveBaseStatus status, uint8_t qual)
{
    m_ringbuffer_hist.increase_hist_storage(get_acutal_start() + offset, status, qual);
}