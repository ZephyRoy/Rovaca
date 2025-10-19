/**
 * @file HcActiveBase.h
 * @author yinlonghui
 * @brief
 * HcActiveBase是ActiveBase派生类，它计算逻辑如下：
 * 1.先计算所有基因型的似然值。单倍型只有两种突变和不突变，例如：二倍体基因型只有 （AA，AB AB）
 * 2.所有基因型中算出最大值，做Normalization PL，计算PL ；
 * 3.根据SNP和Indel先验概率，计算Non-Ref的后验概率.
 * Extension逻辑会记录高质量soft-clip的长度，计算均值和方差，用于该位点扩展的长度
 * @version 0.1
 * @date 2023-05-18
 *
 * @copyright Copyright (c) 2023
 *
 */
#ifndef HC_ACTIVE_BASE_H
#define HC_ACTIVE_BASE_H
#include <math.h>

#include "ActiveBase.h"
#include "math_utils.h"
#include "quality_utils.h"

using namespace rovaca;
using namespace math_utils;
using namespace QualityUtils;

enum HCActiveBaseStatus { VariantIDX = 0, NONREF = 1 };
const int HCBASEMAXQUAL = 64;
const int HCBASEMIXQUAL = 0;
struct HCHistStorage
{
    HCHistStorage(std::pmr::memory_resource *pool)
    {
        int *mem = (int *)pool->allocate(sizeof(int) * HCBASEMAXQUAL * 2);
        memset(mem, 0, sizeof(int) * HCBASEMAXQUAL * 2);
        hist_count[0] = mem;
        hist_count[1] = mem + HCBASEMAXQUAL;
        min_idx[0] = min_idx[1] = HCBASEMAXQUAL;
        max_idx[0] = max_idx[1] = HCBASEMIXQUAL;
    }
    void insert(HCActiveBaseStatus status, uint8_t qual)
    {
        hist_count[status][qual]++;
        if (__glibc_unlikely(qual < min_idx[status])) min_idx[status] = qual;
        if (__glibc_unlikely(qual > max_idx[status])) max_idx[status] = qual;
    }
    int *hist_count[2];
    uint8_t max_idx[2];
    uint8_t min_idx[2];
};

struct HighCount
{
    /*不确定要不要写多线程算这个，先来构造*/
    HighCount(double mean_, double s_, long obs_count_)
        : mean(mean_)
        , s(s_)
        , obs_count(obs_count_){};
    void add_hq_soft_clips(double obs)
    {
        obs_count++;
        double old_mean = mean;
        mean += (obs - mean) / obs_count;
        s += (obs - old_mean) * (obs - mean);
    }
    void clear()
    {
        mean = s = 0.0;
        obs_count = 0;
    }
    double mean;
    double s;
    long obs_count;
};

struct ringBufferHistStorage
{
    static constexpr hts_pos_t buffer_size = 1024;
    static constexpr hts_pos_t mask = buffer_size - 1;
    ringBufferHistStorage(std::pmr::memory_resource *pool, hts_pos_t actual_start, hts_pos_t actual_end)
        : m_hist_count(pool)
        , m_hq_count(pool)
    {
        hts_pos_t actual_size = std::min(actual_end - actual_start, buffer_size);
        m_hist_count.reserve(actual_size);
        m_hq_count.reserve(actual_size);
        for (hts_pos_t iter = 0; iter < actual_size; iter++) {
            m_hist_count.emplace_back(pool);
            m_hq_count.emplace_back(0.0, 0.0, 0);
        }
        m_start = actual_start;
        start_index = 0;
    }

    void increase_hist_storage(hts_pos_t pos, HCActiveBaseStatus status, uint8_t qual)
    {
        int offset = (pos - m_start + start_index) & mask;
        m_hist_count[offset].insert(status, qual);
    }
    void add_hq_soft_clips(hts_pos_t pos, double obs)
    {
        int offset = (pos - m_start + start_index) & mask;
        m_hq_count[offset].add_hq_soft_clips(obs);
    }

    HCHistStorage &get_hist_inst(hts_pos_t pos)
    {
        int offset = (pos - m_start + start_index) & mask;
        return m_hist_count[offset];
    }
    HighCount &get_high_count(hts_pos_t pos)
    {
        int offset = (pos - m_start + start_index) & mask;
        return m_hq_count[offset];
    }

    /*
     *  Note:在Reset之前要确保数据清空
     */
    void reset(hts_pos_t pos)
    {
        int offset = (pos - m_start + start_index) & mask;
        m_start = pos;
        start_index = offset;
    }

    hts_pos_t m_start;
    int start_index;
    std::pmr::vector<HCHistStorage> m_hist_count;
    // 计算High Quality的均值和方差
    std::pmr::vector<HighCount> m_hq_count;
};

class HcActiveBase : public ActiveBase
{
public:
    HcActiveBase(std::pmr::memory_resource *pool, int ploidy, int tid, hts_pos_t start, hts_pos_t end, hts_pos_t actual_start,
                 hts_pos_t actual_end, boost::dynamic_bitset<> *bit, int work_id, char *ref, int ref_len)
        : ActiveBase(pool, tid, start, end, actual_start, actual_end, bit, work_id)
        , m_ringbuffer_hist(pool, actual_start, actual_end)
        , m_likehood(pool)
        , cache(pool)
        // , m_hist_count(pool)
        // , m_hq_count(pool)
        , m_ref(ref)
        , m_ref_len(ref_len)
        , m_ploidy(ploidy)
    {
        // m_hist_count.reserve(actual_end - actual_start);
        m_likehood.resize(ploidy + 1);
        m_log10_ploidy = log10(m_ploidy);
        for (int i = 0; i < 2; i++) {
            bool is_alt = i == HCActiveBaseStatus::VariantIDX;
            for (int qual = 0; qual < HCBASEMAXQUAL; qual++) {
                if (is_alt) {
                    double referenceLikelihood = qual_to_error_prob_log10((uint8_t)qual) + s_log10_one_third;
                    double nonRefLikelihood = qual_to_prob_log10((uint8_t)qual);
                    cache.push_back(referenceLikelihood + m_log10_ploidy);

                    for (int iter = 1; iter < m_ploidy; iter++) {
                        cache.push_back(MathUtils::approximate_log10sum_log10(referenceLikelihood + MathUtils::log10(m_ploidy - iter),
                                                                              nonRefLikelihood + log10(iter)));
                    }
                    cache.push_back(nonRefLikelihood + m_log10_ploidy);
                }
                else {
                    double referenceLikelihood = qual_to_prob_log10((uint8_t)qual);
                    double nonRefLikelihood = qual_to_error_prob_log10((uint8_t)qual) + s_log10_one_third;
                    cache.push_back(referenceLikelihood + m_log10_ploidy);

                    for (int iter = 1; iter < m_ploidy; iter++) {
                        cache.push_back(MathUtils::approximate_log10sum_log10(referenceLikelihood + MathUtils::log10(m_ploidy - iter),
                                                                              nonRefLikelihood + log10(iter)));
                    }
                    cache.push_back(nonRefLikelihood + m_log10_ploidy);
                }
            }
        }
    }
    /**
     * @brief Construct a new Hc Active Base object ，初始化ActiveStroage 1维度
     *
     * @param tid 处理数据集区间的染色体坐标
     * @param start 处理数据集的开始坐标
     * @param end 处理数据集的结束坐标
     * @param bit 处理区间bitset，插入先需要判断是否在bed
     * @param work_id ActiveBase的标号，用于多线程处理时候标识ActiveBase的顺序。
     */
    HcActiveBase(std::pmr::memory_resource *pool, int ploidy, int tid, hts_pos_t start, hts_pos_t end, boost::dynamic_bitset<> *bit,
                 int work_id, char *ref, int ref_len)
        : HcActiveBase(pool, ploidy, tid, start, end, start, end, bit, work_id, ref, ref_len)
    {}

    ~HcActiveBase() {}

    /**
     * @brief 处理reads（bam1_t)：
     * 通过解析reads的cigar结构，判断该碱基是突变位点还是非突变位点。
     * 把该碱基的状态和质量存储到ActiveStorage，存储改位点的突变和非突变的直方图（ActiveStroage）。
     *
     *
     * @param bam  待处理的bam文件
     * @param idx  处理bam来自于第几个bam
     * @param min_quality 过滤碱基质量小于min_quality的碱基
     * @return int reads成功处理返回0，否则返回1 （reads在区间外）
     */
    virtual int process_bam_to_slot(bam1_t *bam, int idx, int min_quality) override;

protected:
    /**
     * @brief 计算偏移位点的基因型的ActiveResult
     *
     * @param offset 保存ActiveStroage的偏移量
     */
    virtual void compute_slot_likelihood(int offset) override;

public:
    /*为了单元测试，这个些接口只能在单元测试用*/
    static unsigned int get_high_quality_soft_clips(bam1_t *bam);

    static int get_adaptor_boundary(bam1_t *bam);

    void increase_hist(int offset, HCActiveBaseStatus status, uint8_t qual);

    std::pmr::vector<double> get_current_likelihood() { return m_likehood; }

    void insert_high_count(int offset, int hq_len)
    {
        m_ringbuffer_hist.add_hq_soft_clips(offset + get_acutal_start(), hq_len);
        // m_hq_count[offset].add_hq_soft_clips(hq_len);
    }

public:
    /**
     * @brief 使用记录高质量soft-clip的长度，计算均值和方差，计算该位点扩展的长度
     *
     * @param offset 偏移量
     */
    void compute_extension_length(int offset, double acitve_value);
    /**
     * @brief 计算所有基因型（突变和非突变的）的Normalize PL值
     *
     * @param offset 偏移量
     */

    void compute_genotype_PL(int offset);

    /**
     * @brief 根据SNP和Indel先验概率，计算Non-Ref的后验概率
     *
     * @param offset  偏移量
     * @return double
     */
    double compute_biallelic_non_ref_posterior();

    void likelihood_and_count(int qual, bool is_alt, int count);

    void likelihood_and_count(int qual, HCActiveBaseStatus status, int count);

private:
    ringBufferHistStorage m_ringbuffer_hist;
    std::pmr::vector<double> m_likehood;
    std::pmr::vector<double> cache;

    // 用于HC定义直方图结果，它只用存储 NON_REF REF的碱基直方图
    // std::pmr::vector<HCHistStorage> m_hist_count;
    // 计算High Quality的均值和方差
    // std::pmr::vector<HighCount> m_hq_count;
    char *m_ref;
    int m_ref_len;
    int m_ploidy;
    double m_log10_ploidy;
};
#endif
