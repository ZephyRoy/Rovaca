/**
 * @file ActiveBase.h
 * @author yinlonghui
 * @brief
 * ActiveBase是ActiveRegion主要耗时模块，设计思想先计算好一堆reads的区间，然后每个区间放到线程池里面并发计算每个位点的后验概率以及扩展长度。
 *         ActiveRegionEngine处理ActiveBase计算出来每个位点的结果，再将似然值向左右扩展，根据高斯核的值累加其AcitveScore。最后根据GATK极小值逻辑输出激活区或非激活区间。
 * @version 0.1
 * @date 2023-04-26
 *
 * @copyright Copyright (c) 2023
 *
 */
#ifndef ACTIVBASE_H
#define ACTIVBASE_H
#include <boost/dynamic_bitset.hpp>
#include <iostream>
#include <memory_resource>

#include "htslib/sam.h"
// #include <vector>

enum BaseState { HC_REGION_BASES_ACTIVE_NONE = 0, HC_REGION_BASES_ACTIVE_HQ_SOFT_CLIPS = 1, HC_REGION_BASES_ACTIVE_ONE_SLOT = 2 };

struct ActiveResult
{
    double acitve_value;
    BaseState state;
    int extend_length;
    ActiveResult(double acitve_value_, BaseState state_, double extend_length_)
        : acitve_value(acitve_value_)
        , state(state_)
        , extend_length(int(extend_length_))
    {}

    ActiveResult(double acitve_value_, BaseState state_, int extend_length_)
        : acitve_value(acitve_value_)
        , state(state_)
        , extend_length(extend_length_)
    {}
};
static ActiveResult out_of_result(0, BaseState::HC_REGION_BASES_ACTIVE_NONE, 0);
class ActiveBase
{
public:
    /**
     * @brief Construct a new Active Base object
     *
     * @param tid 处理数据集区间的染色体坐标
     * @param start 处理数据集的开始坐标
     * @param end 处理数据集的结束坐标  [需要把考虑没有覆盖区间的问题]
     * @param actual_start 实际处理有效开始坐标
     * @param actual_end 实际处理有效结束坐标
     * @param bit 处理区间bitset，插入先需要判断是否在bed
     * @param work_id ActiveBase的标号，用于多线程处理时候标识ActiveBase的顺序。
     */
    ActiveBase(std::pmr::memory_resource *pool, int tid, hts_pos_t start, hts_pos_t end, hts_pos_t actual_start, hts_pos_t actual_end,
               boost::dynamic_bitset<> *bit, int work_id)
        : m_work_id(work_id)
        , m_tid(tid)
        , m_start(start)
        , m_end(end)
        , m_actual_start(actual_start)
        , m_actual_end(actual_end)
        , m_active(pool)
        , m_bit(bit)
    {
        for (hts_pos_t iter = actual_start; iter < actual_end; iter++) {
            m_active.emplace_back(0.0, BaseState::HC_REGION_BASES_ACTIVE_NONE, 1);
        }
    }

    /**
     * @brief Construct a new Active Base object
     *
     * @param tid 处理数据集区间的染色体坐标
     * @param start 处理数据集的开始坐标
     * @param end 处理数据集的结束坐标  [需要把考虑没有覆盖区间的问题]
     * @param bit 处理区间bitset，插入先需要判断是否在bed
     * @param work_id ActiveBase的标号，用于多线程处理时候标识ActiveBase的顺序。
     */
    ActiveBase(std::pmr::memory_resource *pool, int tid, hts_pos_t start, hts_pos_t end, boost::dynamic_bitset<> *bit, int work_id)
        : ActiveBase(pool, tid, start, end, start, end, bit, work_id)
    {}

    virtual ~ActiveBase(){};

    /**
     * @brief 处理reads：
     * 通过解析reads的cigar结构，判断该碱基是突变位点还是非突变位点。（HC/Mutect逻辑不一样子类判断）
     * 把该碱基的状态和质量存储到ActiveStorage（由子类实现，HC/Mutect的直方图结构不一样），大体还是用直方图
     *
     * @param bam  待处理的bam文件
     * @param idx  处理bam来自于第几个bam
     * @param min_quality 过滤碱基质量小于min_quality的碱基
     * @return int reads成功处理返回0，否则返回1 （reads在区间外）
     */
    virtual int process_bam_to_slot(bam1_t *bam, int idx, int min_quality) = 0;

    /**
     * @brief 计算区间内所有位点的似然值，调用compute_slot_likelihood后
     *
     * @return int 成功 0  失败 1（无reads覆盖此区间)
     */
    bool compulte_all_likelihood();
    /**
     * @brief 偏移对应的位点是target预取，wgs应该返回1
     *
     * @param offset 偏移量 [非actual offset]
     * @return true
     * @return false
     */

    bool test_in_target(int offset) { return !m_bit || m_bit->test(offset); }

    bool test_in_actual_target(int offset) { return !m_bit || m_bit->test(offset); }

    int get_extend_length(int offset)
    {
        int actual_offset = offset_transform_actual(offset);
        auto result = get_active_result(actual_offset);
        if (result.state != BaseState::HC_REGION_BASES_ACTIVE_HQ_SOFT_CLIPS)
            return 1;
        else
            return result.extend_length * 2 + 1;
    }

    int get_extend_actual_length(int offset)
    {
        auto result = get_active_actual_result(offset);
        if (result.state != BaseState::HC_REGION_BASES_ACTIVE_HQ_SOFT_CLIPS)
            return 1;
        else
            return result.extend_length * 2 + 1;
    }
    /**
     * @brief Get the active result object
     *
     * @param offset  偏移量
     * @return ActiveResult  偏移量的ActiveResult
     */
    ActiveResult &get_active_result(int offset)
    {
        int actual_offset = offset_transform_actual(offset);
        if (actual_offset == -1) return out_of_result;
        return m_active[actual_offset];
    }

    ActiveResult &get_active_actual_result(int actual_offset) { return m_active[actual_offset]; }

    /**
     * @brief Get the work id
     *
     * @return int  小于0是存在错误
     */
    int get_work_id() { return m_work_id; }

    int get_offset(int tid, hts_pos_t pos);

    int get_actual_offset(int tid, hts_pos_t pos);

    int offset_transform_actual(int offset);

    int get_tid() { return m_tid; }

    hts_pos_t get_acutal_start() { return m_actual_start; }

    hts_pos_t get_acutal_end() { return m_actual_end; }

    hts_pos_t get_actual_size() { return m_actual_end - m_actual_start; }

    hts_pos_t get_start() { return m_start; }

    hts_pos_t get_stop() { return m_end; }

    hts_pos_t get_size() { return m_end - m_start; }

protected:
    /**
     * @brief 计算偏移位点的基因型的ActiveResult
     *
     * @param offset 保存ActiveStroage的偏移量
     * @return ActiveResult
     */
    virtual void compute_slot_likelihood(int offset) = 0;

private:
    int m_work_id;
    int m_tid;
    hts_pos_t m_start;
    hts_pos_t m_end;
    hts_pos_t m_actual_start;
    hts_pos_t m_actual_end;
    std::pmr::vector<ActiveResult> m_active;
    // std::pmr::vector<double> m_extension_value;
    boost::dynamic_bitset<> *m_bit;
    char *ref;
};
#endif