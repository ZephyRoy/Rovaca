#ifndef DOWNSAMPLER_H
#define DOWNSAMPLER_H

#include <iostream>
#include <memory>
#include <vector>

#include "htslib/sam.h"

class Downsampler
{
private:
    int discarded_items_count;

public:
    Downsampler()
        : discarded_items_count(0){};
    virtual ~Downsampler(){};

    /**
     * @brief downsample主接口，向downsample池子中提交一条read。
     * @param item 一条read.
     * @return bam1_t* 正常为NULL，若不为空，则为需要抛弃的read，需做内存处理。
     */
    virtual bam1_t* submit(bam1_t* item) = 0;

    /**
     * @brief 一个downsample周期是否结束
     * @return true
     * @return false
     */
    virtual bool has_finalized_items() = 0;

    /**
     * @brief downsample结果获取接口。
     * @param cache
     * @return int
     */
    virtual int consume_finalized_items(std::list<bam1_t*>& cache) = 0;

    /**
     * @brief 数据流结束信号。
     */
    virtual void input_end_signal() = 0;
    int get_number_of_discarded_items() { return discarded_items_count; }
    void increment_number_of_discarded_items() { discarded_items_count++; }
    void reset_stats() { discarded_items_count = 0; }
};

#endif  // DOWNSAMPLER_H