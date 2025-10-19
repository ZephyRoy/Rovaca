#ifndef REGION_INTERFACE_H
#define REGION_INTERFACE_H
#include <memory>
#include <memory_resource>
#include <vector>

#include "assemble/assemble_interface.h"
#include "htslib/vcf.h"
#include "utils/bam_data_pool.hpp"

struct RegionResource
{
    RegionResource(){};
    RegionResource(int resource_id, int pool_size, int buffer_size, size_t consumer_result_size)
        : id_(resource_id)
        , buffer_size_(buffer_size)
        , result_buffer_size_(consumer_result_size)
    {
        pool_ = new rovaca::BamDataPool(pool_size);
        buffer_ = new uint8_t[buffer_size]{};
        assemble_result_buffer_ = new uint8_t[result_buffer_size_]{};
    }

    void finalize() { pool_->finalize(); }

    ~RegionResource()
    {
        delete pool_;
        delete[] buffer_;
        delete[] assemble_result_buffer_;
    }
    rovaca::BamDataPool *pool_;
    uint8_t *assemble_result_buffer_;
    uint8_t *buffer_;  // 用于内存池，消费者多个模块内部去构建memory_resource对象，memory_resource对象周期只在对应模块
                       // 如果一个模块有单独的空间，需要常驻传给后面的模块，可以考虑多加几个buffer
    int id_;
    int buffer_size_;
    size_t result_buffer_size_;
};

struct RegionResult
{
    RegionResult(std::vector<bam1_t *> reads, p_hc_region_active_storage new_region, int id)
        : padding_reads(reads)
        , region(new_region)
        , region_id(id)
    {}
    // reads周期贯穿整个消费者周期，存的是指针,指针对象由BamDataPool，暂时无性能瓶颈
    std::vector<bam1_t *> padding_reads;
    p_hc_region_active_storage region;
    // 记录 RegionResult的id
    int region_id;
};

// Queue的模板必须是shared_ptr<RegionSource>
struct RegionSource
{
    RegionSource() {}

    RegionSource(std::shared_ptr<char> ref, int tid, int ref_len, std::vector<RegionResult> regions, int source_id, bool finish,
                 std::shared_ptr<RegionResource> resource)
        : ref_(ref)
        , tid_(tid)
        , ref_len_(ref_len)
        , regions_(regions)
        , source_id_(source_id)
        , finish_(finish)
        , resource_(resource)
    {}

    ~RegionSource() = default;

    void release()
    {
        for (const auto &result : regions_) {
            for (auto read : result.padding_reads) {
                bam_destroy1(read);
            }
            free(result.region);
        }

        db_offset_ = 0;
        db_data_ = nullptr;
    }

    std::shared_ptr<char> ref_;  // ref碱基(整条染色体)  不需要释放
    int tid_;                    //
    // ref长度
    int ref_len_;
    /*
    VCF模式下：RegionResult的padding reads的空间占resource的BamDataPool总空间的一半，另外一半留给Genotype
    GVCF模式下：RegionResult的padding reads的空间占resource的BamDataPool总空间的1/4，另外3/4留给Genotype
    */
    std::vector<RegionResult> regions_;
    int source_id_;
    bool finish_;
    // 消费者通过队列，把resource传递给ActiveRegion的Reduce主线程
    std::shared_ptr<RegionResource> resource_;

    // 若存在db，db_offset_指示当前tid_对应的数据db_data_
    int32_t db_offset_{0};
    std::shared_ptr<std::vector<bcf1_t *>> db_data_{nullptr};
};
#endif