#include <algorithm>
#include <boost/asio.hpp>
#include <boost/asio/thread_pool.hpp>
#include <boost/dynamic_bitset.hpp>
#include <boost/lockfree/queue.hpp>
#include <boost/thread/mutex.hpp>
#include <fstream>
#include <iostream>
#include <list>
#include <mutex>
#include <queue>

#include "ActiveBase.h"
#include "ActiveRegionEngine.h"
#include "BlockingQueue.h"
#include "RegionManager.h"
#include "bam_data_pool.hpp"
#include "bed_loader.h"
#include "bqsr_read_transformer.h"
#include "fasta_loader.h"
#include "htslib/sam.h"
#include "reads_stream.h"
#include "reference_manager.h"
#include "region_interface.h"

struct ActiveBaseResource
{
    ActiveBaseResource()
    {
        m_target = 0;
        pool = nullptr;
        target = nullptr;
        buffer = nullptr;
    };
    ActiveBaseResource(int target_num, int buffer_sz)
    {
        m_buffer_sz = buffer_sz;
        buffer = new uint8_t[buffer_sz]{};
        pool = nullptr;
        // pool = new std::pmr::monotonic_buffer_resource(buffer, buffer_sz, std::pmr::new_delete_resource());
        target = new boost::dynamic_bitset<>(target_num, 0x0);
        target->set();
        m_target = target_num;
    }

    void new_mempool()
    {
        if (pool) delete pool;
        pool = new std::pmr::monotonic_buffer_resource(buffer, m_buffer_sz, std::pmr::new_delete_resource());
    }

    ~ActiveBaseResource()
    {
        delete[] buffer;
        if (pool) delete pool;
        delete target;
    }

    std::pmr::memory_resource *pool;
    int m_buffer_sz;
    boost::dynamic_bitset<> *target;
    // capablity of target
    int m_target;
    // new_delete_resource
    uint8_t *buffer;
};

struct BaseSource
{
    std::shared_ptr<ActiveBaseResource> base_source;
    ActiveBase *base;
    bool finish;
};

struct BamSource
{
    static constexpr int buffer_sz = 64 * 1024 * 8;
    BamSource()
    {
        buffer = new uint8_t[buffer_sz]{};
        pool = nullptr;
    }

    void new_mempool()
    {
        if (pool) delete pool;
        pool = new std::pmr::monotonic_buffer_resource(buffer, buffer_sz, std::pmr::new_delete_resource());
    }

    ~BamSource()
    {
        delete[] buffer;
        if (pool) delete pool;
    }
    std::pmr::memory_resource *pool;
    uint8_t *buffer;
};

class ActiveMainThreadDispatchTasks
{
    static constexpr int64_t skip_bp_length = 300;
    static constexpr int batch_reads = 50000;

public:
    ActiveMainThreadDispatchTasks(ReadStream *stream, BedLoader *loader, ReferenceManager *m_fasta_loader, contig_info_t *fasta_info,
                                  ActiveRegionBamBlockListSource *block_resource, boost::asio::thread_pool *thread_pool,
                                  BlockingQueue<std::shared_ptr<ActiveBaseResource>> *base_resource, BlockingQueue<BaseSource> *queue,
                                  BlockingQueue<std::shared_ptr<BamSource>> *bam_queue, BQSRReadTransformer *apply_bqsr)
        : m_stream(stream)
        , m_bed_loader(loader)
        , m_fasta_loader(m_fasta_loader)
        , m_fasta_info(fasta_info)
        , m_block_resource(block_resource)
        , m_run_pool(thread_pool)
        , m_base_resource(base_resource)
        , m_reduce_queue(queue)
        , m_bam_resource(bam_queue)
        , m_apply_bqsr(apply_bqsr)
    {
        // m_bam_resource = new BlockingQueue<std::shared_ptr<BamSource>>(32);
        // for (int i = 0; i < 32; i++) {
        //     std::shared_ptr<BamSource> bam_source = std::make_shared<BamSource>();
        //     m_bam_resource->push(bam_source);
        // }
    }

    ~ActiveMainThreadDispatchTasks() {}
    /**
     * @brief ActiveRegion的分发主线程：
     * 先看m_idle_pool_queue是否有可用内存池，没有则等待，有则开始如下面的操作：
     * 1.从InputReadsStream接口读入Reads,按照策略把Reads拆封成一个区间的一堆Reads（调用split_reads方法）。
     * 2.如果有BED文件，根据区间算出该区间是否属于target。（调用get_target_biset）
     * 3.将ActiveBase计算逻辑丢入线程池；
     * Note：计算区间的时间需要考虑那些没有覆盖的reads，要求ActiveBase的区间都是连续的。
     */

    void run();

    /**
     * @brief
     *
     * @param tid
     * @param start
     * @param end
     * @param actual_start
     * @param actual_end
     * @param cover_reads
     * @param repeat
     * @return true
     * @return false
     */
    bool split_reads(int &tid, hts_pos_t &start, hts_pos_t &end, hts_pos_t &actual_start, hts_pos_t &actual_end,
                     std::pmr::vector<bam1_t *> &cover_reads);

    /**
     * @brief Get the target biset object
     *
     * @param tid  染色体id
     * @param start  区间开始
     * @param end  区间结束
     * @param target  bitset WGS:全是1，WES:里面的值为1就表示在bed区间里，0表示不在
     * @param m  bitset的大小
     * @return true
     * @return false
     */
    bool get_target_biset(int tid, hts_pos_t start, hts_pos_t end, boost::dynamic_bitset<> *target, int &m);

private:
    bool fetch_block_reads();

    void split_block_resource(BamBlockList &block, int &tid, hts_pos_t &end, hts_pos_t &actual_start, hts_pos_t &actual_end,
                              std::pmr::vector<bam1_t *> &cover_reads, bool cover_upstream);

private:
    ReadStream *m_stream;
    BedLoader *m_bed_loader;
    ReferenceManager *m_fasta_loader;
    contig_info_t *m_fasta_info;
    ActiveRegionBamBlockListSource *m_block_resource;
    boost::asio::thread_pool *m_run_pool;
    BlockingQueue<std::shared_ptr<ActiveBaseResource>> *m_base_resource;
    BlockingQueue<BaseSource> *m_reduce_queue;
    BlockingQueue<std::shared_ptr<BamSource>> *m_bam_resource;
    BQSRReadTransformer *m_apply_bqsr;
    // int iter_current_
    p_bed_intervals current_interval{nullptr};
    int interval_tid{-1};
    int interval_offset{0};
};

class ActiveMainThreadReduce

{
    static constexpr int read_padded_span = 100;

public:
    ActiveMainThreadReduce(ReferenceManager *fasta_loader, contig_info_t *fasta_info, BedLoader *bed_loader,
                           BlockingQueue<std::shared_ptr<RegionResource>> *region_resource,
                           BlockingQueue<std::shared_ptr<ActiveBaseResource>> *base_resource,
                           ActiveRegionBamBlockListSource *block_resource, BlockingQueue<BaseSource> *dispatch_queue,
                           BlockingQueue<std::shared_ptr<RegionSource>> *region_queue, bool force_non_active, DbsnpManager *db)
        : m_fasta_loader(fasta_loader)
        , m_region_resource(region_resource)
        , m_base_resource(base_resource)
        , m_block_resource(block_resource)
        , m_dispatch_queue(dispatch_queue)
        , m_region_queue(region_queue)
        , force_non_active_(force_non_active)
    {
        engine = new ActiveRegionEngine(50, 300, bed_loader, fasta_info, force_non_active);
        region_manager = new RegionManager(engine, fasta_loader, region_queue, block_resource, region_resource, force_non_active, db);
    }

    /**
     * @brief  ActiveRegion归并任务的主线程：
     * 1.调用ActiveRegionEngine，合并多个ActiveBase结果，输出Region区间
     * 2.往结果队列里push RegionActiveResult。
     *
     */
    void run();

    ~ActiveMainThreadReduce()
    {
        delete engine;
        delete region_manager;
    }

    bool get_region_reads_block(rovaca::BamDataPool *pool, p_hc_region_active_storage region, std::pmr::vector<bam1_t *> &reads);

    bool remove_non_overlapping_reads_region_block(p_hc_region_active_storage region);

private:
    void region_handle(int &region_id);

    void region_finish();

private:
    ActiveRegionEngine *engine;

    // Non-use
    std::map<int, std::shared_ptr<ActiveBaseResource>> m_cache;
    ReferenceManager *m_fasta_loader;
    BlockingQueue<std::shared_ptr<RegionResource>> *m_region_resource;
    BlockingQueue<std::shared_ptr<ActiveBaseResource>> *m_base_resource;
    ActiveRegionBamBlockListSource *m_block_resource;
    BlockingQueue<BaseSource> *m_dispatch_queue;
    BlockingQueue<std::shared_ptr<RegionSource>> *m_region_queue;
    RegionManager *region_manager;
    int last_tid{0};
    bool force_non_active_;
    std::list<BamBlockList>::iterator block_itertor;
    std::list<RegionReadRecord>::iterator bam_itertor;
};

class FakeConsumer
{
public:
    FakeConsumer(BlockingQueue<std::shared_ptr<RegionResource>> *region_source, BlockingQueue<std::shared_ptr<RegionSource>> *region_queue,
                 std::string outfilename, bam_hdr_t *h, ReferenceManager *reference_manager)
        : m_region_source(region_source)
        , m_region_queue(region_queue)
        , out(outfilename)
        , hdr(h)
        , m_reference_manager(reference_manager)
    {}

    void run();

    ~FakeConsumer() {}

private:
    BlockingQueue<std::shared_ptr<RegionResource>> *m_region_source;
    BlockingQueue<std::shared_ptr<RegionSource>> *m_region_queue;
    std::ofstream out;
    bam_hdr_t *hdr;
    ReferenceManager *m_reference_manager;
};
