
#ifndef REGION_MANAGER_H
#define REGION_MANAGER_H
#include <iostream>
#include <list>

#include "ActiveRegionEngine.h"
#include "BlockingQueue.h"
#include "dbsnp_manager.h"
#include "rovaca_logger.h"
#include "reference_manager.h"
#include "region_interface.h"

struct RegionReadRecord
{
    RegionReadRecord(hts_pos_t pos, bam1_t *r)
        : ref_end_pos(pos)
        , bam(r)
    {}
    hts_pos_t ref_end_pos;
    bam1_t *bam;
};

const int raw_buffer_size = 1024;
struct BamBlockList
{
    BamBlockList() {}

    BamBlockList(RegionReadRecord &record, std::list<RegionReadRecord> &idle, bool out_interval)
    {
        start = record.bam->core.pos;
        end = record.bam->core.pos;
        tid = record.bam->core.tid;
        max_end = record.ref_end_pos;
        is_out_interval = false;
        is_out_chrom = false;
        is_prev_out_interval = out_interval;
        bams.splice(bams.end(), idle, idle.begin());
    }

    void insert(RegionReadRecord &record, std::list<RegionReadRecord> &idle)
    {
        end = record.bam->core.pos;
        max_end = std::max(max_end, record.ref_end_pos);
        bams.splice(bams.end(), idle, idle.begin());
    }
    int tid;
    hts_pos_t start;    // first reads record's coordinate
    hts_pos_t end;      // last reads record's coordinate
    hts_pos_t max_end;  // max reads endPos
    bool is_out_interval;
    bool is_out_chrom;
    bool is_prev_out_interval;
    std::list<RegionReadRecord> bams;
};

struct ActiveRegionBamBlockListSource
{
    static constexpr int skip_bps = 1000;

    ActiveRegionBamBlockListSource(int cap_size)
    {
        for (int i = 0; i < cap_size; i++) {
            bam1_t *bam = bam_init1();
            m_reads_idle.emplace_back(0, bam);
        }
        m_finalize = true;
    }

    ~ActiveRegionBamBlockListSource()
    {
        int num_reads = 0;
        uint64_t cap_size = 0;
        for (auto block : m_reads_buffer) {
            for (auto reads : block.bams) {
                cap_size += sizeof(bam1_t);
                cap_size += reads.bam->m_data;
                bam_destroy1(reads.bam);
                num_reads++;
            }
        }
        for (auto reads : m_reads_idle) {
            cap_size += sizeof(bam1_t);
            cap_size += reads.bam->m_data;
            bam_destroy1(reads.bam);
            num_reads++;
        }
        RovacaLogger::info("read number = {}, cap size = {}", num_reads, cap_size);
    }

    RegionReadRecord &get_back_bam()
    {
        if (m_reads_idle.empty()) {
            bam1_t *bam = bam_init1();
            m_reads_idle.emplace_back(0, bam);
        }
        return m_reads_idle.front();
    }

    bool is_next_adjacent(bam1_t *bam)
    {
        return m_reads_buffer.empty() || m_reads_buffer.back().bams.back().bam->core.pos != bam->core.pos;
    }

    bool insert(bam1_t *bam)
    {
        RegionReadRecord &buffer = get_back_bam();
        buffer.ref_end_pos = bam_endpos(bam);
        // swap
        // swap.data = buffer.bam->data;
        // swap.m_data = buffer.bam->m_data;
        // *buffer.bam = *bam;
        // bam->m_data = swap.m_data;
        // bam->data = swap.data;
        if (!bam_copy1(buffer.bam, bam)) {
            ;
        }
        if (m_finalize || m_reads_buffer.empty()) {
            BamBlockList block(buffer, m_reads_idle, false);

            if (!m_reads_buffer.empty()) {
                m_reads_buffer.back().is_out_interval =
                    m_reads_buffer.back().end + skip_bps < buffer.bam->core.pos || m_reads_buffer.back().tid < buffer.bam->core.tid;
                m_reads_buffer.back().is_out_chrom = m_reads_buffer.back().tid < buffer.bam->core.tid;
            }
            m_reads_buffer.push_back(block);
        }
        else if (m_reads_buffer.back().end + skip_bps < buffer.bam->core.pos || m_reads_buffer.back().tid < buffer.bam->core.tid) {
            BamBlockList block(buffer, m_reads_idle, true);
            m_reads_buffer.back().is_out_interval = true;
            m_reads_buffer.back().is_out_chrom = m_reads_buffer.back().tid < buffer.bam->core.tid;
            m_reads_buffer.push_back(block);
            return true;
        }
        else {
            m_reads_buffer.back().insert(buffer, m_reads_idle);
        }
        m_finalize = false;
        return false;
    }
    void finalize()
    {  // 将无用的数据放回去
        m_finalize = true;
    }

    bam1_t swap;
    bool m_finalize;
    std::mutex m_mutex;
    std::list<BamBlockList> m_reads_buffer;
    std::list<RegionReadRecord> m_reads_idle;
};

class RegionManager
{
public:
    static constexpr int read_padded_span = 100;

    static constexpr int max_reads_size = 300;

    RegionManager();

    RegionManager(ActiveRegionEngine *engine_, ReferenceManager *fasta_loader, BlockingQueue<std::shared_ptr<RegionSource>> *region_queue,
                  ActiveRegionBamBlockListSource *block_resource, BlockingQueue<std::shared_ptr<RegionResource>> *region_resource,
                  bool force_non_active, DbsnpManager *db)
        : engine(engine_)
        , m_fasta_loader(fasta_loader)
        , m_region_queue(region_queue)
        , m_block_resource(block_resource)
        , m_region_resource(region_resource)
        , force_non_active_(force_non_active)
        , db_(db)
    {}

    void poll(bool force);

    void flush(bool force);

    void update_finish_sig();

    void calculate_db(std::shared_ptr<RegionSource>& source);

private:
    uint32_t buffer_threshold;
    ActiveRegionEngine *engine;
    ReferenceManager *m_fasta_loader;
    BlockingQueue<std::shared_ptr<RegionSource>> *m_region_queue;
    ActiveRegionBamBlockListSource *m_block_resource;
    BlockingQueue<std::shared_ptr<RegionResource>> *m_region_resource;
    std::shared_ptr<RegionResource> m_resource{nullptr};
    std::vector<bam1_t *> reads;
    std::vector<RegionResult> result;
    int last_tid{0};
    int region_id{0};
    int source_id{0};
    uint32_t max_used_region{0};
    bool force_non_active_;
    std::list<BamBlockList>::iterator block_itertor;
    std::list<RegionReadRecord>::iterator bam_itertor;

    DbsnpManager *db_;
    int32_t db_tid_{0};  // 若存在db，指示上一次使用的tid，若此次与上次不一致，代表上一条chr结束，置零db_offset同时调用db.notify()；
    int32_t db_offset_{0};  // 对db_tid_指示的数据有序递增索引，避免重复便利
};
#endif