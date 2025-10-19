#ifndef ROVACA_HC_WRITER_H_
#define ROVACA_HC_WRITER_H_
#include <chrono>
#include <cstddef>
#include <cstdint>
#include <memory>
#include <memory_resource>
#include <mutex>
#include <thread>
#include <vector>

#include "common/utils/object_pool.h"
#include "genotype/genotype_argument.h"
#include "genotype/genotype_macors.h"
#include "htslib/hts.h"
#include "htslib/kstring.h"

typedef struct bcf1_t bcf1_t;
typedef struct htsFile htsFile;
typedef struct bcf_hdr_t bcf_hdr_t;
typedef struct sam_hdr_t bam_hdr_t;

template <typename T>
class BlockingQueue;

namespace rovaca
{

typedef struct WriterTask
{
    int32_t tid;
    int32_t source_id;
    kstring_t* cache;
    size_t del_line_offset;
    int64_t next_available_start;
    std::shared_ptr<char> ref;
} WriterTask, *pWriterTask;

struct OutputFile
{
    htsFile* out_file;
    bcf_hdr_t* hdr;
};

#define CONTIG_CACHE_SIZE (64)
#define KS_CACHE_SIZE     (512)

class Writer
{
    /* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
private:
    bool index_;
    int32_t compression_level_;
    pHCArgs args_;
    bam_hdr_t* bam_hdr_;
    htsThreadPool* hts_pool_;
    BlockingQueue<pWriterTask>* result_queue_;

    OutputFile out_;
    int32_t last_id_back_;
    std::thread* merge_thread_;

    ObjectPool<kstring_t*> ks_pool_;

    char contig_cache[CONTIG_CACHE_SIZE]{};
    char pos_cache[CONTIG_CACHE_SIZE]{};

    bool pre_is_del_{false};
    int32_t pre_tid_{INVALID_INT};
    int64_t next_available_start_{INVALID_INT};
    kstring_t del_cache_{0, 0, nullptr};

    /* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
public:
    Writer(bool index, int32_t compression_level, pHCArgs args, bam_hdr_t* bam_hdr, htsThreadPool* tp,
           BlockingQueue<pWriterTask>* result_queue);

    ~Writer() = default;

    kstring_t* pop_bcf_cache();
    void push_bcf_cache(kstring_t* s);
    void close_file();
    htsFile* vcf_file_ptr() { return out_.out_file; }

    /*!
     * @brief 创建 Writer 对象后必须主动调用此函数才会创建相关线程
     */
    bool start();

    /*!
     * @brief 所有任务结束时，main 线程必须手动调用此函数回收线程资源
     */
    void stop();

    void update_last_id_back(int32_t id) { last_id_back_ = id; }

    bcf_hdr_t* bcf_header() const { return out_.hdr; }

    /*!
     * @brief 创建header，设置为静态，方便测试时不创建Writer直接创建hdr
     */
    static bcf_hdr_t* init_vcf_header(pHCArgs args, bam_hdr_t* bam_hdr);
    static bcf_hdr_t* init_gvcf_header(pHCArgs args, bam_hdr_t* bam_hdr);

    /* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
private:
    void merge_task_thread_call();
    void resolve_overrides_and_write_task(pWriterTask task);
    bool check_deletion_variant(const std::shared_ptr<char>& ref, kstring_t* s, size_t* current_offset);

    static bool is_gz_file(const char* filename);
    static std::string generate_gvcf_block_line(int32_t min_gq, int32_t max_gq);
    static std::string generate_contig_line(const char* chr_name, uint32_t len);
    static std::string generate_command_line(const char* tool_name, const char* command_line);

    static void generate_gvcf_block(bcf_hdr_t* hdr, pHCArgs args);
    static void generate_contig(bcf_hdr_t* hdr, bam_hdr_t* bam_hdr);
    static void generate_samples(bcf_hdr_t* hdr, bam_hdr_t* bam_hdr);
};

}  // namespace rovaca

#endif  // ROVACA_HC_WRITER_H_
