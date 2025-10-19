#ifndef BQSR_READ_TRANSFORMER_H
#define BQSR_READ_TRANSFORMER_H

#include <htslib/sam.h>

#include <atomic>
#include <boost/asio.hpp>
#include <boost/asio/thread_pool.hpp>
#include <boost/lockfree/queue.hpp>
#include <condition_variable>
#include <list>
#include <memory>
#include <memory_resource>
#include <mutex>
#include <string>
#include <thread>
#include <vector>

#include "BlockingQueue.h"
#include "bqsr_read_covarivates.h"
#include "bqsr_recal_tables.h"
#include "reads_stream.h"
#include "ring_mem_pool.hpp"

typedef struct TransformTask
{
    std::shared_ptr<std::list<bam1_t*>> batch_reads_;
    uint32_t task_id;
} TransformTask, *pTransformTask;
class BQSRReadTransformer
{
private:
    int preserve_q_less_than;
    double global_qual_score_prior;
    bool emit_original_quals;
    int total_covariate_count;
    int special_covariate_count;
    bool use_origin_base_quals;
    uint8_t* static_quantized_mapping;
    std::vector<std::string> read_groups_;
    std::unique_ptr<RecalArguments> recal_arguments_;
    std::unique_ptr<QuantizationInfo> quantization_info_;
    std::unique_ptr<RecalibrationTables> recalibration_tables_;
    std::unique_ptr<StandardCovariates> covariates_;
    BlockingQueue<std::shared_ptr<TransformTask>>* tasks_resource_;
    std::atomic_bool reads_done_{false};
    RingMemPool<bam1_t>* mem_;

public:
    BQSRReadTransformer(const std::vector<std::string>& read_groups, const std::string& recal_table_file, RingMemPool<bam1_t>* mem);
    ~BQSRReadTransformer();
    void apply(bam1_t* read);
    void load_report(const std::string& recal_table_file);
    void launch(ReadStream* stream, boost::asio::thread_pool* thread_pool, BlockingQueue<std::shared_ptr<TransformTask>>* task_queue);
    void reduce(BlockingQueue<std::shared_ptr<TransformTask>>* task_queue, ReducedQueue* processed_reads);

private:
    void compute_covariates(bam1_t* read, ReadCovariates* target_covariates);
    double hierarchical_bayesian_quality_estimate(double epsilon, const RecalDatum& empiricalQualRG, const RecalDatum& empiricalQualQS,
                                                  const RecalDatum& empirical_qual_cycle, const RecalDatum& empirical_qual_context);
    void initialize_argument_table(const BQSRReport& report);
    void initialize_quantization_table(const BQSRReport& report);
    void initialize_standard_covariate_list();
    void initialize_recalibration_tables(const std::map<std::string, BQSRReport>& tables);
    void parse_read_group_table(RecalDatumArray& array, const BQSRReport& table);
    void parse_quality_score_table(RecalDatumArray& array, const BQSRReport& table);
    void parse_external_table(const BQSRReport& table);
};

#endif  // BQSR_READ_TRANSFORMER_H
