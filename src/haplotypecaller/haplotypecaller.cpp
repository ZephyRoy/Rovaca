#include "haplotypecaller.h"

#include <iostream>
#include <memory>
#include <string>
#include <thread>
#include <unordered_map>

#include "ActiveMainThread.h"
#include "common/enum.h"
#include "dbsnp_manager.h"
#include "downsampler_hc.h"
#include "pairhmm/pairhmm_engine.h"
#include "reads_filter_hc.h"
#include "reference_manager.h"
#include "ring_mem_pool.hpp"
#include "transformer.h"

// Register Haplotypecaller.
static bool HC = RovacaToolRegister::instance().register_tool("HaplotypeCaller", [] { return std::make_unique<HaplotypeCaller>(); });

static constexpr const size_t k_reads_buffer_mem = 8 * 1024 * 1024;
static constexpr const size_t k_wes_bamdata_pool_mem = 8 * 1024 * 1024;
static constexpr const size_t k_wgs_bamdata_pool_mem = 32 * 1024 * 1024;
static constexpr const size_t k_target_mem = 100 * 1024 * 1024;
static constexpr const size_t k_max_bqsr_tasks = 5;

HaplotypeCaller::HaplotypeCaller()
    : hc_args(nullptr)
    , hc_engine(nullptr)
    , assemble_output_thread_(nullptr)
    , debug_assemble_stream_(false)
{}

HaplotypeCaller::~HaplotypeCaller() {}

void HaplotypeCaller::clear_and_exit()
{
    // TODO: Close output files when signal received.
    RovacaLogger::error("exiting");
    writer_->close_file();
    _exit(-1);
}

void HaplotypeCaller::do_work()
{
    init_args();
    // 使用参数控制
    rovaca::init_pairhmm_ptr(rovaca_args_->old_pairhmm_engine());
    rovaca::GenotypeArgument *genotype_args = &hc_args->genotype_args_;

    int n_threads = rovaca_args_->run_pool_size();
    int n_resource = rovaca_args_->run_pool_size() * 2;
    bool force = genotype_args->emit_reference_confidence();
    bool has_db = !rovaca_args_->dbsnp_path().empty();

    ReferenceManager fasta_manager_inst(reference_path_, bed_loader_.get());
    ActiveRegionBamBlockListSource block_resource(10);

    RingMemPool<bam1_t> *bqsr_mem = nullptr;

    if (apply_transformer_) {
        bqsr_mem = new RingMemPool<bam1_t>();
        bqsr_mem->init(BAM_DATA_CAPACITY, 1024 * 1024 * 5);
        apply_bqsr_ = new BQSRReadTransformer(read_group_id_, rovaca_args_->recal_table(), bqsr_mem);
        RovacaLogger::info("Apply BQSR transformer");
    }

    DbsnpManager db_manager{};
    if (has_db) {
        db_manager.initialization(rovaca_args_->dbsnp_path(), bed_loader_.get(), &fasta_manager_inst.get_contig(), stream_pool_.get(), 10);
        RovacaLogger::info("Apply DBSNP");
    }

    BlockingQueue<std::shared_ptr<ActiveBaseResource>> *base_resource = new BlockingQueue<std::shared_ptr<ActiveBaseResource>>(n_resource);
    BlockingQueue<BaseSource> *base_source = new BlockingQueue<BaseSource>(n_resource);
    BlockingQueue<std::shared_ptr<RegionResource>> *region_source = new BlockingQueue<std::shared_ptr<RegionResource>>(n_resource);
    BlockingQueue<std::shared_ptr<RegionSource>> *region_queue = new BlockingQueue<std::shared_ptr<RegionSource>>(n_resource);
    BlockingQueue<pWriterTask> *result_queue = new BlockingQueue<pWriterTask>(2048);
    BlockingQueue<std::shared_ptr<BamSource>> *m_bam_resource = new BlockingQueue<std::shared_ptr<BamSource>>(n_resource + 20);

    std::thread resource_thread([&]() {
        for (int i = 0; i < n_threads; i++) {
            m_bam_resource->push(std::make_shared<BamSource>());
            std::shared_ptr<ActiveBaseResource> active_base_resource = std::make_shared<ActiveBaseResource>(1024, 8 * 1024 * 1024);
            base_resource->push(active_base_resource);
        }
        for (int i = 0; i < n_threads; i++) {
            m_bam_resource->push(std::make_shared<BamSource>());
        }
        for (int i = 0; i < n_resource; i++) {
            std::shared_ptr<RegionResource> region_resource = std::make_shared<RegionResource>(
                i, wes_ ? k_wes_bamdata_pool_mem : k_wgs_bamdata_pool_mem, k_reads_buffer_mem, k_target_mem);
            region_source->push(region_resource);
        }
    });
    boost::asio::thread_pool run_pool(n_threads);
    std::thread *bqsr_executor = nullptr;
    std::thread *bqsr_reducer = nullptr;
    BlockingQueue<std::shared_ptr<TransformTask>> *bqsr_task_queue = new BlockingQueue<std::shared_ptr<TransformTask>>(512);
    ReducedQueue *bqsr_processed_reads = new ReducedQueue();
    ReadStream *bqsr_streamer = nullptr;

    if (apply_transformer_) {
        bqsr_executor = new std::thread(&BQSRReadTransformer::launch, apply_bqsr_, streamer_.get(), &run_pool, bqsr_task_queue);
        bqsr_reducer = new std::thread(&BQSRReadTransformer::reduce, apply_bqsr_, bqsr_task_queue, bqsr_processed_reads);
        bqsr_streamer = new ReadsTransformerIterator(streamer_.release(), bqsr_processed_reads, bqsr_mem);
    }
    auto streamer = apply_bqsr_ ? bqsr_streamer : streamer_.get();
    std::thread fasta_process(&ReferenceManager::run, &fasta_manager_inst);
    ActiveMainThreadDispatchTasks main_dispatch_process(streamer, bed_loader_.get(), &fasta_manager_inst, &fasta_manager_inst.get_contig(),
                                                        &block_resource, &run_pool, base_resource, base_source, m_bam_resource,
                                                        apply_bqsr_);
    ActiveMainThreadReduce main_reduce_process(&fasta_manager_inst, &fasta_manager_inst.get_contig(), bed_loader_.get(), region_source,
                                               base_resource, &block_resource, base_source, region_queue, force,
                                               has_db ? &db_manager : nullptr);
    std::thread thread_dispatch(&ActiveMainThreadDispatchTasks::run, &main_dispatch_process);
    std::thread thread_reduce(&ActiveMainThreadReduce::run, &main_reduce_process);

    std::unique_ptr<std::thread> db_process;
    if (has_db) {
        db_process = std::make_unique<std::thread>(&DbsnpManager::run, &db_manager);
    }

    bool index = rovaca_args_->create_output_index();
    int32_t compression_level = rovaca_args_->compression_level();
    writer_ = std::make_unique<Writer>(index, compression_level, genotype_args, merged_header_, stream_pool_.get(), result_queue);
    writer_->start();
    k_outfile = writer_->vcf_file_ptr();  // TODO: writer 初始化放框架层
    hc_engine =
        std::make_unique<HaplotypeCallerEngine>(hc_args.get(), region_source, region_queue, result_queue, &run_pool, &fasta_manager_inst,
                                                merged_header_, writer_.get(), samplelist_, debug_assemble_stream_);

    std::thread thread_consumer(&HaplotypeCallerEngine::call_region, hc_engine.get());

    pthread_setname_np(fasta_process.native_handle(), "Reference");
    pthread_setname_np(thread_dispatch.native_handle(), "RegionDispatch");
    pthread_setname_np(thread_reduce.native_handle(), "RegionReduce");
    pthread_setname_np(thread_consumer.native_handle(), "CallRegion");
    if (has_db) {
        pthread_setname_np(db_process->native_handle(), "DBSNP");
    }
    if (apply_transformer_) {
        pthread_setname_np(bqsr_executor->native_handle(), "BQSRLaunch");
        pthread_setname_np(bqsr_reducer->native_handle(), "BQSRReduce");
    }

    if (apply_transformer_) {
        bqsr_executor->join();
        bqsr_reducer->join();
    }

    fasta_process.join();
    if (has_db) {
        db_process->join();
    }
    thread_dispatch.join();
    thread_reduce.join();
    thread_consumer.join();
    resource_thread.join();
    run_pool.join();
    writer_->stop();

    if (has_db) {
        db_manager.terminalization();
    }

    delete base_resource;
    delete base_source;
    delete region_source;
    delete region_queue;
    delete result_queue;
    if (apply_transformer_) {
        delete apply_bqsr_;
        delete bqsr_task_queue;
        delete bqsr_processed_reads;
        delete bqsr_streamer;
        bqsr_mem->term();
        delete bqsr_mem;
    }
}

const UniqueStream &HaplotypeCaller::custom_streamer()
{
    // return a customted reads stream.
    if (apply_filter_) {
        filter_ = std::make_unique<HCReadFilter>(merged_header_, rovaca_args_->inspect_reads());
        streamer_.reset(new ReadsFilterIterator(streamer_.release(), filter_.get()));
    }
    // 先降采样再重校准质量.
    if (apply_downsampler_) {
        downsampler_ = std::make_unique<HCDownsampler>(downsample_threshold_);
        streamer_.reset(new ReadsDownsampleIterator(streamer_.release(), downsampler_.get()));
    }

    return streamer_;
}

void HaplotypeCaller::init_args()
{
    hc_args = std::make_unique<HaplotypeCallerArgs>();
    hc_args->initGQBands(rovaca_args_->gq_bands());
    hc_args->pcrIndelModel = s_str2pcr.at(rovaca_args_->pcr_model());
    hc_args->referenceConfidenceMode = s_str2erc.at(rovaca_args_->erc_model());
    hc_args->pairhmm_arguments_.base_quality_score_threshold = rovaca_args_->base_quality_score_threshold();
    hc_args->pairhmm_arguments_.pcr_option = hc_args->pcrIndelModel;

    debug_assemble_stream_ = rovaca_args_->assemble_output_path() != nullptr;

    rovaca::GenotypeArgument *genotype_args = &hc_args->genotype_args_;
    genotype_args->output = rovaca_args_->vcf_path().c_str();
    genotype_args->tool_name = rovaca_args_->tool().c_str();
    genotype_args->command_line = rovaca_args_->command_line();
    genotype_args->gvcf_gq_bands = hc_args->GQBands_;
    genotype_args->init_reference_confidence_mode(hc_args->referenceConfidenceMode);

    if (hc_args->referenceConfidenceMode == ReferenceConfidenceMode::GVCF) {
        std::string gvcf_gq_bands;
        if (rovaca_args_->gq_bands().empty()) {
            gvcf_gq_bands = "[1,2,3,4,...59,60,70,80,90,99]";
        }
        else {
            gvcf_gq_bands.append("[");
            for (const auto &band : rovaca_args_->gq_bands()) {
                gvcf_gq_bands.append(std::to_string(band)).append(1, ',');
            }
            gvcf_gq_bands.back() = ']';
        }

        fprintf(stdout,
                "%s startup parameters:\n"
                "    version                        = %s\n"
                "    threads                        = %d\n"
                "    gvcf-gq-bands                  = %s\n"
                "    max-reads-depth                = %d\n"
                "    pcr-indel-model                = %s\n"
                "    inspect-reads                  = %s\n"
                "    pairhmm-engine                 = %s\n"
                "    compression-level              = %d\n"
                "    interval-padding               = %d\n"
                "    emit-ref-confidence            = %s\n"
                "    base-quality-score-threshold   = %d\n\n",
                genotype_args->tool_name, MAIN_VERSION, rovaca_args_->run_pool_size(), gvcf_gq_bands.c_str(), rovaca_args_->max_reads_depth(),
                pcr2str(hc_args->pairhmm_arguments_.pcr_option).c_str(), (rovaca_args_->inspect_reads() ? "true" : "false"),
                (rovaca_args_->old_pairhmm_engine() ? "intel" : "rovaca"), rovaca_args_->compression_level(), rovaca_args_->interval_padding(),
                erc2str(hc_args->referenceConfidenceMode).c_str(), rovaca_args_->base_quality_score_threshold());
    }
    else if (hc_args->referenceConfidenceMode == ReferenceConfidenceMode::NONE) {
        fprintf(stdout,
                "%s startup parameters:\n"
                "    version                        = %s\n"
                "    threads                        = %d\n"
                "    max-reads-depth                = %d\n"
                "    pcr-indel-model                = %s\n"
                "    inspect-reads                  = %s\n"
                "    pairhmm-engine                 = %s\n"
                "    compression-level              = %d\n"
                "    interval-padding               = %d\n"
                "    emit-ref-confidence            = %s\n"
                "    base-quality-score-threshold   = %d\n\n",
                genotype_args->tool_name, MAIN_VERSION, rovaca_args_->run_pool_size(), rovaca_args_->max_reads_depth(),
                pcr2str(hc_args->pairhmm_arguments_.pcr_option).c_str(), (rovaca_args_->inspect_reads() ? "true" : "false"),
                (rovaca_args_->old_pairhmm_engine() ? "intel" : "rovaca"), rovaca_args_->compression_level(), rovaca_args_->interval_padding(),
                erc2str(hc_args->referenceConfidenceMode).c_str(), rovaca_args_->base_quality_score_threshold());
    }
    fflush(stdout);
}
