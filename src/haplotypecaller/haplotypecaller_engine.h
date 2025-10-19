#ifndef HAPLOTYPECALLER_ENGINE_H
#define HAPLOTYPECALLER_ENGINE_H

#include <boost/asio/thread_pool.hpp>
#include <fstream>
#include <future>

#include "BlockingQueue.h"
#include "assemble_engine.h"
#include "common/utils/object_pool.h"
#include "genotype/forward.h"
#include "haplotypecaller_args.h"
#include "reference_manager.h"
#include "region_interface.h"
#include "writer/writer.h"

class HaplotypeCallerEngine
{
private:
    struct AssembleDebugOutput
    {
        AssembleResult *result;
        int id;
    };
    BlockingQueue<std::shared_ptr<RegionResource>> *region_source_;
    BlockingQueue<std::shared_ptr<RegionSource>> *region_queue_;
    BlockingQueue<pWriterTask> *result_queue_;
    BlockingQueue<AssembleDebugOutput> *assemble_results_;
    boost::asio::thread_pool *run_pool_;
    Writer *writer_;
    bam_hdr_t *hdr_;
    bcf_hdr_t *vcf_hdr_;
    ReferenceManager *reference_manager_;
    std::fstream assemble_debug_streamer_;
    HaplotypeCallerArgs *hc_args_;
    AssembleArgument assemble_argumets_;
    bool debug_assemble_;
    pInterfaceSampleList sample_list_;
    pInterfacePloidyModel ploidy_model_;
    ObjectPool<pGermlineGenotyingEngine> *genotype_engine_pool_;
    ObjectPool<p_hc_apply> *assembler_pool_;

public:
    HaplotypeCallerEngine(HaplotypeCallerArgs *hc_args, BlockingQueue<std::shared_ptr<RegionResource>> *region_source,
                          BlockingQueue<std::shared_ptr<RegionSource>> *region_queue, BlockingQueue<pWriterTask> *result_queue,
                          boost::asio::thread_pool *thread_pool, ReferenceManager *reference_manager, bam_hdr_t *hdr, Writer *writer,
                          const std::vector<std::string> &samples, bool debug_output);
    ~HaplotypeCallerEngine();
    void call_region();
    void assemble_debug_stream(const char *assemble_output_path);

private:
    AssembleResult *run_local_assemble(p_hc_apply assembler, const RegionResult &result, std::shared_ptr<RegionSource> source,
                                       pMemoryPool target_mem);
    p_hc_region_active_storage get_assemble_region(p_hc_region_active_storage region);

    pSimpleInterval get_original_region(p_hc_region_active_storage region, pMemoryPool pool);
    pSimpleInterval get_original_padded_region(pSimpleInterval original, pMemoryPool pool);
    pSimpleInterval get_ref_padded_region(pSimpleInterval original_padded, pMemoryPool pool);

    Int32ToReadVectorMap filter_non_passing_reads(ReadHashSet &reads, pMemoryPool pool);
};

#endif  // HAPLOTYPECALLER_ENGINE_H