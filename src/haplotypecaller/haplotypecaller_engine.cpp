#include "haplotypecaller_engine.h"

#include <algorithm>
#include <boost/asio.hpp>
#include <cstdio>
#include <iostream>
#include <map>
#include <unordered_map>
#include <vector>

#include "common/utils/rovaca_memory_pool.h"
#include "genotype/allele_likelihoods.hpp"
#include "genotype/block_combiner.h"
#include "genotype/genotype.h"
#include "genotype/genotypes_context.hpp"
#include "genotype/germline_genotying_engine.h"
#include "genotype/homogeneous_ploidy_model.hpp"
#include "genotype/indexed_allele_list.hpp"
#include "genotype/indexed_sample_list.hpp"
#include "genotype/read_record.h"
#include "genotype/utils/adapter_utils.h"
#include "genotype/utils/assembly_based_caller_utils.h"
#include "genotype/utils/debug_utils.h"
#include "genotype/variant.h"
#include "rovaca_logger.h"
#include "pairhmm/pairhmm_engine.h"

static constexpr int k_region_padding = 100;
static constexpr int k_reference_padding = 500;
static constexpr int32_t k_minimum_read_length_after_trimming = 10;
static constexpr int32_t k_read_length_filter_threshold = 10;

#define CLEAR_ERSOURCE()                                                                                         \
    do {                                                                                                         \
        hc_apply_reset(assembler);                                                                               \
        resource->pool_->reset(mark_size);                                                                       \
        std::for_each(extra_memory_reads.begin(), extra_memory_reads.end(), [](bam1_t *b) { bam_destroy1(b); }); \
    } while (false)

HaplotypeCallerEngine::HaplotypeCallerEngine(HaplotypeCallerArgs *hc_args, BlockingQueue<std::shared_ptr<RegionResource>> *region_source,
                                             BlockingQueue<std::shared_ptr<RegionSource>> *region_queue,
                                             BlockingQueue<pWriterTask> *result_queue, boost::asio::thread_pool *thread_pool,
                                             ReferenceManager *reference_manager, bam_hdr_t *hdr, Writer *writer,
                                             const std::vector<std::string> &samples, bool debug_output)
    : region_source_(region_source)
    , region_queue_(region_queue)
    , result_queue_(result_queue)
    , run_pool_(thread_pool)
    , writer_(writer)
    , hdr_(hdr)
    , vcf_hdr_(writer->bcf_header())
    , reference_manager_(reference_manager)
    , hc_args_(hc_args)
    , assemble_argumets_(hc_args->assemble_arguments_)
    , debug_assemble_(debug_output)
    , sample_list_(nullptr)
    , ploidy_model_(nullptr)
    , genotype_engine_pool_(new ObjectPool<pGermlineGenotyingEngine>{})
{
    assemble_argumets_.read_threading_argument.kmer = {10, 25};
    assemble_argumets_.debugAssembly = false;
    AssembleEngine::init_assemble_argument(&assemble_argumets_);
    assemble_results_ = new BlockingQueue<AssembleDebugOutput>(1024);
    assembler_pool_ = new ObjectPool<p_hc_apply>();
    sample_list_ = IndexedSampleList::create(samples);
    ploidy_model_ = HomogeneousPloidyModel::create(hc_args->genotype_args_.sample_ploidy, sample_list_);
}

HaplotypeCallerEngine::~HaplotypeCallerEngine()
{
    for (auto &assembler : assembler_pool_->obj_list()) {
        AssembleEngine::destory_assembler(assembler);
    }
    assembler_pool_->obj_list().clear();
    delete sample_list_;
    delete ploidy_model_;
    delete genotype_engine_pool_;
    delete assembler_pool_;
    AssembleEngine::finit_assemble_argument();
    delete assemble_results_;
}

// todo: 下个版本优化空read的区间直接出结果
static pVariant emptyVariant(pSimpleInterval original, uint8_t ref, pMemoryPool pool)
{
    pVariant vc = Variant::create(pool);
    vc->set_tid(original->get_tid());
    vc->set_start(original->get_start());
    vc->set_stop(original->get_stop());
    vc->set_end_key(original->get_stop());
    pAllele ref_a = Allele::create_allele(ref, 1);
    AlleleVector vc_alleles{{ref_a, StaticAllele::get_instance()->_non_ref_allele.get()}, pool};
    vc->set_alleles(vc_alleles);

    pGenotype g = Genotype::create(pool);
    g->set_alleles({{ref_a, ref_a}, pool});
    g->set_dp(0);
    g->set_gq(0);
    g->set_pl({{0, 0, 0}, pool});
    pGenotypesContext gc = GenotypesContext::create(pool);
    gc->add(g);
    vc->set_genotype(gc);

    return vc;
}

void HaplotypeCallerEngine::call_region()
{
    bool is_avx512 = avx512_supported();

    while (true) {
        std::shared_ptr<RegionSource> source = region_queue_->pop();

        if (source->finish_) {
            writer_->update_last_id_back(source->source_id_);
            break;
        }

        boost::asio::post(*run_pool_, [=]() mutable {
            pthread_setname_np(pthread_self(), "Assemble");
            std::shared_ptr<RegionResource> resource = source->resource_;

            uint8_t *buffer = resource->assemble_result_buffer_;
            size_t buffer_size = resource->result_buffer_size_;
            p_hc_apply assembler = assembler_pool_->pop();
            if (assembler == nullptr) assembler = AssembleEngine::creat_assembler();
            pGermlineGenotyingEngine genotype_engine = genotype_engine_pool_->pop();
            if (nullptr == genotype_engine) genotype_engine = new GermlineGenotyingEngine{};
            genotype_engine->clear_upstream_deletions_loc();
            genotype_engine->set_dbsnp(source->db_offset_, source->db_data_.get());

            p_lib_sw_avx sw = genotype_engine->sw();
            PcrIndelModel pcr = hc_args_->pairhmm_arguments_.pcr_option;
            bool emit_reference_confidence = hc_args_->genotype_args_.emit_reference_confidence();
            int32_t base_quality = hc_args_->pairhmm_arguments_.base_quality_score_threshold;

            uint32_t mark_size;
            rovaca::RefFragment ref_bases;
            rovaca::pBlockCombiner block_combiner{nullptr};
            pSimpleInterval ref_loc, original, original_padded, variant, variant_padded;
            GenotypeArgument *gargs = &hc_args_->genotype_args_;

            kstring_t *bcf_cache = source->regions_.empty() ? nullptr : writer_->pop_bcf_cache();

            for (const auto &result : source->regions_) {
                rovaca::RovacaMemoryPool target_mem{buffer, buffer_size};
                std::pmr::list<bam1_t *> extra_memory_reads{&target_mem};

                VariantVector region_result;
                mark_size = resource->pool_->get_used_data_size();
                genotype_engine->init_engine_per_loop(gargs, &target_mem, resource->pool_, hdr_, vcf_hdr_, sample_list_, ploidy_model_);
                block_combiner = genotype_engine->block_combiner();

                original = get_original_region(result.region, &target_mem);
                original_padded = get_original_padded_region(original, &target_mem);
                ref_loc = get_ref_padded_region(original_padded, &target_mem);
                ref_bases.data = (uint8_t *)source->ref_.get() + ref_loc->get_start() - 1;
                ref_bases.len = uint32_t(ref_loc->get_length());

                // todo: 下个版本优化空read的区间直接出结果
                if (emit_reference_confidence && result.padding_reads.empty()) {
                    uint8_t ref_byte = source->ref_.get()[original->get_start() - 1];
                    pVariant vc = emptyVariant(original, ref_byte, &target_mem);
                    block_combiner->submit(vc, bcf_cache);
                    CLEAR_ERSOURCE();
                    continue;
                }

                /*!
                 * @brief 几种返回情况
                 * assemble 之前
                 *      non active          -> reads to finalized
                 *      reads empty         -> reads to finalized
                 * assemble 之后
                 *      null variantSpan    -> trim region
                 *      haplotypes empty    -> trim haplotypes
                 *      reads empty         -> 经过两轮过滤以后
                 */
                AssembleResult *untrimed_result = run_local_assemble(assembler, result, source, &target_mem);
                ReadHashSet original_reads{{untrimed_result->get_reads().begin(), untrimed_result->get_reads().end()}, &target_mem};
                if (!result.region->active || untrimed_result->get_reads().empty() || untrimed_result->get_haplotypes().size() <= 1) {
                    if (emit_reference_confidence) {
                        region_result = genotype_engine->reference_model_for_no_variation(&ref_bases, ref_loc, original, original_padded,
                                                                                          gargs->sample_ploidy, original_reads);
                        block_combiner->submit_vector(region_result, bcf_cache);
                    }
                    CLEAR_ERSOURCE();
                    continue;
                }

                // trim_region 中包含了计算 all_variation_events
                HaplotypeVector &untrim_h = const_cast<HaplotypeVector &>(untrimed_result->get_haplotypes());
                auto trim_result = AdapterUtils::trim_region(untrim_h, &ref_bases, ref_loc, original, original_padded, gargs, &target_mem);
                if (trim_result.first == nullptr) {
                    if (emit_reference_confidence) {
                        region_result = genotype_engine->reference_model_for_no_variation(&ref_bases, ref_loc, original, original_padded,
                                                                                          gargs->sample_ploidy, original_reads);
                        block_combiner->submit_vector(region_result, bcf_cache);
                    }
                    CLEAR_ERSOURCE();
                    continue;
                }

                variant = trim_result.first;
                variant_padded = trim_result.second;
                ReadHashSet trimed_reads =
                    AdapterUtils::trim_reads_by_region(original_reads, variant_padded, &target_mem, resource->pool_, extra_memory_reads);
                HaplotypeVector trimed_haps = AdapterUtils::trim_haplotype_by_region(untrim_h, trim_result.second, &target_mem);
                if (trimed_haps.size() <= 1) {
                    if (emit_reference_confidence) {
                        region_result = genotype_engine->reference_model_for_no_variation(&ref_bases, ref_loc, original, original_padded,
                                                                                          gargs->sample_ploidy, original_reads);
                        block_combiner->submit_vector(region_result, bcf_cache);
                    }
                    CLEAR_ERSOURCE();
                    continue;
                }

                // 对assembly_result中的reads进行第一次过滤，此次过滤掉的reads不再需要
                ReadHashSet passing_reads{&target_mem};
                std::for_each(trimed_reads.begin(), trimed_reads.end(), [&](pReadRecord r) {
                    if (r->unclipped_read_length() >= k_minimum_read_length_after_trimming) {
                        passing_reads.insert(r);
                    }
                });

                // 对passing_reads进行第二次过滤
                // 通过条件的per_sample_passing_read_list传递给GermlineGenotypingEngine
                // 未通过条件的per_sample_filtered_read_list传递给pairhmm
                Int32ToReadVectorMap per_sample_passing_read_list = filter_non_passing_reads(passing_reads, &target_mem);
                Int32ToReadVectorMap per_sample_filtered_read_list{{{0, {passing_reads.begin(), passing_reads.end()}}}, &target_mem};
                if (per_sample_filtered_read_list.at(0).empty()) {
                    if (emit_reference_confidence) {
                        region_result = genotype_engine->reference_model_for_no_variation(&ref_bases, ref_loc, original, original_padded,
                                                                                          gargs->sample_ploidy, original_reads);
                        block_combiner->submit_vector(region_result, bcf_cache);
                    }
                    CLEAR_ERSOURCE();
                    continue;
                }

                // pairhmm
                pHaplotype ref_haplotype = trimed_haps.front();
                ReadVector &trimed_reads2 = per_sample_filtered_read_list.at(0);
                if (is_avx512) {
                    std::sort(trimed_reads2.begin(), trimed_reads2.end(),
                              [](rovaca::pReadRecord l, rovaca::pReadRecord r) { return l->seq_length() < r->seq_length(); });
                }

                DoubleVector2D likelihoods = rovaca::call_pairhmm(trimed_haps, trimed_reads2, base_quality, pcr, &target_mem);

                // realignReadsToTheirBestHaplotype
                auto *alleles = IndexedAlleleList<pHaplotype>::create(trimed_haps, &target_mem);
                ReadVector2D evidence_by_sample{{trimed_reads2}, &target_mem};
                ReadVector2D filtered_evidence_by_sample{&target_mem};
                std::for_each(per_sample_passing_read_list.begin(), per_sample_passing_read_list.end(),
                              [&](Int32ToReadVectorMap::value_type itr) { filtered_evidence_by_sample.push_back(std::move(itr.second)); });
                DoubleVector3D likelihoods_by_sample{{likelihoods}, &target_mem};

                auto *rh_likelihoods = RHLikelihoods::create<pReadRecord, pHaplotype>(
                    &target_mem, sample_list_, alleles, std::move(evidence_by_sample), std::move(filtered_evidence_by_sample),
                    std::move(likelihoods_by_sample));

                AssemblyBasedCallerUtils::realign_reads_to_their_best_haplotype(rh_likelihoods, ref_haplotype, ref_loc->get_start(), sw,
                                                                                &target_mem, resource->pool_);
                // genotype
                auto variants = genotype_engine->assign_genotype_likelihoods(rh_likelihoods, &ref_bases, ref_loc, variant,
                                                                             per_sample_passing_read_list);
                if (emit_reference_confidence) {
                    if (!GermlineGenotyingEngine::contains_calls(variants.first)) {
                        region_result = genotype_engine->reference_model_for_no_variation(&ref_bases, ref_loc, original, original_padded,
                                                                                          gargs->sample_ploidy, original_reads);
                    }
                    else {
                        const ReadVector &realign_reads = rh_likelihoods->sample_evidence(0);
                        ReadHashSet genotype_reads{{realign_reads.begin(), realign_reads.end()}, &target_mem};
                        region_result = genotype_engine->call_non_active_site(ref_haplotype, &ref_bases, ref_loc, original, original_padded,
                                                                              variant, variant_padded, variants, original_reads,
                                                                              genotype_reads, extra_memory_reads);
                    }
                }
                const VariantVector &res = emit_reference_confidence ? region_result : variants.first;
                block_combiner->submit_vector(res, bcf_cache);
                CLEAR_ERSOURCE();
            }

            pWriterTask task = new WriterTask{};
            if (block_combiner) {
                task->del_line_offset = block_combiner->del_line_offset();
                task->next_available_start = block_combiner->next_variant_pos();
                block_combiner->force_output(bcf_cache);
            }
            task->source_id = source->source_id_;
            task->cache = bcf_cache;
            task->ref = source->ref_;
            task->tid = source->tid_;
            result_queue_->push(task);

            genotype_engine_pool_->push(genotype_engine);
            hc_apply_reset(assembler);
            assembler_pool_->push(assembler);
            source->release();
            resource->finalize();
            region_source_->push(resource);
        });
    }

    reference_manager_->assign_finish();
}

AssembleResult *HaplotypeCallerEngine::run_local_assemble(p_hc_apply assembler, const RegionResult &result,
                                                          std::shared_ptr<RegionSource> source, pMemoryPool target_mem)
{
    p_hc_region_active_storage region = get_assemble_region(result.region);
    const std::vector<bam1_t *> &region_reads = result.padding_reads;
    // printf("ActiveRegion\t%s:%lu-%lu\tReads\t%ld\n", contig, region->start_index, region->end_index, region_reads.size());
    // Call local asemble.
    const uint8_t *ref = reinterpret_cast<const uint8_t *>(source->ref_.get());
    AssembleReadsBuffer reads_mem{buffer_ : source->resource_->buffer_, used_ : 0, capacity_ : (uint32_t)source->resource_->buffer_size_};
    AssembleResult *untrimed_result =
        AssembleEngine::local_assemble(assembler, region, ref, source->ref_len_, region_reads, target_mem, &reads_mem);
    if (debug_assemble_) {
        assemble_results_->push({untrimed_result ? new AssembleResult(*untrimed_result) : nullptr, result.region_id});
        if (__glibc_unlikely(source->finish_)) assemble_results_->push({nullptr, INT_MAX});
    }
    return untrimed_result;
}
p_hc_region_active_storage HaplotypeCallerEngine::get_assemble_region(p_hc_region_active_storage region_result)
{
    p_hc_region_active_storage region = region_result;
    region->start_index += 1;
    region->end_index += 1;
    region->activeSpan = {region->start_index, region->end_index};
    region->paddedSpan = {std::max(int64_t(region->start_index - k_region_padding), int64_t(1)),
                          std::min(int64_t(region->end_index + k_region_padding), int64_t(hdr_->target_len[region->tid]))};
    return region;
}

void HaplotypeCallerEngine::assemble_debug_stream(const char *assemble_output_path)
{
    std::fstream assemble_debug_streamer(assemble_output_path, std::ofstream::out);
    std::unordered_map<int, AssembleResult *> assemble_wait_hash;
    using HashItr = std::unordered_map<int, AssembleResult *>::iterator;
    int region_count = 0;
    while (true) {
        AssembleDebugOutput assemble_debug_outputs = assemble_results_->pop();
        if (assemble_debug_outputs.id == INT_MAX) {
            break;
        };
        while (true) {
            HashItr itr = assemble_wait_hash.find(region_count);
            if (itr != assemble_wait_hash.end()) {
                auto *untrimed_result = itr->second;
                if (!untrimed_result) {
                    region_count++;
                    assemble_wait_hash.erase(itr);
                    continue;
                }
                // Get Assemble results.
                auto &reads = const_cast<std::pmr::vector<rovaca::pReadRecord> &>(untrimed_result->get_reads());
                auto &haplotypes = const_cast<std::pmr::vector<rovaca::pHaplotype> &>(untrimed_result->get_haplotypes());
                auto region = untrimed_result->get_region();
                // Sort in Lexicographic order for convenient compair.
                std::sort(haplotypes.begin(), haplotypes.end(), [](rovaca::pHaplotype a, rovaca::pHaplotype b) {
                    return strcmp((const char *)(a->get_display_string()->data), (const char *)(b->get_display_string()->data)) < 0;
                });
                std::sort(reads.begin(), reads.end(), [](rovaca::pReadRecord a, rovaca::pReadRecord b) {
                    return a->assemble_read()->pos_start < b->assemble_read()->pos_start;
                });
                if (assemble_debug_streamer.good()) {
                    assemble_debug_streamer << sam_hdr_tid2name(hdr_, region->activeSpan->get_tid()) << ":"
                                            << region->activeSpan->get_start() << "-" << region->activeSpan->get_stop() << "\t";
                    assemble_debug_streamer << "reads: " << reads.size() << "\t";
                    assemble_debug_streamer << "haplotypes: " << haplotypes.size() << std::endl;
                    for (auto r : reads) {
                        assemble_debug_streamer << r->assemble_read()->read_data << "\t" << r->assemble_read()->pos_start << "\n";
                    }
                    for (auto hap : haplotypes) {
                        assemble_debug_streamer << hap->get_display_string()->data << "\n";
                    }
                    assemble_debug_streamer << std::endl;
                    assemble_debug_streamer.flush();
                }
                region_count++;
                delete itr->second;
                assemble_wait_hash.erase(itr);
            }
            else {
                break;
            }
        }
        assemble_wait_hash.insert({assemble_debug_outputs.id, assemble_debug_outputs.result});
    }
    assemble_debug_streamer.close();
}

pSimpleInterval HaplotypeCallerEngine::get_original_region(p_hc_region_active_storage region, pMemoryPool pool)
{
    int64_t stop = std::min(region->end_index + 1, int64_t(hdr_->target_len[region->tid]));
    return SimpleInterval::create(region->tid, region->start_index + 1, stop, pool);
}

pSimpleInterval HaplotypeCallerEngine::get_original_padded_region(pSimpleInterval original, pMemoryPool pool)
{
    int32_t tid = original->get_tid();
    int64_t start = std::max(original->get_start() - k_region_padding, int64_t(1));
    int64_t stop = std::min(original->get_stop() + k_region_padding, int64_t(hdr_->target_len[original->get_tid()]));
    return SimpleInterval::create(tid, start, stop, pool);
}

pSimpleInterval HaplotypeCallerEngine::get_ref_padded_region(pSimpleInterval original_padded, pMemoryPool pool)
{
    int32_t tid = original_padded->get_tid();
    int64_t start = std::max(original_padded->get_start() - k_reference_padding, int64_t(1));
    int64_t stop = std::min(original_padded->get_stop() + k_reference_padding, int64_t(hdr_->target_len[original_padded->get_tid()]));
    return SimpleInterval::create(tid, start, stop, pool);
}

Int32ToReadVectorMap HaplotypeCallerEngine::filter_non_passing_reads(ReadHashSet &reads, pMemoryPool pool)
{
    ReadVector filtered_reads{pool};
    for (pReadRecord r : reads) {
        if (r->unclipped_read_length() < k_read_length_filter_threshold || r->mapping_quality() < hc_args_->mappingQualityThreshold ||
            !r->mate_on_same_contig_or_no_mapped_mate_read()) {
            filtered_reads.push_back(r);
            reads.erase(r);
        }
    }

    // hc 单样本, 此处直接指定为0
    rovaca::Int32ToReadVectorMap per_sample_passing_read_list{{{int32_t(0), {}}}, pool};
    for (pReadRecord r : filtered_reads) {
        per_sample_passing_read_list.at(0).push_back(r);
    }

    return per_sample_passing_read_list;
}
