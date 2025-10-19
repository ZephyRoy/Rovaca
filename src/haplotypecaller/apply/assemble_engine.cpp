#include "assemble_engine.h"

#include <iostream>

#include "haplotype.h"
#include "simple_interval.h"
#include "read_record.h"
#include "utils/base_utils.h"

using namespace rovaca;

p_hc_apply AssembleEngine::creat_assembler() { return hc_apply_init(); }
void AssembleEngine::destory_assembler(p_hc_apply assembler) { hc_apply_finit(assembler); }
AssembleResult* AssembleEngine::local_assemble(p_hc_apply assembler, p_hc_region_active_storage region, const uint8_t* chr_ref,
                                               uint32_t chr_ref_len, const std::vector<bam1_t*>& region_reads, pMemoryPool mem,
                                               AssembleReadsBuffer* reads_mem)
{
    assembler->reads_buffer_mem = reads_mem;
    hc_apply_set_target(assembler, region, const_cast<uint8_t*>(chr_ref), chr_ref_len);
    for (auto one_read : region_reads) {
        p_hc_apply_one_read assemble_read = create_one_assemble_read(one_read, reads_mem);
        if (!hc_assemble_finalize_region(assembler, assemble_read)) {
            destory_one_assemble_read(assemble_read);
            continue;
        }
        adjust_assemble_read_mem(assemble_read, reads_mem);
        CDL_APPEND(assembler->reads, assemble_read);
        assembler->reads_count++;
    }
    hc_apply_main(assembler);
    AssembleRegion* assemble_region = new (mem->allocate(sizeof(AssembleRegion))) AssembleRegion(region, mem);
    pSimpleInterval interval = assemble_region->paddedSpan;
    std::pmr::vector<pReadRecord> reads{mem};
    std::pmr::vector<pHaplotype> haplotypes{mem};
    if (!extract_reads(assembler, reads, mem)) {
        // std::cerr << "No reads." << std::endl;
    }
    if (!extract_haplotypes(assembler, haplotypes, interval, mem)) {
        // std::cerr << "No haplotypes" << std::endl;
    }
    AssembleResult* ret = create_one_assemble_result(mem);
    pSimpleInterval refLoc = SimpleInterval::create(region->tid, assembler->ref_start + 1, assembler->ref_end, mem);
    ret->set_reads(reads);
    ret->set_haplotypes(haplotypes);
    ret->set_padded_ref(const_cast<uint8_t*>(assembler->ref));
    ret->set_padded_refloc(refLoc);
    ret->set_region(assemble_region);

    return ret;
}

bool AssembleEngine::extract_reads(p_hc_apply assembler, std::pmr::vector<pReadRecord>& reads, pMemoryPool mem)
{
    p_hc_apply_one_read assemble_read;
    CDL_FOREACH(assembler->reads, assemble_read)
    {
        assemble_read->read.core.mpos += (assemble_read->read.core.flag & BAM_FPAIRED) ? 1 : 0;
        auto read = ReadRecord::create(mem, nullptr, assemble_read);
        reads.emplace_back(read);
    }
    return !reads.empty();
}

bool AssembleEngine::extract_haplotypes(p_hc_apply assembler, std::pmr::vector<pHaplotype>& haplotypes, pSimpleInterval region,
                                        pMemoryPool mem)
{
    p_hc_shared_lib_list_head path = assembler->read_graph.dijkstra_path_finder.result;
    p_hc_shared_lib_list_item node;
    CDL_FOREACH(path->nodes, node)
    {
        OneNodeData bestHaplotype = static_cast<OneNodeData>(node->data);
        int64_t wrt_ref = assembler->read_thread_graph->origin_ref_start <= HC_APPLY_REFERENCE_PADDING
                              ? assembler->read_thread_graph->origin_ref_start - 1
                              : HC_APPLY_REFERENCE_PADDING;
        pHaplotype hap = Haplotype::create(mem);
        // several classes. allVariationEvents 放上层吧！
        hap->init_haplotype(build_bases_from_apply(bestHaplotype, mem), bestHaplotype->is_reference);
        hap->set_cigar(build_cigar_from_apply(bestHaplotype, mem));
        hap->set_score(bestHaplotype->score);
        hap->set_alignment_start_hap_wrt_ref(wrt_ref);
        hap->set_interval(region);
        hap->set_kmer_size(bestHaplotype->is_reference ? 0 : bestHaplotype->kmer_size);
        haplotypes.emplace_back(hap);
    }
    return !haplotypes.empty();
}

pCigar AssembleEngine::build_cigar_from_apply(OneNodeData apply_node, pMemoryPool mem)
{
    using NodeMem = p_hc_assemble_dijkstra_out_path_storage;
    NodeMem node_storage = apply_node->st;
    pCigar ret = new ALLOC_FLEXIBLE_IN_POOL(mem, Cigar, node_storage->cigar_len, uint32_t) Cigar{};
    ret->num = node_storage->cigar_len;
    memcpy(ret->data, node_storage->cigar, sizeof(uint32_t) * ret->num);
    return ret;
}

pBases AssembleEngine::build_bases_from_apply(OneNodeData apply_node, pMemoryPool mem)
{
    static const uint32_t k_sw_padding_len = 10;
    std::string basesStr((char*)apply_node->st->seq + k_sw_padding_len, apply_node->st->seq_len - k_sw_padding_len * 2);
    char* bases = basesStr.data();
    return BaseUtils::bases_create(bases, mem);
}

void AssembleEngine::init_assemble_argument(const AssembleArgument* assemble_argument)
{
    const auto& graph_argument = assemble_argument->read_threading_argument;
    // const auto& sw_argument = assemble_argument->sw_argument;

    auto kmer_size = graph_argument.kmer.size();
    apply_arguments = (ApplyArgument*)malloc(sizeof(ApplyArgument));
    apply_arguments->kmers = (struct ApplyArgument::Kmer*)malloc(sizeof(struct ApplyArgument::Kmer) + sizeof(uint32_t) * kmer_size);
    apply_arguments->kmers->nkmer = kmer_size;
    for (size_t i = 0; i < kmer_size; i++) {
        apply_arguments->kmers->kmer[i] = graph_argument.kmer[i];
    }
    apply_arguments->debug.debug_print = false;
    apply_arguments->debug.debugGraphTransformations = false;
    apply_arguments->generateSeqGraph = true;
    apply_arguments->removePathsNotConnectedToRef = true;
    apply_arguments->pruneBeforeCycleCounting = true;
    apply_arguments->recover_dangling_branches = true;
    apply_arguments->debug.debug_print = assemble_argument->debugAssembly;
    apply_arguments->soft_clip_low_quality_ends = assemble_argument->softClipLowQualityEnds;
    apply_arguments->min_dangling_branch_length = graph_argument.minDanglingBranchLength;
    apply_arguments->max_unpruned_variants = graph_argument.maxUnprunedVariants;
    apply_arguments->recover_all_dangling_branches = graph_argument.recoverAllDanglingBranches;
    apply_arguments->allow_non_unique_kmer_in_ref = graph_argument.allowNonUniqueKmersInRef;
    apply_arguments->num_pruning_samples = graph_argument.numPruningSamples;
    apply_arguments->max_num_haplotypes_inpopulation = graph_argument.maxNumHaplotypesInPopulation;
    apply_arguments->min_pruning = graph_argument.minPruneFactor;
    apply_arguments->kmer_length_for_read_error_correction = graph_argument.kmerLengthForReadErrorCorrection;
    apply_arguments->min_observations_for_kmer_tobesolid = graph_argument.minObservationsForKmerToBeSolid;
    apply_arguments->min_matchingbases_to_dangling_end_recovery = graph_argument.minMatchingBasesToDanglingEndRecovery;
    apply_arguments->initial_error_rate_for_pruning = graph_argument.initialErrorRateForPruning;
    apply_arguments->pruning_logodds_threshold = graph_argument.pruningLogOddsThreshold;
    apply_arguments->pruning_seeding_logodds_threshold = graph_argument.pruningSeedingLogOddsThreshold;
    apply_arguments->pileup_error_correction_logodds = graph_argument.pileupErrorCorrectionLogOdds;
}

void AssembleEngine::finit_assemble_argument()
{
    free(apply_arguments->kmers);
    free(apply_arguments);
}
AssembleResult* AssembleEngine::create_one_assemble_result(pMemoryPool resource) { return AssembleResult::create(resource); }
void AssembleEngine::destory_one_assemble_result(AssembleResult* result) { result->~AssembleResult(); }
p_hc_apply_one_read AssembleEngine::create_one_assemble_read(bam1_t* one_read, pAssembleReadsBuffer reads_mem)
{
    if (!one_read) {
        return NULL;
    }

    p_hc_apply_one_read run_read = nullptr;
    uint32_t want_size = ((sizeof(hc_apply_one_read) + one_read->core.l_qseq + 7) & (~7U)) + k_assemble_cigar_string_mem +
                         sizeof(uint32_t) * k_assemble_cigar_cache_mem * 2;
    want_size = (want_size + 7) & (~7U);
    uint32_t this_read_used = 0;
    uint8_t* this_read_mem = reads_mem->buffer_ + reads_mem->used_;
    if (reads_mem->used_ + want_size < reads_mem->capacity_) {
        run_read = (p_hc_apply_one_read)(this_read_mem);
        this_read_used += sizeof(hc_apply_one_read);
        run_read->mempolicy = 1;
        run_read->seq = (uint8_t*)(this_read_mem + this_read_used);
        this_read_used += one_read->core.l_qseq;
        this_read_used = (this_read_used + 7) & (~7U);
        run_read->cigar_cache = (uint32_t*)(this_read_mem + this_read_used);
        this_read_used += (k_assemble_cigar_cache_mem * sizeof(uint32_t));
        run_read->cigar_cache_back = (uint32_t*)(this_read_mem + this_read_used);
        this_read_used += (k_assemble_cigar_cache_mem * sizeof(uint32_t));
        run_read->cigar_str = (char*)(this_read_mem + this_read_used);
        this_read_used += k_assemble_cigar_string_mem;
    }
    else {
        run_read = (p_hc_apply_one_read)calloc(1, sizeof(hc_apply_one_read));
        run_read->mempolicy = 0;
        run_read->cigar_str = (char*)malloc(sizeof(char) * k_assemble_cigar_string_mem);
        run_read->cigar_cache = (uint32_t*)malloc(sizeof(uint32_t) * k_assemble_cigar_cache_mem);
        run_read->cigar_cache_back = (uint32_t*)malloc(sizeof(uint32_t) * k_assemble_cigar_cache_mem);
        run_read->seq = (uint8_t*)malloc(sizeof(uint8_t) * one_read->core.l_qseq);
    }
    run_read->seq_cache = run_read->seq;
    // initialize
    run_read->read = *one_read;
    run_read->read_data = one_read->data;
    run_read->pos_start = one_read->core.pos + 1;  // transfer to 1base.
    run_read->pos_end = bam_endpos(one_read);
    run_read->cigar_len = one_read->core.n_cigar;
    run_read->read_len = (uint32_t)one_read->core.l_qseq;
    run_read->ref_len = 0;
    run_read->insert_size = one_read->core.isize;
    run_read->cigar = bam_get_cigar(one_read);
    run_read->cigar_back = bam_get_cigar(one_read);
    run_read->qual = bam_get_qual(one_read);
    run_read->prev = NULL;
    run_read->next = NULL;
    run_read->node = {};

    // Transfer seq to 8byte one base from 4byte.
    uint8_t* seq = bam_get_seq(one_read);
    for (uint32_t i = 0; i < run_read->read_len; ++i) {
        run_read->seq[i] = seq_nt16_str[bam_seqi(seq, i)];
    }

    return run_read;

    for (uint32_t i = 0; i < one_read->core.n_cigar; ++i) {
        char op = bam_cigar_opchr(run_read->cigar[i]);
        int len = bam_cigar_oplen(run_read->cigar[i]);
        if (i == 0) {
            sprintf(run_read->cigar_str, "%d%c", len, op);
            continue;
        }
        sprintf(run_read->cigar_str + strlen(run_read->cigar_str), "%d%c", len, op);
    }
}

void AssembleEngine::destory_one_assemble_read(p_hc_apply_one_read read)
{
    if (__glibc_unlikely(read && !read->mempolicy)) {
        if (read->cigar_str) free(read->cigar_str);
        if (read->cigar_cache) free(read->cigar_cache);
        if (read->cigar_cache_back) free(read->cigar_cache_back);
        if (read->seq_cache) free(read->seq_cache);
        free(read);
    }
}

void AssembleEngine::adjust_assemble_read_mem(p_hc_apply_one_read read, pAssembleReadsBuffer reads_mem)
{
    if (!read->mempolicy) return;
    // Use memmove maybe better?
    uint32_t real_size = ((sizeof(hc_apply_one_read) + read->read.core.l_qseq + 7) & (~7U)) + sizeof(uint32_t) * read->cigar_len;
    reads_mem->used_ = (reads_mem->used_ + real_size + 7) & (~7U);
    if (read->cigar == read->cigar_cache_back) {
        memmove(read->cigar_cache, read->cigar_cache_back, sizeof(uint32_t) * read->cigar_len);
        read->cigar = read->cigar_cache;
    }
    read->cigar_cache_back = NULL;
    read->cigar_str = NULL;
    read->cigar_back = NULL;
}