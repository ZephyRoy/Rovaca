/**
 * @file hc_assemble_graph_wrapper.c
 * @author Cao Ce (caoce@genomics.cn)
 * @brief HC Assemble Graph
 * @version 0.1
 * @date 2021-05-26
 *
 * @copyright Copyright (c) 2021 ROVACA SDK
 *
 */

#include "assemble_interface.h"
#include "debug.h"
#include "hc_apply.h"
#include "hc_assemble.h"
#include "hc_assemble_debug_print.h"
#include "hc_assemble_graph.h"
#include "hc_func.h"
#include "htslib/sam.h"
#include "mem_pool_auto_enlarge.h"
#include "rbtree.h"
#include "rbtree_shared.h"

#define HC_ASSEMBLE_GRAPH_LOW_QUAL_MUL (4)  // 一个固定值，不依赖参数。

static inline int hc_assemble_graph_usable_for_assembly(char base, uint8_t qual);
static void hc_assemble_graph_add_ref_item(p_hc_assemble_graph_read graph_item, p_hc_apply apply);
static void hc_assemble_graph_add_read(p_hc_apply_one_read read, p_hc_apply apply);
static void hc_assemble_graph_build_preprocess_reads(p_hc_apply apply);
static int hc_assemble_build_graph_necessarily(p_hc_apply apply);
static int hc_assemble_graph_finder_node_down(p_hc_apply apply, p_assemble_vertex source, p_assemble_vertex sink);
static int hc_assemble_graph_finder_node_up(p_hc_apply apply, p_assemble_vertex source, p_assemble_vertex sink);
int hc_assemble_graph_remove_paths_not_connected_to_ref(p_hc_apply apply);
static inline int hc_assemble_graph_is_low_quality_graph(p_hc_apply apply);

/**
 * @brief ReadThreadingAssembler::createGraph
 * @param apply apply
 * @return int 1:succeed, 0:failed
 * @note 无内存直接退出
 */
int hc_assemble_graph_build(p_hc_apply apply)
{
    p_hc_apply_one_read read;
    struct rb_node* node;
    int res;

    // Happens in cases where the assembled region is just too small.
    if (apply->ref_hap_len < (uint32_t)apply->read_graph.kmer) {
        return 0;
    }

    assemble_graph_to_normal_read_thread_graph(apply->read_thread_graph);

    // Check non-unique kmer in ref.
    if (__glibc_likely(apply->read_graph.kmer != HC_ASSEMBLE_MAX_KMER) && !apply_arguments->allow_non_unique_kmer_in_ref &&
        !hc_assemble_graph_determine_ref_non_unique_kmers(apply)) {
        return 0;
    }

    // Add reference and reads to graph.
    hc_assemble_graph_add_ref_item(&apply->read_graph.ref, apply);
    for (node = rb_first(&apply->read_sorted); node; node = rb_next(node)) {
        read = rb_entry(node, hc_apply_one_read, node);
        hc_assemble_graph_add_read(read, apply);
    }
    if (apply_arguments->debug.debug_print) {
        printf("\nReads for assembling graph:\treads count: %d\n", apply->read_graph.graph_reads_count);
        p_hc_assemble_graph_read graph_item;
        CDL_FOREACH(apply->read_graph.graph_reads, graph_item)
        {
            printf("%s\t%d-%d\n", graph_item->name, graph_item->start, graph_item->end);
        }
    }
    if (!hc_assemble_build_graph_necessarily(apply)) {
        res = 0;
        goto end;
    }

    if (apply_arguments->debug.debugGraphTransformations) {
        hc_assemble_graph_layout(apply->read_thread_graph, "0.0.raw_readthreading_graph.dot");
    }

    // 修整图
    assemble_graph_to_working_read_thread_graph(apply->read_thread_graph);
    if (__glibc_likely(apply_arguments->recover_dangling_branches)) {
        hc_assemble_base_graph_recover_dangling_tails(apply);
        hc_assemble_base_graph_recover_dangling_heads(apply);
    }

    if (__glibc_likely(apply_arguments->removePathsNotConnectedToRef)) {
        hc_assemble_graph_remove_paths_not_connected_to_ref(apply);
    }

    if (__glibc_unlikely(apply_arguments->debug.debugGraphTransformations)) {
        hc_assemble_graph_layout(apply->read_thread_graph, "0.2.cleaned_readthreading_graph.dot");
    }

    if (__glibc_likely(apply_arguments->generateSeqGraph)) {
        hc_assemble_seq_graph_to_sequence_graph(apply->read_thread_graph);
        if (apply_arguments->debug.debugGraphTransformations) {
            hc_assemble_graph_layout(apply->read_thread_graph, "0.3.initial_seqgraph.dot");
        }

        {
            // the very first thing we need to do is zip up the graph, or pruneGraph will be too aggressive
            hc_assemble_seq_graph_zip_linear_chains(apply);
            // now go through and prune the graph, removing vertices no longer connected to the reference chain
            hc_assemble_seq_graph_remove_singleton_orphan_vertices(apply);
            hc_assemble_seq_graph_remove_vertices_not_connected_to_ref_regardless_of_edge_direction(apply);
            hc_assemble_seq_graph_simplify_graph(apply);

            // The graph has degenerated in some way, so the reference source and/or sink cannot be id'd.  Can
            // happen in cases where for example the reference somehow manages to acquire a cycle, or
            // where the entire assembly collapses back into the reference sequence.
            if (!hc_assemble_utils_get_reference_source_vertex_with_head(apply->read_thread_graph->vertics_hash->node_list,
                                                                         apply->read_thread_graph->vertics_hash->node_list,
                                                                         apply->read_thread_graph) ||
                !hc_assemble_utils_get_reference_sink_vertex_with_head(apply->read_thread_graph->vertics_hash->node_list,
                                                                       apply->read_thread_graph->vertics_hash->node_list,
                                                                       apply->read_thread_graph)) {
                res = 0;
                goto end;
            }
            hc_assemble_seq_graph_simplify_graph(apply);
        }  // cleanupSeqGraph.
    }

    hc_assemble_seq_finder_find_best_paths(apply);
    res = apply->read_graph.dijkstra_path_finder.result ? 1 : 0;

end:
    // reset
    assemble_graph_reset(apply->read_thread_graph);
    ring_mempool_auto_enlarge_reset(apply->read_graph.kmer_buffer);
    ring_mempool_auto_enlarge_reset(apply->read_graph.item_buffer);
    hc_assemble_thread_mem_reset(apply->read_graph.cache_pool);
    apply->read_graph.graph_reads = NULL;
    apply->read_graph.graph_reads_count = 0;

    return res;
}

/**
 * @brief   read的seq及qual是否合适做assemble
 *
 * @param base  seq
 * @param qual  qual
 *
 * @return int 1:True, 0:False
 */
static inline int hc_assemble_graph_usable_for_assembly(char base, uint8_t qual)
{
    return base != 'N' && qual >= HC_ASSEMBLE_GRAPH_MIN_QUAL;
}

/**
 * @brief 将ref以graph read的形式添加进去
 *
 * @param graph_item    graph read
 * @param apply         apply
 */
static void hc_assemble_graph_add_ref_item(p_hc_assemble_graph_read graph_item, p_hc_apply apply)
{
    graph_item->seq = apply->ref_haplotype;
    graph_item->qual = NULL;
    graph_item->len = apply->ref_hap_len;

    graph_item->dup_seq = graph_item->seq;
    graph_item->dup_stop = graph_item->len;

    graph_item->is_ref = 1;
    graph_item->is_seq = 0;
    graph_item->is_dup = 0;
    graph_item->start = 0;
    graph_item->end = apply->ref_hap_len;

    graph_item->name = "ref";
    graph_item->next = NULL;
    graph_item->prev = NULL;

    CDL_APPEND(apply->read_graph.graph_reads, graph_item);
    apply->read_graph.graph_reads_count++;
}

/**
 * @brief
 *
 * @param read
 * @param apply
 * @param kmerSize
 */
static void hc_assemble_graph_add_read(p_hc_apply_one_read read, p_hc_apply apply)
{
    int lastGood = -1;
    uint32_t end;
    uint8_t* sequence = read->seq;
    uint8_t* qualities = read->qual;
    p_hc_assemble_graph_read graph_item;

    for (end = 0; end <= read->read_len; end++) {
        if (end == read->read_len || !hc_assemble_graph_usable_for_assembly(sequence[end], qualities[end])) {
            // the first good base is at lastGood, can be -1 if last base was bad
            int start = lastGood;
            // the stop base is end - 1 (if we're not at the end of the sequence)
            int len = end - start;

            if (start != -1 && len >= apply->read_graph.kmer) {
                graph_item = ring_mempool_auto_enlarge_malloc(apply->read_graph.item_buffer);
                graph_item->start = start;
                graph_item->end = end;
                graph_item->dup_seq = sequence;
                graph_item->dup_stop = end;
                graph_item->seq = sequence + start;
                graph_item->qual = qualities + start;
                graph_item->len = end - start;
                graph_item->is_ref = 0;
                graph_item->is_seq = 1;
                graph_item->is_dup = 0;
                graph_item->name = bam_get_qname(&read->read);
                CDL_APPEND(apply->read_graph.graph_reads, graph_item);
                apply->read_graph.graph_reads_count++;
            }

            lastGood = -1;  // reset the last good base
        }
        else if (lastGood == -1) {
            lastGood = end;  // we're at a good base, the last good one is us
        }
    }
}
static int hc_assemble_build_graph_necessarily(p_hc_apply apply)
{
    hc_assemble_graph_build_preprocess_reads(apply);
    hc_assemble_graph_thread_sequence(apply);

    // It's important to prune before recovering dangling ends so that we don't waste time recovering bad ends.
    // It's also important to prune before checking for cycles so that sequencing errors don't create false cycles
    // and unnecessarily abort assembly.
    if (__glibc_likely(apply_arguments->pruneBeforeCycleCounting)) {
        hc_assemble_chain_pruner_prune_low_weight_chains(apply);
    }
    if (hc_assemble_base_graph_cycle_finder(apply) != 0 ||
        (apply->read_graph.kmer != HC_ASSEMBLE_MAX_KMER && hc_assemble_graph_is_low_quality_graph(apply))) {
        return 0;
    }
    return 1;
}

/**
 * Does the graph not have enough complexity?  We define low complexity as a situation where the number
 * of non-unique kmers is more than 20% of the total number of kmers.
 *
 * @return true if the graph has low complexity, false otherwise
 */
static inline int hc_assemble_graph_is_low_quality_graph(p_hc_apply apply)
{
    return apply->read_graph.kmer_non_unique_sum * HC_ASSEMBLE_GRAPH_LOW_QUAL_MUL > apply->read_thread_graph->vertices_sum;
}

static void hc_assemble_graph_build_preprocess_reads(p_hc_apply apply)
{
    p_hc_assemble_graph_read graph_item;

    apply->read_graph.kmer_non_unique_sum = 0;

    hc_assemble_graph_read_determine_non_unique_kmers(&apply->read_graph.ref, apply);
    CDL_FOREACH(apply->read_graph.graph_reads, graph_item) { hc_assemble_graph_read_determine_non_unique_kmers(graph_item, apply); }
}

/**
 * @brief 搜索出一条路径后执行的操作
 *
 * @param last_vertex   上一个点
 * @param iter_st       遍历储存
 *
 * @return int 0
 */
static int hc_assemble_graph_finder_down_travel_one_path_done(p_assemble_vertex last_vertex, p_assemble_graph_iter iter_st)
{
    (void)(last_vertex);
    p_hc_shared_lib_list_item path_head = iter_st->run_path;
    p_hc_shared_lib_list_item path_el;

    CDL_FOREACH(path_head, path_el)
    {
        p_assemble_vertex search_v = path_el->data;
        search_v->color |= ASSEMBLE_GRAPH_VERTEX_COLOR_RED;
    }
    return 0;
}

/**
 * @brief graph finder down
 *
 * @param apply apply
 *
 * @return int 0:succeed, 1:cycle
 */
static int hc_assemble_graph_finder_node_down(p_hc_apply apply, p_assemble_vertex source, p_assemble_vertex sink)
{
    int ret;
    p_assemble_graph_iter travel_config = &apply->read_thread_graph->graph_iter;

    travel_config->input_layer.from = source;
    travel_config->input_layer.to = sink;
    travel_config->input_layer.apply = apply;
    travel_config->cycle_detected = hc_assemble_utils_default_graph_func_cycle;
    travel_config->one_path_done = hc_assemble_graph_finder_down_travel_one_path_done;
    travel_config->all_path_done = NULL;

    ret = hc_assemble_base_graph_iter_dfs_non_recursion(travel_config);

    return ret;
}

/**
 * @brief 搜索出一条路径后执行的操作
 *
 * @param last_vertex   上一个点
 * @param iter_st       遍历储存
 *
 * @return int 0
 */
static int hc_assemble_graph_finder_down_travel_one_path_done_up(p_assemble_vertex last_vertex, p_assemble_graph_iter iter_st)
{
    (void)last_vertex;
    p_hc_shared_lib_list_item path_head = iter_st->run_path;
    p_hc_shared_lib_list_item path_el;

    CDL_FOREACH(path_head, path_el)
    {
        p_assemble_vertex search_v = path_el->data;
        search_v->color |= ASSEMBLE_GRAPH_VERTEX_COLOR_YELLOW;
    }

    return 0;
}

/**
 * @brief graph finder up
 *
 * @param apply apply
 *
 * @return int 0:succeed, 1:cycle
 */
static int hc_assemble_graph_finder_node_up(p_hc_apply apply, p_assemble_vertex source, p_assemble_vertex sink)
{
    int ret;
    p_assemble_graph_iter travel_config = &apply->read_thread_graph->graph_iter;

    travel_config->input_layer.from = source;
    travel_config->input_layer.to = sink;
    travel_config->input_layer.apply = apply;
    travel_config->cycle_detected = hc_assemble_utils_default_graph_func_cycle;
    travel_config->one_path_done = hc_assemble_graph_finder_down_travel_one_path_done_up;
    travel_config->all_path_done = NULL;

    ret = hc_assemble_base_graph_iter_dfs_non_recursion_up(travel_config);

    return ret;
}

/**
 * Remove all vertices in the graph that aren't on a path from the reference source vertex to the reference sink vertex
 *
 * More aggressive reference pruning algorithm than removeVerticesNotConnectedToRefRegardlessOfEdgeDirection,
 * as it requires vertices to not only be connected by a series of directed edges but also prunes away
 * paths that do not also meet eventually with the reference sink vertex
 */
int hc_assemble_graph_remove_paths_not_connected_to_ref(p_hc_apply apply)
{
    p_assemble_graph_hash_node vertex_node, vertex_el1, vertex_el2;
    p_assemble_vertex ref_source, ref_sink;
    p_assemble_vertex vertex;

    if (apply->read_thread_graph->vertices_sum == 1) {
        return 0;
    }

    // cache
    CDL_FOREACH(apply->read_thread_graph->vertics_hash->node_list, vertex_node)
    {
        vertex = vertex_node->value;
        vertex->color = ASSEMBLE_GRAPH_VERTEX_COLOR_NONE;
    }

    // append
    ref_source = apply->read_graph.ref_source;
    ref_sink = apply->read_graph.ref_end;
    if (ref_source == NULL || ref_sink == NULL) {
        return 0;
    }

    if (ref_source->out_edges) {
        hc_assemble_graph_finder_node_down(apply, ref_source, ref_sink);
    }
    if (ref_sink->from_edges) {
        hc_assemble_graph_finder_node_up(apply, ref_sink, ref_source);
    }

    // we want to remove anything that's not in both the sink and source sets
    CDL_FOREACH_SAFE(apply->read_thread_graph->vertics_hash->node_list, vertex_node, vertex_el1, vertex_el2)
    {
        vertex = vertex_node->value;

        if (vertex->color == (ASSEMBLE_GRAPH_VERTEX_COLOR_RED | ASSEMBLE_GRAPH_VERTEX_COLOR_YELLOW)) {
            vertex->color = ASSEMBLE_GRAPH_VERTEX_COLOR_NONE;
        }
        else {
            assemble_graph_remove_vertex(apply->read_thread_graph, vertex_node->value);
        }
    }

    // simple sanity checks that this algorithm is working.
    vertex_el1 = hc_assemble_utils_get_reference_source_vertex_with_head(
        apply->read_thread_graph->vertics_hash->node_list, apply->read_thread_graph->vertics_hash->node_list, apply->read_thread_graph);
    vertex_el2 = hc_assemble_utils_get_reference_sink_vertex_with_head(
        apply->read_thread_graph->vertics_hash->node_list, apply->read_thread_graph->vertics_hash->node_list, apply->read_thread_graph);
    if (vertex_el1 == NULL || vertex_el2 == NULL) {
        return 0;
    }
    vertex_el1 = hc_assemble_utils_get_reference_source_vertex_with_head(vertex_el1, apply->read_thread_graph->vertics_hash->node_list,
                                                                         apply->read_thread_graph);
    vertex_el2 = hc_assemble_utils_get_reference_sink_vertex_with_head(vertex_el2, apply->read_thread_graph->vertics_hash->node_list,
                                                                       apply->read_thread_graph);
    if (vertex_el1 == NULL || vertex_el2 == NULL) {
        return 0;
    }
    return 1;
}