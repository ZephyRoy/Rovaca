/**
 * @file hc_assemble_read_threading_graph.c
 * @author Cao Ce (caoce@genomics.cn)
 * @brief Read Thread Graph
 * @version 0.1
 * @date 2021-10-22
 *
 * @copyright Copyright (c) 2021 ROVACA SDK
 *
 */

#include "debug.h"
#include "graph.h"
#include "hc_apply.h"
#include "hc_assemble.h"
#include "hc_func.h"
#include "mem_pool_auto_enlarge.h"

static inline void hc_assemble_graph_reset_thread_graph(p_hc_apply apply);
static void hc_assemble_graph_thread_sequence_read(p_hc_assemble_graph_read read, p_hc_apply apply);
static int hc_assemble_graph_read_find_seq_kmer_start(p_hc_assemble_graph_read read, p_hc_apply apply);
static p_assemble_vertex hc_assemble_graph_get_or_create_kmer_vertex(p_hc_assemble_graph_read read, int startPos, p_hc_apply apply);
static void hc_assemble_graph_increase_counts_in_matched_kmers(p_hc_assemble_graph_read read, p_assemble_vertex vertex,
                                                               uint8_t* originalKmer, int offset);
static inline p_assemble_vertex hc_assemble_graph_get_kmer_vertex(p_assemble_vertex prevVertex, int kmerStart,
                                                                  p_hc_assemble_graph_read read, p_hc_apply apply);

/**
 * @brief 在指定Kmer下检测Ref是否有dup
 *
 * @param apply         apply
 * @param kmer_size     kmer size
 *
 * @return int 1:成功, 0:失败, -1:异常
 */
int hc_assemble_graph_determine_ref_non_unique_kmers(p_hc_apply apply)
{
    int kmer_size = apply->read_graph.kmer;
    const uint8_t* ref_haplotype = apply->ref_haplotype;  // ref_haplotype
    int i, stop_position = apply->ref_hap_len - kmer_size;
    int res = 0;
    p_assemble_graph_hash_table head = hc_assemble_utils_check_hash_work(apply->read_graph.hash_buffer_1);

    for (i = 0; i <= stop_position; i++) {
        unsigned int key;

        key = assemble_graph_hash_get_key(head, (void*)ref_haplotype, kmer_size);
        if (assemble_graph_hash_search_with_key(head, (void*)ref_haplotype, kmer_size, key)) {
            goto err;
        }
        assemble_graph_hash_insert_with_key(head, (void*)ref_haplotype, kmer_size, (void*)ref_haplotype, key);
        ref_haplotype += 1;
    }

    res = 1;
err:
    assemble_graph_hash_reset(head);
    return res;
}

/**
 * @brief 重置graph
 *
 * @param apply apply
 */
static inline void hc_assemble_graph_reset_thread_graph(p_hc_apply apply)
{
    ring_mempool_auto_enlarge_reset(apply->read_graph.item_buffer);
    assemble_graph_reset(apply->read_thread_graph);
}

/**
 * @brief 检测Read的Unique Kmer
 *
 * @param read      read
 * @param apply     apply
 */
void hc_assemble_graph_read_determine_non_unique_kmers(p_hc_assemble_graph_read read, p_hc_apply apply)
{
    int kmer_size = apply->read_graph.kmer;
    int i, stop_position = read->dup_stop - kmer_size;
    p_assemble_graph_hash_node hash_item, read_hash_item;
    p_assemble_graph_hash_table read_hash = hc_assemble_utils_check_hash_work(apply->read_graph.hash_buffer_1);
    uint32_t kmer_non_unique_sum = 0;

    for (i = 0; i <= stop_position; i++) {
        uint32_t key = assemble_graph_hash_get_key(apply->read_thread_graph->non_unique_kmers, (void*)(read->dup_seq + i), kmer_size);

        // read内非独立kmer
        read_hash_item = assemble_graph_hash_search_with_key(read_hash, (void*)(read->dup_seq + i), kmer_size, key);
        if (!read_hash_item) {
            assemble_graph_hash_insert_with_key(read_hash, (void*)(read->dup_seq + i), kmer_size, (void*)(read->dup_seq + i), key);
        }
        else {
            // 总非独立kmer
            hash_item =
                assemble_graph_hash_search_with_key(apply->read_thread_graph->non_unique_kmers, (void*)(read->dup_seq + i), kmer_size, key);
            if (!hash_item) {
                kmer_non_unique_sum += 1;
                assemble_graph_hash_insert_with_key(apply->read_thread_graph->non_unique_kmers, (void*)(read->dup_seq + i), kmer_size,
                                                    (void*)(read->dup_seq + i), key);
            }
        }
    }
    apply->read_graph.kmer_non_unique_sum += kmer_non_unique_sum;
    assemble_graph_hash_reset(read_hash);
}

/**
 * @brief Thread sequence seqForKmers through the current graph, updating the graph as appropriate
 *
 * @param apply apply
 */
int hc_assemble_graph_thread_sequence(p_hc_apply apply)
{
    p_hc_assemble_graph_read graph_item;
    CDL_FOREACH(apply->read_graph.graph_reads, graph_item) { hc_assemble_graph_thread_sequence_read(graph_item, apply); }
    return 1;
}

/**
 * Find vertex and its position in seqForKmers where we should start assembling seqForKmers
 *
 * @param seqForKmers the sequence we want to thread into the graph
 * @return the position of the starting vertex in seqForKmer, or -1 if it cannot find one
 */
static int hc_assemble_graph_read_find_seq_kmer_start(p_hc_assemble_graph_read read, p_hc_apply apply)
{
    int i = 0, stop = read->len - apply->read_graph.kmer;
    p_assemble_graph_hash_node hash_item;

    if (__glibc_unlikely(read->is_ref)) {
        return 0;
    }

    for (; i < stop; i++) {
        uint32_t key =
            assemble_graph_hash_get_key(apply->read_thread_graph->non_unique_kmers, (void*)(read->seq + i), apply->read_graph.kmer);
        hash_item = assemble_graph_hash_search_with_key(apply->read_thread_graph->non_unique_kmers, (void*)(read->seq + i),
                                                        apply->read_graph.kmer, key);
        if (!hash_item) {
            return i;
        }
    }

    return -1;
}

/**
 * Get the vertex for the kmer in sequence starting at start
 *
 * @param sequence the sequence
 * @param start    the position of the kmer start
 * @return a non-null vertex
 */
static p_assemble_vertex hc_assemble_graph_get_or_create_kmer_vertex(p_hc_assemble_graph_read read, int startPos, p_hc_apply apply)
{
    p_assemble_vertex vertex;
    p_hc_assemble_graph_kmer kmer_item;

    vertex = assemble_graph_find_vertex(apply->read_thread_graph, (void*)(read->seq + startPos), apply->read_graph.kmer);
    if (vertex) {
        return vertex;
    }

    kmer_item = hc_assemble_utils_fill_kmer_by_read(read, startPos, apply);
    vertex = hc_assemble_utils_add_graph_vertex_by_kmer(kmer_item, apply);

    return vertex;
}

/**
 * Get the suffix byte of this DeBruijnVertex
 *
 * The suffix byte is simply the last byte of the kmer sequence, so if this is holding sequence ACT
 * hc_assemble_graph_get_suffix would return T
 *
 * @return a byte
 */
uint8_t hc_assemble_graph_get_suffix(p_assemble_vertex vertex)
{
    uint8_t* data = vertex->data;

    return data[vertex->data_len - 1];
}

static void hc_assemble_graph_increase_counts_in_matched_kmers(p_hc_assemble_graph_read read, p_assemble_vertex vertex,
                                                               uint8_t* originalKmer, int offset)
{
    p_assemble_edge edge;
    p_assemble_vertex prev_vertex;
    uint8_t suffix, seqBase;

    if (offset <= -1) {
        return;
    }

    CDL_FOREACH2(vertex->from_edges, edge, link.to.next)
    {
        prev_vertex = edge->link.from.vertex;
        suffix = hc_assemble_graph_get_suffix(prev_vertex);
        seqBase = originalKmer[offset];
        if (suffix == seqBase && vertex->in_degree == 1) {
            edge->weight += 1;
            hc_assemble_graph_increase_counts_in_matched_kmers(read, prev_vertex, originalKmer, offset - 1);
        }
    }
}

static p_assemble_vertex hc_assemble_graph_extend_chain_by_one_ref(p_assemble_vertex prevVertex, int kmerStart,
                                                                   p_hc_assemble_graph_read read, p_hc_apply apply)
{
    p_hc_assemble_graph_kmer kmer_item;
    p_assemble_vertex nextVertex;

    kmer_item = hc_assemble_utils_fill_kmer_by_read(read, kmerStart, apply);
    nextVertex = hc_assemble_utils_add_graph_vertex_by_kmer(kmer_item, apply);

    // either use our merge vertex, or create a new one in the chain
    if (!assemble_graph_vertex_add_edge_weight_between_vertex(prevVertex, nextVertex, 1) &&
        !assemble_graph_vertex_add_edge_to_vertex(apply->read_thread_graph, prevVertex, nextVertex, 1)) {
        statics_print_error("No Mem");
        exit(-1);
    }

    return nextVertex;
}

/**
 * Get the unique vertex for kmer, or null if not possible.
 *
 * @param allowRefSource if true, we will allow kmer to match the reference source vertex
 * @return a vertex for kmer, or null (either because it doesn't exist or is non-unique for graphs that have such a distinction)
 */
static inline p_assemble_vertex hc_assemble_graph_get_kmer_vertex(p_assemble_vertex prevVertex, int kmerStart,
                                                                  p_hc_assemble_graph_read read, p_hc_apply apply)
{
    (void)prevVertex;
    p_assemble_vertex vertex;

    if (__glibc_unlikely(memcmp(apply->read_graph.ref_source->data, (void*)(read->seq + kmerStart), apply->read_graph.kmer) == 0)) {
        return NULL;
    }

    vertex = assemble_graph_find_vertex(apply->read_thread_graph, (void*)(read->seq + kmerStart), apply->read_graph.kmer);
    if (!vertex || vertex->is_dup) {
        return NULL;
    }
    return vertex;
}

#if 0
/**
 * @brief 链接尾部至source vertex
 * 
 * @param prev_vertex       prev vertex
 * @param next_vertex       next vertex
 * @param graph             graph
 * 
 * @return p_assemble_vertex  next vertex
 */
static inline p_assemble_vertex hc_assemble_graph_extend_source_vertex(p_assemble_vertex prev_vertex, p_assemble_vertex next_vertex, p_assemble_graph graph)
{
    p_assemble_edge from_edges, outgoing_edges = NULL;
    p_assemble_vertex working_vertex;
    uint32_t added_edge_weight = 0;
    int32_t i;

    if (__glibc_unlikely(next_vertex->in_degree)) {
        return next_vertex;
    }

    // 首先链接两个vertex
    if (!assemble_graph_vertex_add_edge_weight_between_vertex(prev_vertex, next_vertex, 1) &&
        !(outgoing_edges = assemble_graph_vertex_add_edge_to_vertex(graph, prev_vertex, next_vertex, 1))) {
        statics_print_error("No Mem");
        exit(-1);
    }
    if (outgoing_edges) {
        outgoing_edges->is_ref = 0;
    }

    // 获取新增edge权
    CDL_FOREACH2(next_vertex->out_edges, outgoing_edges, link.from.next) {
        added_edge_weight += outgoing_edges->weight;
    }
    
    for (i = next_vertex->data_len - 1, working_vertex = next_vertex; i > 0; i --) {
        int32_t target_data_start = next_vertex->data_len - i;
        p_assemble_vertex target = NULL;
        int get_target = 0;

        CDL_FOREACH2(working_vertex->from_edges, from_edges, link.to.next) {
            target = from_edges->link.from.vertex;
            char* target_seq = target->data;

            if (memcmp(target_seq + target_data_start, next_vertex->data, i) == 0) {
                from_edges->weight += added_edge_weight;
                working_vertex = target;
                get_target = 1;
                break;
            }
        }
        if (!get_target) {
            return next_vertex;
        }
    }
    return next_vertex;
}
#endif
/**
 * Workhorse routine of the assembler.  Given a sequence whose last vertex is anchored in the graph, extend
 * the graph one bp according to the bases in sequence.
 *
 * @param prevVertex a non-null vertex where sequence was last anchored in the graph
 * @param sequence   the sequence we're threading through the graph
 * @param kmerStart  the start of the current kmer in graph we'd like to add
 * @param count      the number of observations of this kmer in graph (can be > 1 for GGA)
 * @param isRef      is this the reference sequence?
 * @return a non-null vertex connecting prevVertex to in the graph based on sequence
 */
static p_assemble_vertex hc_assemble_graph_extend_chain_by_one(p_assemble_vertex prevVertex, int kmerStart, p_hc_assemble_graph_read read,
                                                               p_hc_apply apply)
{
    p_assemble_edge outgoingEdges = prevVertex->out_edges, outgoingEdge;
    p_hc_assemble_graph_kmer kmer_item;
    p_assemble_vertex mergeVertex = NULL, nextVertex;
    int nextPos = kmerStart + apply->read_graph.kmer - 1;

    if (__glibc_unlikely(read->is_ref)) {
        return hc_assemble_graph_extend_chain_by_one_ref(prevVertex, kmerStart, read, apply);
    }

    CDL_FOREACH2(outgoingEdges, outgoingEdge, link.from.next)
    {
        p_assemble_vertex target = outgoingEdge->link.to.vertex;
        if (hc_assemble_graph_get_suffix(target) == read->seq[nextPos]) {
            // we've got a match in the chain, so simply increase the count of the edge by 1 and continue
            outgoingEdge->weight += 1;
            return target;
        }
    }

    // none of our outgoing edges had our unique suffix base, so we check for an opportunity to merge back in
    mergeVertex = hc_assemble_graph_get_kmer_vertex(prevVertex, kmerStart, read, apply);
    if (mergeVertex == NULL) {
        assemble_graph_to_working_read_thread_graph(apply->read_thread_graph);
        kmer_item = hc_assemble_utils_fill_kmer_by_read(read, kmerStart, apply);
        nextVertex = hc_assemble_utils_add_graph_vertex_by_kmer(kmer_item, apply);
        assemble_graph_to_normal_read_thread_graph(apply->read_thread_graph);
    }
    else {
        nextVertex = mergeVertex;
        // if (__glibc_unlikely(!nextVertex->from_edges)) {
        //     return hc_assemble_graph_extend_source_vertex(prevVertex, nextVertex, apply->read_thread_graph);
        // }
    }

    // either use our merge vertex, or create a new one in the chain
    if (!assemble_graph_vertex_add_edge_weight_between_vertex(prevVertex, nextVertex, 1) &&
        !(outgoingEdge = assemble_graph_vertex_add_edge_to_vertex(apply->read_thread_graph, prevVertex, nextVertex, 1))) {
        statics_print_error("No Mem");
        exit(-1);
    }
    if (outgoingEdge) {
        outgoingEdge->is_ref = 0;
    }

    return nextVertex;
}

/**
 * @brief read thread sequence
 *
 * @param read      read
 * @param apply     apply
 */
static void hc_assemble_graph_thread_sequence_read(p_hc_assemble_graph_read read, p_hc_apply apply)
{
    p_assemble_vertex startingVertex, vertex;
    int startPos, i;

    startPos = hc_assemble_graph_read_find_seq_kmer_start(read, apply);
    if (startPos == -1) {
        return;
    }
    if (read->is_ref) {
        startingVertex = hc_assemble_utils_add_graph_vertex_by_kmer(hc_assemble_utils_fill_kmer_by_read(read, startPos, apply), apply);
    }
    else {
        startingVertex = hc_assemble_graph_get_or_create_kmer_vertex(read, startPos, apply);
    }
    // increase the counts of all edges incoming into the starting vertex supported by going back in sequence
    hc_assemble_graph_increase_counts_in_matched_kmers(read, startingVertex, startingVertex->data, apply->read_graph.kmer - 2);

    // keep track of information about the reference source
    if (read->is_ref) {
        apply->read_thread_graph->origin_ref_start = apply->span_start;
        memcpy(apply->read_thread_graph->origin_ref, startingVertex->data, startingVertex->data_len);
        apply->read_thread_graph->origin_ref_len = startingVertex->data_len;
        apply->read_graph.ref_source = startingVertex;
    }

    // loop over all of the bases in sequence, extending the graph by one base at each point, as appropriate
    vertex = startingVertex;
    for (i = startPos + 1; i <= read->len - apply->read_graph.kmer; i++) {
        vertex = hc_assemble_graph_extend_chain_by_one(vertex, i, read, apply);
        if (read->is_ref && vertex && apply->read_thread_graph->origin_ref_len < ASSEMBLE_GRAPH_MAX_SEQ_LEN) {
            apply->read_thread_graph->origin_ref[apply->read_thread_graph->origin_ref_len] = hc_assemble_graph_get_suffix(vertex);
            apply->read_thread_graph->origin_ref_len += 1;
            apply->read_graph.ref_end = vertex;
        }
        if (__glibc_unlikely(!vertex)) {
            vertex = startingVertex;
        }
    }
}