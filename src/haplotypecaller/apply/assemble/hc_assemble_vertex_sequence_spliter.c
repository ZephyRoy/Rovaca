/**
 * @file hc_assemble_vertex_sequence_spliter.c
 * @author Cao Ce (caoce@genomics.cn)
 * @brief 将图中的一组中间节点拆分为其共享的前缀和后缀值
 * @version 0.1
 * @date 2021-11-26
 *
 * @copyright Copyright (c) 2021 ROVACA SDK
 *
 */

#include <limits.h>

#include "graph.h"
#include "hc_apply.h"
#include "hc_assemble.h"
#include "hc_func.h"
#include "hc_marco.h"
#include "list.h"

#define MIN_COMMON_SEQUENCE_TO_MERGE_SOURCE_SINK_VERTICES (10)

static int hc_assemble_vertex_sequence_spliter_can_merge(p_assemble_vertex vertex, p_hc_shared_lib_list_head incoming_vertices);

static p_assemble_vertex hc_assemble_vertex_sequence_spliter_common_suffix_one(p_hc_shared_lib_list_head middle_vertices, p_hc_apply apply);
static int hc_assemble_vertex_sequence_spliter_would_eliminate_ref_source(p_assemble_vertex common_suffix,
                                                                          p_hc_shared_lib_list_head to_splits, p_assemble_graph graph);
static int hc_assemble_vertex_sequence_spliter_all_vertices_are_the_common_suffix(p_assemble_vertex common_suffix,
                                                                                  p_hc_shared_lib_list_head to_splits);
static inline int hc_assemble_vertex_sequence_spliter_vertex_out_degree_contians_vertex(p_assemble_vertex search, p_assemble_vertex target);
static int hc_assemble_vertex_sequence_spliter_safe_to_split(p_assemble_vertex bot, p_hc_shared_lib_list_head toMerge);
static p_assemble_vertex hc_assemble_vertex_sequence_spliter_common_suffix(p_assemble_vertex vertex, p_hc_shared_lib_list_head to_split,
                                                                           p_hc_apply apply);
static p_assemble_vertex hc_assemble_vertex_sequence_spliter_without_suffix(p_assemble_vertex sequence, p_assemble_vertex suffix,
                                                                            p_hc_apply apply);
static int hc_assemble_vertex_sequence_spliter_common_suffix_splitter_split(p_hc_apply apply, p_assemble_vertex vertex);

static inline int hc_assemble_vertex_sequence_spliter_meets_min_mergable_sequence_for_either_prefix_or_suffix(
    int minCommonSequence, p_hc_assemble_read_graph graph);
static inline int hc_assemble_vertex_sequence_spliter_meets_min_mergable_sequence_for_prefix(int minCommonSequence,
                                                                                             p_hc_assemble_read_graph graph);
static inline int hc_assemble_vertex_sequence_spliter_meets_min_mergable_sequence_for_suffix(int minCommonSequence,
                                                                                             p_hc_assemble_read_graph graph);
static void hc_assemble_vertex_sequence_spliter_common_prefix_and_suffix_of_vertices(p_hc_shared_lib_list_head middle_vertices,
                                                                                     p_hc_apply apply);
static int hc_assemble_vertex_sequence_spliter_common_maximum_prefix_length(p_hc_shared_lib_list_head listOfBytes, int minLength);
static int hc_assemble_vertex_sequence_spliter_common_maximum_suffix_length(p_hc_shared_lib_list_head listOfBytes, int minLength);
static p_assemble_vertex hc_assemble_vertex_sequence_spliter_without_prefix_and_suffix(p_assemble_graph split_graph,
                                                                                       p_assemble_vertex sequence, p_assemble_vertex prefix,
                                                                                       p_assemble_vertex suffix);
static int hc_assemble_vertex_sequence_spliter_split_and_update(p_hc_shared_lib_list_head to_splits, p_assemble_vertex top,
                                                                p_assemble_vertex bottom, p_hc_apply apply);
static void hc_assemble_vertex_sequence_spliter_split(p_hc_shared_lib_list_head to_splits, p_hc_apply apply);
static void hc_assemble_vertex_sequence_spliter_update_graph(p_hc_apply apply, p_hc_shared_lib_list_head to_splits, p_assemble_vertex top,
                                                             p_assemble_vertex bot);
static void hc_assemble_vertex_sequence_spliter_add_prefix_node_and_edges(p_hc_apply apply, p_assemble_vertex top);
static void hc_assemble_vertex_sequence_spliter_add_suffix_node_and_edges(p_hc_apply apply, p_assemble_vertex bot);
static void hc_assemble_vertex_sequence_spliter_add_edges_from_top_node(p_hc_apply apply, p_assemble_vertex top_for_connect,
                                                                        p_assemble_vertex bot_for_connect);
static void hc_assemble_vertex_sequence_spliter_add_edges_to_bottom_node(p_hc_apply apply, p_assemble_vertex bot_for_connect);

/**
 * Merge diamond configurations:
 *
 * Performance the transformation:
 *
 * { A -> x + S_i + y -> Z }
 *
 * goes to:
 *
 * { A -> x -> S_i -> y -> Z }
 *
 * for all nodes that match this configuration.
 */
int hc_assemble_vertex_sequence_spliter_merge_diamonds_try_to_transform(p_assemble_vertex top, p_hc_apply apply)
{
    p_assemble_vertex mt, mi, bottom;
    p_assemble_edge edge, el;
    uint32_t middles = 0;
    p_hc_shared_lib_list_head middles_head;
    int ret;

    // we can only merge if there's at least two middle nodes
    if (top->out_degree <= 1) {
        return 0;
    }

    middles_head = hc_assemble_utils_check_list_work(apply->read_graph.list_buffer_2);
    CDL_FOREACH2(top->out_edges, edge, link.from.next)
    {
        mi = edge->link.to.vertex;
        hc_shared_list_insert(middles_head, mi);
        middles += 1;
    }

    bottom = NULL;
    CDL_FOREACH2(top->out_edges, edge, link.from.next)
    {
        mi = edge->link.to.vertex;
        // all nodes must have at least 1 connection
        if (mi->out_degree < 1) {
            hc_shared_list_reset(middles_head);
            return 0;
        }

        // can only have 1 incoming node, the root vertex
        if (mi->in_degree != 1) {
            hc_shared_list_reset(middles_head);
            return 0;
        }

        // make sure that all outgoing vertices of mi go only to the bottom node
        CDL_FOREACH2(mi->out_edges, el, link.from.next)
        {
            mt = el->link.to.vertex;
            if (bottom == NULL) {
                bottom = mt;
            }
            else if (bottom != mt) {
                hc_shared_list_reset(middles_head);
                return 0;
            }
        }
    }

    // bottom has some connections coming in from other nodes, don't allow
    if (bottom->in_degree != (int)middles) {
        hc_shared_list_reset(middles_head);
        return 0;
    }

    hc_assemble_vertex_sequence_spliter_common_prefix_and_suffix_of_vertices(middles_head, apply);
    // actually do the merging, returning true if at least 1 base was successfully split
    ret = hc_assemble_vertex_sequence_spliter_meets_min_mergable_sequence_for_either_prefix_or_suffix(1, &apply->read_graph) &&
          hc_assemble_vertex_sequence_spliter_split_and_update(middles_head, top, bottom, apply);

    hc_shared_list_reset(middles_head);
    hc_shared_pair_list_reset(apply->read_graph.vertex_spliter.middles);
    hc_shared_list_reset(apply->read_graph.vertex_spliter.edges_to_remove);
    assemble_graph_reset(apply->read_graph.vertex_spliter.seq_graph);
    apply->read_graph.vertex_spliter.prefix_vertex = apply->read_graph.vertex_spliter.suffix_vertex =
        apply->read_graph.vertex_spliter.prefix_vertex_main_graph = apply->read_graph.vertex_spliter.suffix_vertex_main_graph = NULL;

    return ret;
}

/**
 * Merge tail configurations:
 *
 * Performs the transformation:
 *
 * { A -> x + S_i + y }
 *
 * goes to:
 *
 * { A -> x -> S_i -> y }
 *
 * for all nodes that match this configuration.
 *
 * Differs from the diamond transform in that no bottom node is required
 */
int hc_assemble_vertex_sequence_spliter_merge_tails_try_to_transform(p_assemble_vertex top, p_hc_apply apply)
{
    p_assemble_vertex vertex;
    p_assemble_edge edge;
    p_hc_shared_lib_list_head tails;
    int ret;

    if (top->out_degree <= 1) {
        return 0;
    }
    tails = hc_assemble_utils_check_list_work(apply->read_graph.list_buffer_2);

    CDL_FOREACH2(top->out_edges, edge, link.from.next)
    {
        vertex = edge->link.to.vertex;
        if (vertex->out_degree == 0 || vertex->in_degree > 1) {
            hc_shared_list_reset(tails);
            return 0;
        }
        hc_shared_list_insert(tails, vertex);
    }

    hc_assemble_vertex_sequence_spliter_common_prefix_and_suffix_of_vertices(tails, apply);

    ret = hc_assemble_vertex_sequence_spliter_meets_min_mergable_sequence_for_suffix(MIN_COMMON_SEQUENCE_TO_MERGE_SOURCE_SINK_VERTICES,
                                                                                     &apply->read_graph) &&
          hc_assemble_vertex_sequence_spliter_split_and_update(tails, top, NULL, apply);

    hc_shared_list_reset(tails);
    hc_shared_pair_list_reset(apply->read_graph.vertex_spliter.middles);
    assemble_graph_reset(apply->read_graph.vertex_spliter.seq_graph);
    apply->read_graph.vertex_spliter.prefix_vertex = apply->read_graph.vertex_spliter.suffix_vertex =
        apply->read_graph.vertex_spliter.prefix_vertex_main_graph = apply->read_graph.vertex_spliter.suffix_vertex_main_graph = NULL;

    return ret;
}

/**
 * Performs the transformation:
 *
 * { x + S_i + y -> Z }
 *
 * goes to:
 *
 * { x -> S_i -> y -> Z }
 *
 * for all nodes that match this configuration.
 *
 * Differs from the diamond transform in that no top node is required
 */
int hc_assemble_vertex_sequence_spliter_split_common_suffices_try_to_transform(p_assemble_vertex bottom, p_hc_apply apply)
{
    p_assemble_graph_hash_table alreadySplit = apply->read_graph.hash_buffer_1;
    uint32_t key = bottom->key;

    if (assemble_graph_hash_search_pt_with_key(alreadySplit, bottom, key)) {
        return 0;
    }
    else {
        assemble_graph_hash_insert_with_key(alreadySplit, &bottom, sizeof(void*), bottom, key);
        return hc_assemble_vertex_sequence_spliter_common_suffix_splitter_split(apply, bottom);
    }
}

/**
 * Merge headless configurations:
 *
 * Performs the transformation:
 *
 * { x + S_i -> y -> Z }
 *
 * goes to:
 *
 * { x -> S_i -> y + Z }
 *
 * for all nodes that match this configuration.
 */
int hc_assemble_vertex_sequence_spliter_split_merge_common_suffices(p_assemble_vertex vertex, p_hc_apply apply)
{
    p_hc_shared_lib_list_head prevs = hc_assemble_utils_check_list_work(apply->read_graph.list_buffer_1);
    p_hc_shared_lib_list_head edges_to_remove = hc_assemble_utils_check_list_work(apply->read_thread_graph->alt_path);
    p_hc_shared_lib_list_item item;
    p_assemble_edge edge, cached_edge;
    p_assemble_vertex cached_vertex, prev_seq_v, new_v;

    CDL_FOREACH2(vertex->from_edges, edge, link.to.next)
    {
        cached_vertex = edge->link.from.vertex;
        if (!cached_vertex) {
            continue;
        }
        hc_shared_list_insert(prevs, cached_vertex);
    }

    if (!hc_assemble_vertex_sequence_spliter_can_merge(vertex, prevs)) {
        return 0;
    }

    prev_seq_v = prevs->nodes->data;
    new_v = hc_assemble_utils_merge_seq_vertex_from_two_seq_vertex(prev_seq_v, vertex, apply);
    assemble_graph_add_vertex(apply->read_thread_graph, new_v);

    CDL_FOREACH(prevs->nodes, item)
    {
        cached_vertex = item->data;

        CDL_FOREACH2(cached_vertex->from_edges, edge, link.to.next)
        {
            hc_shared_list_insert(edges_to_remove, edge);
            if (__glibc_unlikely(!edge->link.from.vertex)) {
                continue;
            }
            cached_edge = assemble_graph_vertex_get_edge_between_vertex(edge->link.from.vertex, new_v);
            if (__glibc_likely(!cached_edge)) {
                cached_edge =
                    assemble_graph_vertex_add_edge_to_vertex(apply->read_thread_graph, edge->link.from.vertex, new_v, edge->weight);
                cached_edge->is_ref = edge->is_ref;
            }
        }
    }

    CDL_FOREACH2(vertex->out_edges, edge, link.from.next)
    {
        if (__glibc_unlikely(!edge->link.to.vertex)) {
            continue;
        }
        cached_edge = assemble_graph_vertex_get_edge_between_vertex(new_v, edge->link.to.vertex);
        if (__glibc_likely(!cached_edge)) {
            cached_edge = assemble_graph_vertex_add_edge_to_vertex(apply->read_thread_graph, new_v, edge->link.to.vertex, edge->weight);
            cached_edge->is_ref = edge->is_ref;
        }
    }

    CDL_FOREACH(prevs->nodes, item)
    {
        cached_vertex = item->data;
        assemble_graph_remove_vertex(apply->read_thread_graph, cached_vertex);
    }
    hc_shared_list_reset(prevs);
    assemble_graph_remove_vertex(apply->read_thread_graph, vertex);
    CDL_FOREACH(edges_to_remove->nodes, item)
    {
        edge = item->data;
        assemble_graph_remove_edge(apply->read_thread_graph, edge);
    }
    hc_shared_list_reset(edges_to_remove);

    return 1;
}

/**
 * Can we safely merge the incoming vertices of v
 *
 * @param graph the graph containing v and incoming_vertices
 * @param v the vertex we want to merge into
 * @param incoming_vertices the incoming vertices of v
 * @return true if we can safely merge incoming_vertices
 */
static int hc_assemble_vertex_sequence_spliter_can_merge(p_assemble_vertex vertex, p_hc_shared_lib_list_head incoming_vertices)
{
    p_hc_shared_lib_list_item item;
    p_assemble_vertex first, prev;

    if (!incoming_vertices->nodes) {
        return 0;
    }

    first = incoming_vertices->nodes->data;
    CDL_FOREACH(incoming_vertices->nodes, item)
    {
        prev = item->data;
        // cannot merge if our sequence isn't the same as the first sequence
        if (!(prev->data_len == first->data_len && (memcmp(prev->data, first->data, prev->data_len) == 0))) {
            return 0;
        }
        // prev -> vertex must be the only edge from prev
        if (prev->out_degree != 1) {
            return 0;
        }
        // don't allow cyles
        if (prev->out_edges && prev->out_edges->link.to.vertex != vertex) {
            return 0;
        }
        // cannot merge when any of the incoming nodes are sources
        if (prev->in_degree == 0) {
            return 0;
        }
    }

    return 1;
}

/**
 * Return the longest suffix of bases shared among all provided vertices
 *
 * For example, if the vertices have sequences AC, CC, and ATC, this would return
 * a single C.  However, for ACC and TCC this would return CC.  And for AC and TG this
 * would return null;
 *
 * @param middle_vertices a non-empty set of vertices
 * @return a single vertex that contains the common suffix of all middle vertices
 */
static p_assemble_vertex hc_assemble_vertex_sequence_spliter_common_suffix_one(p_hc_shared_lib_list_head middle_vertices, p_hc_apply apply)
{
    p_hc_assemble_graph_kmer suffix_kmer;
    p_hc_shared_lib_list_item el;
    p_assemble_vertex vertex;
    uint32_t min_len = UINT_MAX;
    int suffix_len;

    CDL_FOREACH(middle_vertices->nodes, el)
    {
        vertex = el->data;
        min_len = assemble_min(min_len, vertex->data_len);
    }
    suffix_len = hc_assemble_vertex_sequence_spliter_common_maximum_suffix_length(middle_vertices, min_len);

    vertex = middle_vertices->nodes->data;
    suffix_kmer = ring_mempool_auto_enlarge_malloc(apply->read_graph.kmer_buffer);
    memset(suffix_kmer, 0, sizeof(hc_assemble_graph_kmer));

    suffix_kmer->seq = ((uint8_t*)(vertex->data)) + vertex->data_len - suffix_len;
    suffix_kmer->len = suffix_len;
    suffix_kmer->item = ((p_hc_assemble_graph_kmer)(vertex->storage))->item;

    vertex = assemble_graph_vertex_create(apply->read_graph.vertex_spliter.seq_graph, suffix_kmer->seq, suffix_kmer->len, suffix_kmer);
    // assemble_graph_add_vertex(apply->read_graph.vertex_spliter.seq_graph, apply->read_graph.vertex_spliter.suffix_vertex);
    return vertex;
}

/**
 * Would factoring out this suffix result in elimating the reference source vertex?
 * @param graph the graph
 * @param common_suffix the common suffix of all to_splits
 * @param to_splits the list of vertices we're are trying to split
 * @return true if to_split contains the reference source and this ref source has all and only the bases of common_suffix
 */
static int hc_assemble_vertex_sequence_spliter_would_eliminate_ref_source(p_assemble_vertex common_suffix,
                                                                          p_hc_shared_lib_list_head to_splits, p_assemble_graph graph)
{
    p_hc_shared_lib_list_item to_split;
    p_assemble_vertex vertex;

    CDL_FOREACH(to_splits->nodes, to_split)
    {
        vertex = to_split->data;

        if (hc_assemble_base_graph_is_ref_source(vertex, graph)) {
            return vertex->data_len == common_suffix->data_len;
        }
    }
    return 0;
}

/**
 * Would all vertices that we'd split just result in the common suffix?
 *
 * That is, suppose we have prefix nodes ABC and ABC.  After splitting all of the vertices would
 * just be ABC again, and we'd enter into an infinite loop.
 *
 * @param common_suffix the common suffix of all vertices in to_splits
 * @param to_splits the collection of vertices we want to split
 * @return true if all of the vertices are equal to the common suffix
 */
static int hc_assemble_vertex_sequence_spliter_all_vertices_are_the_common_suffix(p_assemble_vertex common_suffix,
                                                                                  p_hc_shared_lib_list_head to_splits)
{
    p_hc_shared_lib_list_item to_split;
    p_assemble_vertex vertex;

    CDL_FOREACH(to_splits->nodes, to_split)
    {
        vertex = to_split->data;

        if (vertex->data_len != common_suffix->data_len) {
            return 0;
        }
    }
    return 1;
}

/**
 * @brief 搜索out vertex是否包括target
 *
 * @param search out vertex
 * @param target target
 *
 * @return int 0:failed, 1:succeed
 */
static inline int hc_assemble_vertex_sequence_spliter_vertex_out_degree_contians_vertex(p_assemble_vertex search, p_assemble_vertex target)
{
    p_assemble_edge edge;

    CDL_FOREACH2(search->out_edges, edge, link.from.next)
    {
        if (edge->link.to.vertex == target) {
            return 1;
        }
    }

    return 0;
}

/**
 * Can we safely split up the vertices in toMerge?
 *
 * @param graph a graph
 * @param bot a vertex whose incoming vertices we want to split
 * @param toMerge the set of vertices we'd be splitting up
 * @return true if we can safely split up toMerge
 */
static int hc_assemble_vertex_sequence_spliter_safe_to_split(p_assemble_vertex bot, p_hc_shared_lib_list_head toMerge)
{
    p_hc_shared_lib_list_item item;
    p_assemble_vertex m;

    CDL_FOREACH(toMerge->nodes, item)
    {
        m = item->data;

        // m == bot => don't allow self cycles in the graph
        if (m == bot || m->out_degree != 1 || !hc_assemble_vertex_sequence_spliter_vertex_out_degree_contians_vertex(m, bot)) {
            return 0;
        }
        // forbid cycles from bottom -> mid
        if (hc_assemble_vertex_sequence_spliter_vertex_out_degree_contians_vertex(bot, m)) {
            return 0;
        }
    }
    return 1;
}

/**
 * @brief Split a collection of middle nodes in a graph into their shared prefix and suffix values
 *
 * @param vertex        vertex
 * @param to_split      to split
 * @param apply         apply
 *
 * @return p_assemble_vertex    suffix vertex
 */
static p_assemble_vertex hc_assemble_vertex_sequence_spliter_common_suffix(p_assemble_vertex vertex, p_hc_shared_lib_list_head to_split,
                                                                           p_hc_apply apply)
{
    p_assemble_vertex suffix_v_template;

    if (!hc_assemble_vertex_sequence_spliter_safe_to_split(vertex, to_split)) {
        return NULL;
    }

    suffix_v_template = hc_assemble_vertex_sequence_spliter_common_suffix_one(to_split, apply);
    if (suffix_v_template->data_len == 0) {
        return NULL;
    }
    else if (hc_assemble_vertex_sequence_spliter_would_eliminate_ref_source(suffix_v_template, to_split, apply->read_thread_graph)) {
        return NULL;
    }
    else if (hc_assemble_vertex_sequence_spliter_all_vertices_are_the_common_suffix(suffix_v_template, to_split)) {
        return NULL;
    }
    else {
        return suffix_v_template;
    }
}

/**
 * Return a new SeqVertex derived from this one but not including the suffix bases
 *
 * @param suffix the suffix bases to remove from this vertex
 * @return a newly allocated SeqVertex with appropriate prefix, or null if suffix removes all bases from this node
 */
static p_assemble_vertex hc_assemble_vertex_sequence_spliter_without_suffix(p_assemble_vertex sequence, p_assemble_vertex suffix,
                                                                            p_hc_apply apply)
{
    int prefix_size = sequence->data_len - suffix->data_len;

    if (prefix_size <= 0) {
        return NULL;
    }
    return hc_assemble_utils_create_seq_vertex_from_data(sequence->data, prefix_size, apply);
}

/**
 * Simple single-function interface to split and then update a graph
 *
 * @param graph the graph containing the vertices in toMerge
 * @param v The bottom node whose incoming vertices we'd like to split
 * @return true if some useful splitting was done, false otherwise
 */
static int hc_assemble_vertex_sequence_spliter_common_suffix_splitter_split(p_hc_apply apply, p_assemble_vertex vertex)
{
    p_assemble_graph graph = apply->read_thread_graph;
    p_hc_shared_lib_list_head to_split, edges_to_remove;
    p_hc_shared_lib_list_item item;
    p_assemble_vertex suffix_v_template, mid, suffix_v, prefix_v, incoming_target;
    p_assemble_edge edge, out, loop_edge, added_edge;

    // Can only split at least 2 vertices
    if (vertex->in_degree < 2) {
        return 0;
    }

    to_split = hc_assemble_utils_check_list_work(apply->read_graph.list_buffer_1);
    CDL_FOREACH2(vertex->from_edges, edge, link.to.next)
    {
        if (__glibc_likely(edge->link.from.vertex != NULL)) {
            hc_shared_list_insert(to_split, edge->link.from.vertex);
        }
    }

    suffix_v_template = hc_assemble_vertex_sequence_spliter_common_suffix(vertex, to_split, apply);
    if (suffix_v_template == NULL) {
        hc_shared_list_reset(to_split);
        return 0;
    }

    edges_to_remove = hc_assemble_utils_check_list_work(apply->read_thread_graph->alt_path);
    CDL_FOREACH(to_split->nodes, item)
    {
        mid = item->data;

        // create my own copy of the suffix
        suffix_v = hc_assemble_utils_create_seq_vertex_from_data(suffix_v_template->data, suffix_v_template->data_len, apply);
        assemble_graph_add_vertex(graph, suffix_v);
        prefix_v = hc_assemble_vertex_sequence_spliter_without_suffix(mid, suffix_v, apply);
        out = mid->out_edges;

        if (prefix_v == NULL) {
            // this node is entirely explained by suffix
            incoming_target = suffix_v;
        }
        else {
            incoming_target = prefix_v;
            assemble_graph_add_vertex(graph, prefix_v);
            added_edge = assemble_graph_vertex_add_edge_to_vertex(graph, prefix_v, suffix_v, 1);
            added_edge->is_ref = out->is_ref;
            hc_shared_list_insert(edges_to_remove, out);
        }

        added_edge = assemble_graph_vertex_add_edge_to_vertex(graph, suffix_v, out->link.to.vertex, 1);
        added_edge->is_ref = out->is_ref;

        CDL_FOREACH2(mid->from_edges, loop_edge, link.to.next)
        {
            added_edge = assemble_graph_vertex_add_edge_to_vertex(graph, loop_edge->link.from.vertex, incoming_target, loop_edge->weight);
            added_edge->is_ref = loop_edge->is_ref;
            hc_shared_list_insert(edges_to_remove, loop_edge);
        }
    }

    CDL_FOREACH(to_split->nodes, item) { assemble_graph_remove_vertex(graph, item->data); }
    CDL_FOREACH(edges_to_remove->nodes, item) { assemble_graph_remove_edge(graph, item->data); }
    hc_shared_list_reset(to_split);
    hc_shared_list_reset(edges_to_remove);

    return 1;
}

/**
 * Does either the common suffix or prefix have at least minCommonSequence bases in it?
 * @param minCommonSequence a minimum length of the common sequence, must be >= 0
 * @return true if either suffix or prefix length >= minCommonSequence
 */
static inline int hc_assemble_vertex_sequence_spliter_meets_min_mergable_sequence_for_either_prefix_or_suffix(
    int minCommonSequence, p_hc_assemble_read_graph graph)
{
    return hc_assemble_vertex_sequence_spliter_meets_min_mergable_sequence_for_prefix(minCommonSequence, graph) ||
           hc_assemble_vertex_sequence_spliter_meets_min_mergable_sequence_for_suffix(minCommonSequence, graph);
}

/**
 * Does the common prefix have at least minCommonSequence bases in it?
 * @param minCommonSequence a minimum length of the common sequence, must be >= 0
 * @return true if prefix length >= minCommonSequence
 */
static inline int hc_assemble_vertex_sequence_spliter_meets_min_mergable_sequence_for_prefix(int minCommonSequence,
                                                                                             p_hc_assemble_read_graph graph)
{
    return graph->vertex_spliter.prefix_vertex->data_len >= (uint32_t)minCommonSequence;
}

/**
 * Does the common suffix have at least minCommonSequence bases in it?
 * @param minCommonSequence a minimum length of the common sequence, must be >= 0
 * @return true if suffix length >= minCommonSequence
 */
static inline int hc_assemble_vertex_sequence_spliter_meets_min_mergable_sequence_for_suffix(int minCommonSequence,
                                                                                             p_hc_assemble_read_graph graph)
{
    return graph->vertex_spliter.suffix_vertex->data_len >= (uint32_t)minCommonSequence;
}

/**
 * Return the longest suffix of bases shared among all provided vertices
 *
 * For example, if the vertices have sequences AC, CC, and ATC, this would return
 * a single C.  However, for ACC and TCC this would return CC.  And for AC and TG this
 * would return null;
 *
 * @param middle_vertices a non-empty set of vertices
 * @return
 */
static void hc_assemble_vertex_sequence_spliter_common_prefix_and_suffix_of_vertices(p_hc_shared_lib_list_head middle_vertices,
                                                                                     p_hc_apply apply)
{
    p_hc_shared_lib_list_head kmers = hc_assemble_utils_check_list_work(apply->read_graph.list_buffer_1);
    p_hc_shared_lib_list_item el;
    p_assemble_vertex vertex;
    p_hc_assemble_graph_kmer pref_kmer, suffix_kmer;
    uint32_t min_len = UINT_MAX, prefixLen, suffix_len;

    CDL_FOREACH(middle_vertices->nodes, el)
    {
        vertex = el->data;

        hc_shared_list_insert(kmers, vertex->data);
        min_len = assemble_min(min_len, vertex->data_len);
    }

    prefixLen = hc_assemble_vertex_sequence_spliter_common_maximum_prefix_length(kmers, min_len);
    suffix_len = hc_assemble_vertex_sequence_spliter_common_maximum_suffix_length(middle_vertices, min_len - prefixLen);

    // kmer
    vertex = middle_vertices->nodes->data;
    pref_kmer = ring_mempool_auto_enlarge_malloc(apply->read_graph.kmer_buffer);
    suffix_kmer = ring_mempool_auto_enlarge_malloc(apply->read_graph.kmer_buffer);
    memset(pref_kmer, 0, sizeof(hc_assemble_graph_kmer));
    memset(suffix_kmer, 0, sizeof(hc_assemble_graph_kmer));

    pref_kmer->seq = vertex->data;
    pref_kmer->len = prefixLen;
    pref_kmer->item = ((p_hc_assemble_graph_kmer)(vertex->storage))->item;
    suffix_kmer->seq = ((uint8_t*)(vertex->data)) + vertex->data_len - suffix_len;
    suffix_kmer->len = suffix_len;
    suffix_kmer->item = ((p_hc_assemble_graph_kmer)(vertex->storage))->item;

    // vertex
    apply->read_graph.vertex_spliter.prefix_vertex =
        assemble_graph_vertex_create(apply->read_graph.vertex_spliter.seq_graph, pref_kmer->seq, pref_kmer->len, pref_kmer);
    apply->read_graph.vertex_spliter.suffix_vertex =
        assemble_graph_vertex_create(apply->read_graph.vertex_spliter.seq_graph, suffix_kmer->seq, suffix_kmer->len, suffix_kmer);
    assemble_graph_add_vertex(apply->read_graph.vertex_spliter.seq_graph, apply->read_graph.vertex_spliter.prefix_vertex);
    assemble_graph_add_vertex(apply->read_graph.vertex_spliter.seq_graph, apply->read_graph.vertex_spliter.suffix_vertex);

    hc_shared_list_reset(kmers);
}

/**
 * Compute the maximum shared prefix length of list of bytes.
 *
 * @param listOfBytes a list of bytes with at least one element
 * @return the number of shared bytes common at the start of all bytes
 */
static int hc_assemble_vertex_sequence_spliter_common_maximum_prefix_length(p_hc_shared_lib_list_head listOfBytes, int minLength)
{
    p_hc_shared_lib_list_item item;
    int i = 0;
    uint8_t data;

    for (i = 0; i < minLength; i++) {
        item = listOfBytes->nodes;
        data = ((uint8_t*)(item->data))[i];
        CDL_FOREACH(listOfBytes->nodes, item)
        {
            if (item == listOfBytes->nodes) {
                continue;
            }
            if (data != ((uint8_t*)(item->data))[i]) {
                return i;
            }
        }
    }

    return minLength;
}

/**
 * Compute the maximum shared suffix length of list of bytes.
 *
 * @param listOfBytes a list of bytes with at least one element
 * @param minLength the min. length among all byte[] in listOfBytes
 * @return the number of shared bytes common at the end of all bytes
 */
static int hc_assemble_vertex_sequence_spliter_common_maximum_suffix_length(p_hc_shared_lib_list_head listOfBytes, int minLength)
{
    p_hc_shared_lib_list_item item;
    p_assemble_vertex vertex;
    uint8_t data;
    int suffix_len = 0;

    for (; suffix_len < minLength; suffix_len++) {
        item = listOfBytes->nodes;
        vertex = item->data;
        data = ((uint8_t*)(vertex->data))[vertex->data_len - 1 - suffix_len];
        CDL_FOREACH(listOfBytes->nodes, item)
        {
            vertex = item->data;

            if (item == listOfBytes->nodes) {
                continue;
            }
            if (data != ((uint8_t*)(vertex->data))[vertex->data_len - 1 - suffix_len]) {
                return suffix_len;
            }
        }
    }
    return minLength;
}

/**
 * Return a new SeqVertex derived from this one but not including prefix or suffix bases
 *
 * @param prefix the previx bases to remove
 * @param suffix the suffix bases to remove from this vertex
 * @return a newly allocated SeqVertex
 */
static p_assemble_vertex hc_assemble_vertex_sequence_spliter_without_prefix_and_suffix(p_assemble_graph split_graph,
                                                                                       p_assemble_vertex sequence, p_assemble_vertex prefix,
                                                                                       p_assemble_vertex suffix)
{
    int start = prefix->data_len;
    int length = sequence->data_len - suffix->data_len - prefix->data_len;
    int stop = start + length;
    p_assemble_vertex ret;

    if (length <= 0) {
        return NULL;
    }
    ret = assemble_graph_vertex_create(split_graph, (uint8_t*)(sequence->data) + start, stop - start, sequence->storage);
    return ret;
}

/**
 * Simple single-function interface to split and then update a graph
 *
 * @see #hc_assemble_vertex_sequence_spliter_update_graph(SeqVertex, SeqVertex) for a full description of top and bottom
 *
 * @param top the top vertex, may be null
 * @param bottom the bottom vertex, may be null
 * @return true if some useful splitting was done, false otherwise
 */
static int hc_assemble_vertex_sequence_spliter_split_and_update(p_hc_shared_lib_list_head to_splits, p_assemble_vertex top,
                                                                p_assemble_vertex bottom, p_hc_apply apply)
{
    hc_assemble_vertex_sequence_spliter_split(to_splits, apply);
    hc_assemble_vertex_sequence_spliter_update_graph(apply, to_splits, top, bottom);
    return 1;
}

/**
 * Actually do the splitting up of the vertices
 *
 * Must be called before calling hc_assemble_vertex_sequence_spliter_update_graph
 */
static void hc_assemble_vertex_sequence_spliter_split(p_hc_shared_lib_list_head to_splits, p_hc_apply apply)
{
    p_assemble_graph split_graph = apply->read_graph.vertex_spliter.seq_graph;
    p_hc_shared_lib_list_head new_middles = hc_assemble_utils_check_list_work(apply->read_thread_graph->alt_path);
    p_hc_shared_lib_list_item item;
    p_assemble_vertex mid, remaining;
    p_assemble_edge to_mid, from_mid, edge;

    CDL_FOREACH(to_splits->nodes, item)
    {
        mid = item->data;
        to_mid = mid->from_edges;
        from_mid = mid->out_edges;

        hc_shared_list_insert(apply->read_graph.vertex_spliter.edges_to_remove, to_mid);
        hc_shared_list_insert(apply->read_graph.vertex_spliter.edges_to_remove, from_mid);
        remaining = hc_assemble_vertex_sequence_spliter_without_prefix_and_suffix(
            split_graph, mid, apply->read_graph.vertex_spliter.prefix_vertex, apply->read_graph.vertex_spliter.suffix_vertex);
        if (remaining != NULL) {
            // there's some sequence prefix + seq + suffix, so add the node and make edges
            assemble_graph_add_vertex(split_graph, remaining);
            hc_shared_list_insert(new_middles, remaining);
            // update edge from top -> middle to be top -> without suffix
            edge = assemble_graph_vertex_add_edge_to_vertex(split_graph, apply->read_graph.vertex_spliter.prefix_vertex, remaining,
                                                            to_mid->weight);
            edge->is_ref = to_mid->is_ref;
            edge = assemble_graph_vertex_add_edge_to_vertex(split_graph, remaining, apply->read_graph.vertex_spliter.suffix_vertex,
                                                            from_mid->weight);
            edge->is_ref = from_mid->is_ref;
        }
        else {
            // prefix + suffix completely explain this node
            double weight = to_mid->weight + from_mid->weight;
            uint32_t is_ref = to_mid->is_ref || from_mid->is_ref;

            edge = assemble_graph_vertex_get_edge_between_vertex(apply->read_graph.vertex_spliter.prefix_vertex,
                                                                 apply->read_graph.vertex_spliter.suffix_vertex);
            if (!edge) {
                edge = assemble_graph_vertex_add_edge_to_vertex(split_graph, apply->read_graph.vertex_spliter.prefix_vertex,
                                                                apply->read_graph.vertex_spliter.suffix_vertex, weight);
                edge->is_ref = 0;
            }
            else {
                edge->weight += weight;
            }
            edge->is_ref = edge->is_ref || is_ref;
        }
    }
}

/**
 * Update graph outer, replacing the previous middle vertices that were split out with the new
 * graph structure of the split, linking this subgraph into the graph at top and bot (the
 * vertex connecting the middle nodes and the vertex outgoing of all middle node)
 *
 * @param top an optional top node that must have outgoing edges to all split vertices.  If null, this subgraph
 *            will be added without any incoming edges
 * @param bot an optional bottom node that must have incoming edges to all split vertices.  If null, this subgraph
 *            will be added without any outgoing edges to the rest of the graph
 */
static void hc_assemble_vertex_sequence_spliter_update_graph(p_hc_apply apply, p_hc_shared_lib_list_head to_splits, p_assemble_vertex top,
                                                             p_assemble_vertex bot)
{
    p_assemble_graph outer = apply->read_thread_graph;
    p_hc_shared_lib_list_item item;
    p_assemble_vertex vertex, graph_vertex;
    p_hc_shared_lib_list_head new_middles = apply->read_thread_graph->alt_path;
    p_assemble_vertex top_for_connect = NULL, bot_for_connect = NULL;
    p_hc_shared_lib_pair_list_head st_head = hc_assemble_utils_check_pair_list_work(apply->read_graph.vertex_spliter.middles);
    int has_prefix_suffix_edge, has_only_prefix_suffix_edges, need_prefix_node, need_suffix_node;

    CDL_FOREACH(to_splits->nodes, item)
    {
        vertex = item->data;
        assemble_graph_remove_vertex(outer, vertex);
    }
    CDL_FOREACH(apply->read_graph.vertex_spliter.edges_to_remove->nodes, item)
    {
        p_assemble_edge edge = item->data;
        assemble_graph_remove_edge(outer, edge);
    }

    CDL_FOREACH(new_middles->nodes, item)
    {
        vertex = item->data;

        graph_vertex = assemble_graph_vertex_create(outer, vertex->data, vertex->data_len, vertex->storage);
        assemble_graph_add_vertex(outer, graph_vertex);
        hc_shared_pair_list_insert(st_head, graph_vertex, vertex);
    }
    hc_shared_list_reset(new_middles);

    has_prefix_suffix_edge = assemble_graph_vertex_get_edge_between_vertex(apply->read_graph.vertex_spliter.prefix_vertex,
                                                                           apply->read_graph.vertex_spliter.suffix_vertex) != NULL;
    has_only_prefix_suffix_edges = has_prefix_suffix_edge && apply->read_graph.vertex_spliter.prefix_vertex->out_degree == 1;
    need_prefix_node = apply->read_graph.vertex_spliter.prefix_vertex->data_len || (top == NULL && !has_only_prefix_suffix_edges);
    need_suffix_node = apply->read_graph.vertex_spliter.suffix_vertex->data_len || (bot == NULL && !has_only_prefix_suffix_edges);

    // if prefix / suffix are needed, keep them
    if (need_prefix_node) {
        hc_assemble_vertex_sequence_spliter_add_prefix_node_and_edges(apply, top);
        top_for_connect = apply->read_graph.vertex_spliter.prefix_vertex_main_graph;
    }
    else {
        top_for_connect = top;
    }

    if (need_suffix_node) {
        hc_assemble_vertex_sequence_spliter_add_suffix_node_and_edges(apply, bot);
        bot_for_connect = apply->read_graph.vertex_spliter.suffix_vertex_main_graph;
    }
    else {
        bot_for_connect = bot;
    }

    if (top_for_connect != NULL) {
        hc_assemble_vertex_sequence_spliter_add_edges_from_top_node(apply, top_for_connect, bot_for_connect);
    }

    if (bot_for_connect != NULL) {
        hc_assemble_vertex_sequence_spliter_add_edges_to_bottom_node(apply, bot_for_connect);
    }
}

static void hc_assemble_vertex_sequence_spliter_add_prefix_node_and_edges(p_hc_apply apply, p_assemble_vertex top)
{
    p_assemble_vertex vertex, prefix_vertex = apply->read_graph.vertex_spliter.prefix_vertex;
    p_assemble_edge edge, el;

    vertex = assemble_graph_vertex_create(apply->read_thread_graph, prefix_vertex->data, prefix_vertex->data_len, prefix_vertex->storage);
    vertex = assemble_graph_add_vertex(apply->read_thread_graph, vertex);
    apply->read_graph.vertex_spliter.prefix_vertex_main_graph = vertex;

    if (__glibc_unlikely(top == NULL)) {
        return;
    }

    edge = assemble_graph_vertex_add_edge_to_vertex(apply->read_thread_graph, top, vertex, 1);
    if (__glibc_unlikely(!edge)) {
        return;
    }
    CDL_FOREACH2(prefix_vertex->out_edges, el, link.from.next)
    {
        if (el->is_ref) {
            edge->is_ref = 1;
            return;
        }
    }
    edge->is_ref = 0;
}

static void hc_assemble_vertex_sequence_spliter_add_suffix_node_and_edges(p_hc_apply apply, p_assemble_vertex bot)
{
    p_assemble_vertex vertex, suffix_vertex = apply->read_graph.vertex_spliter.suffix_vertex;
    p_assemble_edge edge, el;

    vertex = assemble_graph_vertex_create(apply->read_thread_graph, suffix_vertex->data, suffix_vertex->data_len, suffix_vertex->storage);
    vertex = assemble_graph_add_vertex(apply->read_thread_graph, vertex);
    apply->read_graph.vertex_spliter.suffix_vertex_main_graph = vertex;

    if (__glibc_unlikely(bot == NULL)) {
        return;
    }

    edge = assemble_graph_vertex_add_edge_to_vertex(apply->read_thread_graph, vertex, bot, 1);
    if (__glibc_unlikely(!edge)) {
        return;
    }
    CDL_FOREACH2(suffix_vertex->from_edges, el, link.to.next)
    {
        if (el->is_ref) {
            edge->is_ref = 1;
            return;
        }
    }
    edge->is_ref = 0;
}

static void hc_assemble_vertex_sequence_spliter_add_edges_from_top_node(p_hc_apply apply, p_assemble_vertex top_for_connect,
                                                                        p_assemble_vertex bot_for_connect)
{
    p_assemble_edge edge;
    p_assemble_vertex target;
    p_assemble_edge run_edge;

    CDL_FOREACH2(apply->read_graph.vertex_spliter.prefix_vertex->out_edges, edge, link.from.next)
    {
        target = edge->link.to.vertex;
        if (__glibc_unlikely(!target)) {
            continue;
        }

        // going straight from prefix -> suffix
        if (target == apply->read_graph.vertex_spliter.suffix_vertex) {
            if (bot_for_connect == NULL) {
                continue;
            }
            run_edge = assemble_graph_vertex_get_edge_between_vertex(top_for_connect, bot_for_connect);
            if (run_edge) {
                continue;
            }
            run_edge = assemble_graph_vertex_add_edge_to_vertex(apply->read_thread_graph, top_for_connect, bot_for_connect, edge->weight);
            run_edge->is_ref = edge->is_ref;
        }
        else {
            p_hc_shared_lib_pair_list_head vertex_st = apply->read_graph.vertex_spliter.middles;
            p_hc_shared_lib_pair_list_item vertex_item;

            CDL_FOREACH(vertex_st->nodes, vertex_item)
            {
                if (vertex_item->second_data != target) {
                    continue;
                }
                target = vertex_item->first_data;
                run_edge = assemble_graph_vertex_get_edge_between_vertex(top_for_connect, target);
                if (run_edge) {
                    continue;
                }
                run_edge = assemble_graph_vertex_add_edge_to_vertex(apply->read_thread_graph, top_for_connect, target, edge->weight);
                run_edge->is_ref = edge->is_ref;
                break;
            }
        }
    }
}

static void hc_assemble_vertex_sequence_spliter_add_edges_to_bottom_node(p_hc_apply apply, p_assemble_vertex bot_for_connect)
{
    p_assemble_edge edge, find_edge;
    p_assemble_vertex target;
    p_hc_shared_lib_pair_list_head vertex_st = apply->read_graph.vertex_spliter.middles;
    p_hc_shared_lib_pair_list_item vertex_item;

    CDL_FOREACH2(apply->read_graph.vertex_spliter.suffix_vertex->from_edges, edge, link.to.next)
    {
        target = edge->link.from.vertex;
        if (!target) {
            continue;
        }

        CDL_FOREACH(vertex_st->nodes, vertex_item)
        {
            if (vertex_item->second_data != target) {
                continue;
            }
            find_edge = assemble_graph_vertex_get_edge_between_vertex(vertex_item->first_data, bot_for_connect);
            if (!find_edge) {
                find_edge = assemble_graph_vertex_add_edge_to_vertex(apply->read_thread_graph, vertex_item->first_data, bot_for_connect,
                                                                     edge->weight);
                find_edge->is_ref = edge->is_ref;
            }
            break;
        }
    }
}
