/**
 * @file hc_assemble_seq_graph.c
 * @author Cao Ce (caoce@genomics.cn)
 * @brief A graph that contains base sequence at each node
 * @version 0.1
 * @date 2021-11-16
 *
 * @copyright Copyright (c) 2021 ROVACA SDK
 *
 */

#include "hc_assemble_seq_graph.h"

#include "debug.h"
#include "hc_assemble.h"
#include "hc_func.h"

HC_ASSEMBLE_SEQ_GRAPH_SIMLIFY_DECLARE(merge_diamonds)
// HC_ASSEMBLE_SEQ_GRAPH_SIMLIFY_DECLARE(merge_tails)
HC_ASSEMBLE_SEQ_GRAPH_SIMLIFY_DECLARE(split_common_suffices)
HC_ASSEMBLE_SEQ_GRAPH_SIMLIFY_DECLARE(merge_common_suffices)
HC_ASSEMBLE_SEQ_GRAPH_SIMLIFY(merge_diamonds, hc_assemble_vertex_sequence_spliter_merge_diamonds_try_to_transform)
// HC_ASSEMBLE_SEQ_GRAPH_SIMLIFY(merge_tails, hc_assemble_vertex_sequence_spliter_merge_tails_try_to_transform)
HC_ASSEMBLE_SEQ_GRAPH_SIMLIFY(split_common_suffices, hc_assemble_vertex_sequence_spliter_split_common_suffices_try_to_transform)
HC_ASSEMBLE_SEQ_GRAPH_SIMLIFY(merge_common_suffices, hc_assemble_vertex_sequence_spliter_split_merge_common_suffices)

/**
 * Convert this kmer graph to a simple sequence graph.
 *
 * Each kmer suffix shows up as a distinct SeqVertex, attached in the same structure as in the kmer
 * graph.  Nodes that are sources are mapped to SeqVertex nodes that contain all of their sequence
 *
 * @return a newly allocated SequenceGraph
 */
void hc_assemble_seq_graph_to_sequence_graph(p_assemble_graph graph)
{
    p_assemble_graph_hash_node vertex_node;
    p_assemble_vertex vertex;
    p_hc_shared_lib_list_head list_temp = hc_assemble_utils_check_list_work(graph->alt_path);
    p_hc_shared_lib_list_item list_el;

    // Run through the graph and clean up singular orphaned nodes
    CDL_FOREACH(graph->vertics_hash->node_list, vertex_node)
    {
        vertex = vertex_node->value;
        hc_shared_list_insert(list_temp, vertex);
    }

    CDL_FOREACH(list_temp->nodes, list_el)
    {
        vertex = list_el->data;

        assemble_graph_hash_remove_item(graph->vertics_hash, vertex->data, vertex->data_len, vertex);
        assemble_graph_hash_insert(graph->vertics_hash, &vertex, sizeof(void*), vertex);
        if (vertex->in_degree == 0) {
            continue;
        }
        vertex->data = vertex->data + vertex->data_len - 1;
        vertex->data_len = 1;
    }

    assemble_graph_to_seq_graph(graph);
    hc_shared_list_reset(list_temp);
}

/**
 * @brief
 *
 * @param list
 * @return p_hc_shared_lib_list_head
 */
static inline p_hc_shared_lib_list_head hc_assemble_seq_graph_check_list_work(p_hc_shared_lib_list_head list)
{
    if (list->nodes) {
        hc_shared_list_reset(list);
    }
    return list;
}

/**
 * Is source vertex potentially a start of a linear chain of vertices?
 *
 * We are a start of a zip chain if our out degree is 1 and either the
 * the vertex has no incoming connections or 2 or more (we must start a chain) or
 * we have exactly one incoming vertex and that one has out-degree > 1 (i.e., source's incoming
 * vertex couldn't be a start itself
 *
 * @param source a non-null vertex
 * @return true if source might start a linear chain
 */
static int hc_assemble_seq_graph_is_linear_chain_start(p_assemble_vertex vertex)
{
    p_assemble_edge edge;
    p_assemble_vertex search_vertex;

    if (vertex->out_degree != 1) {
        return 0;
    }
    if (vertex->in_degree != 1) {
        return 1;
    }

    CDL_FOREACH2(vertex->from_edges, edge, link.to.next)
    {
        search_vertex = edge->link.from.vertex;
        if (!search_vertex) {
            continue;
        }

        if (search_vertex->out_degree > 1) {
            return 1;
        }
    }
    return 0;
}

/**
 * @param v the vertex to test
 * @return  true if this vertex is a reference node (meaning that it appears on the reference path in the graph)
 */
static inline int hc_assemble_seq_graph_is_reference_node(p_assemble_vertex vertex)
{
    p_assemble_edge edge;

    CDL_FOREACH2(vertex->out_edges, edge, link.from.next)
    {
        if (edge->is_ref) {
            return 1;
        }
    }
    return 0;
}

/**
 * Get all of the vertices in a linear chain of vertices starting at zip_start
 *
 * Build a list of vertices (in order) starting from zip_start such that each sequential pair of vertices
 * in the chain A and B can be zipped together.
 *
 * @param zip_start a vertex that starts a linear chain
 * @return a list of vertices that comprise a linear chain starting with zip_start.  The resulting
 *         list will always contain at least zip_start as the first element.
 */
static p_hc_shared_lib_list_head hc_assemble_seq_graph_trace_linear_chain(p_assemble_vertex zip_start, p_hc_apply apply)
{
    p_hc_shared_lib_list_head linear_chain;
    p_assemble_vertex last, target;
    int last_is_ref, target_is_ref;

    linear_chain = hc_assemble_seq_graph_check_list_work(apply->read_graph.list_buffer_1);
    hc_shared_list_insert(linear_chain, zip_start);
    last_is_ref = hc_assemble_seq_graph_is_reference_node(zip_start);  // remember because this calculation is expensive
    last = zip_start;
    while (1) {
        if (last->out_degree != 1) {
            // cannot extend a chain from last if last has multiple outgoing branches
            break;
        }

        // there can only be one (outgoing edge of last) by contract
        target = last->out_edges->link.to.vertex;

        if (target->in_degree != 1 || last == target) {
            // cannot zip up a target that has multiple incoming nodes or that's a cycle to the last node
            break;
        }

        target_is_ref = hc_assemble_seq_graph_is_reference_node(target);
        if (last_is_ref != target_is_ref) {  // both our isRef states must be equal
            break;
        }

        hc_shared_list_insert(linear_chain, target);  // extend our chain by one

        // update our last state to be the current state, and continue
        last = target;
        last_is_ref = target_is_ref;
    }

    return linear_chain;
}

/**
 * @brief
 *
 * @param vertices
 * @param apply
 * @return p_assemble_vertex
 */
static p_assemble_vertex hc_assemble_seq_graph_merge_linear_chain_vertices(p_hc_shared_lib_list_head vertices, p_hc_apply apply)
{
    int32_t length = 0;
    p_hc_shared_lib_list_item list_item;
    p_assemble_vertex vertex;
    p_hc_assemble_graph_kmer kmer;

    // vertex长度
    CDL_FOREACH(vertices->nodes, list_item)
    {
        vertex = list_item->data;
        length += vertex->data_len;
    }
    length += 1;

    // kmer
    kmer = ring_mempool_auto_enlarge_malloc(apply->read_graph.kmer_buffer);
    memset(kmer, 0, sizeof(hc_assemble_graph_kmer));
    kmer->seq = hc_assemble_thread_mem_malloc(apply->read_graph.cache_pool, length);
    kmer->len = length - 1;

    // kmer->seq
    length = 0;
    CDL_FOREACH(vertices->nodes, list_item)
    {
        vertex = list_item->data;
        if (vertex->data_len != 1) {
            memcpy(kmer->seq + length, vertex->data, vertex->data_len);
        }
        else {
            kmer->seq[length] = *(uint8_t*)vertex->data;
        }
        length += vertex->data_len;
    }
    kmer->seq[length] = 0;

    // vertex
    vertex = assemble_graph_vertex_create(apply->read_thread_graph, kmer->seq, kmer->len, kmer);
    return vertex;
}

/**
 * Merge a linear chain of vertices into a single combined vertex, and update this graph to such that
 * the incoming edges into the first element of the linear_chain and the outgoing edges from linear_chain.getLast()
 * all point to this new combined vertex.
 *
 * @param linear_chain a non-empty chain of vertices that can be zipped up into a single vertex
 * @return true if we actually merged at least two vertices together
 */
static int hc_assemble_seq_graph_merge_linear_chain(p_hc_shared_lib_list_head linear_chain, p_hc_apply apply)
{
    p_hc_shared_lib_list_item first = linear_chain->nodes, last = linear_chain->nodes->prev, el;
    p_assemble_vertex added_vertex, vertex;
    p_assemble_edge edge, tmp1, tmp2;

    if (first == last) {
        return 0;  // only one element in the chain, cannot be extended
    }

    // create the combined vertex, and add it to the graph
    // TODO -- performance problem -- can be optimized if we want
    added_vertex = hc_assemble_seq_graph_merge_linear_chain_vertices(linear_chain, apply);
    added_vertex = assemble_graph_add_vertex(apply->read_thread_graph, added_vertex);

    // update the incoming and outgoing edges to point to the new vertex
    vertex = last->data;
    CDL_FOREACH_SAFE2(vertex->out_edges, edge, tmp1, tmp2, link.from.prev, link.from.next)
    {
        // edge to new vertex
        CDL_DELETE2(vertex->out_edges, edge, link.from.prev, link.from.next);
        edge->link.from.vertex = added_vertex;
        CDL_APPEND2(added_vertex->out_edges, edge, link.from.prev, link.from.next);
        added_vertex->out_degree++;
    }

    vertex = first->data;
    CDL_FOREACH_SAFE2(vertex->from_edges, edge, tmp1, tmp2, link.to.prev, link.to.next)
    {
        // edge to new vertex
        CDL_DELETE2(vertex->from_edges, edge, link.to.prev, link.to.next);
        edge->link.to.vertex = added_vertex;
        CDL_APPEND2(added_vertex->from_edges, edge, link.to.prev, link.to.next);
        added_vertex->in_degree++;
    }

    CDL_FOREACH(linear_chain->nodes, el) { assemble_graph_remove_vertex(apply->read_thread_graph, el->data); }
    return added_vertex ? 1 : 0;
}

/**
 * Zip up all of the simple linear chains present in this graph.
 * Merges together all pairs of vertices in the graph v1 -> v2 into a single vertex v' containing v1 + v2 sequence
 * Only works on vertices where v1's only outgoing edge is to v2 and v2's only incoming edge is from v1.
 * If such a pair of vertices is found, they are merged and the graph is update.  Otherwise nothing is changed.
 *
 * @return true if any such pair of vertices could be found, false otherwise
 */
int hc_assemble_seq_graph_zip_linear_chains(p_hc_apply apply)
{
    p_assemble_graph_hash_node vertex_node;
    p_assemble_vertex vertex;
    p_hc_shared_lib_list_head zipStarts, linear_chain;
    p_hc_shared_lib_list_item list_el;
    int mergedOne = 0;

    // create the list of start sites [doesn't modify graph yet]
    zipStarts = hc_assemble_seq_graph_check_list_work(apply->read_thread_graph->alt_path);
    CDL_FOREACH(apply->read_thread_graph->vertics_hash->node_list, vertex_node)
    {
        if (hc_assemble_seq_graph_is_linear_chain_start(vertex_node->value)) {
            hc_shared_list_insert(zipStarts, vertex_node->value);
        }
    }
    if (!zipStarts->nodes) {  // nothing to do, as nothing could start a chain
        return 0;
    }

    // At this point, zipStarts contains all of the vertices in this graph that might start some linear
    // chain of vertices.  We walk through each start, building up the linear chain of vertices and then
    // zipping them up with hc_assemble_seq_graph_merge_linear_chain, if possible
    // for ( final SeqVertex zip_start : zipStarts ) {
    CDL_FOREACH(zipStarts->nodes, list_el)
    {
        vertex = list_el->data;

        linear_chain = hc_assemble_seq_graph_trace_linear_chain(vertex, apply);
        if (!linear_chain->nodes) {
            continue;
        }
        // merge the linearized chain, recording if we actually did some useful work
        mergedOne |= hc_assemble_seq_graph_merge_linear_chain(linear_chain, apply);
        hc_shared_list_reset(linear_chain);
    }
    hc_shared_list_reset(zipStarts);

    return mergedOne;
}

/**
 * Remove all vertices in the graph that have in and out degree of 0
 */
void hc_assemble_seq_graph_remove_singleton_orphan_vertices(p_hc_apply apply)
{
    p_assemble_graph_hash_node vertex_node, tmp1, tmp2;
    p_assemble_vertex vertex;

    // Run through the graph and clean up singular orphaned nodes
    CDL_FOREACH_SAFE(apply->read_thread_graph->vertics_hash->node_list, vertex_node, tmp1, tmp2)
    {
        vertex = vertex_node->value;

        if (vertex->in_degree == 0 && vertex->out_degree == 0 && apply->read_thread_graph->vertices_sum != 1) {
            assemble_graph_remove_vertex(apply->read_thread_graph, vertex);
        }
    }
}

/**
 * @brief 搜索出一条路径后执行的操作
 *
 * @param last_vertex   上一个点
 * @param iter_st       遍历储存
 *
 * @return int 0
 */
static int hc_assemble_seq_graph_travel_one_path_done(p_assemble_vertex last_vertex, p_assemble_graph_iter iter_st)
{
    (void)last_vertex;
    p_hc_shared_lib_list_item path_head = iter_st->run_path;
    p_hc_shared_lib_list_item path_el;

    CDL_FOREACH(path_head, path_el)
    {
        p_assemble_vertex search_v = path_el->data;
        search_v->color = ASSEMBLE_GRAPH_VERTEX_COLOR_RED;
    }

    return 0;
}

/**
 * @brief 记录所有从ref vertex开始，dfs遍历的点
 *
 * @param ref_vertex    ref_vertex
 * @param apply         apply
 */
static void hc_assemble_seq_graph_travel_graph(p_assemble_vertex ref_vertex, p_hc_apply apply)
{
    p_assemble_graph_iter travel_config = &apply->read_thread_graph->graph_iter;

    travel_config->input_layer.from = ref_vertex;
    travel_config->input_layer.to = NULL;
    travel_config->input_layer.apply = apply;
    travel_config->cycle_detected = hc_assemble_utils_default_graph_func_cycle;
    travel_config->one_path_done = hc_assemble_seq_graph_travel_one_path_done;
    travel_config->all_path_done = hc_assemble_utils_default_graph_func_all_path_done;

    hc_assemble_base_graph_iter_dfs_non_recursion(travel_config);
}

/**
 * Remove all vertices on the graph that cannot be accessed by following any edge,
 * regardless of its direction, from the reference source vertex
 */
void hc_assemble_seq_graph_remove_vertices_not_connected_to_ref_regardless_of_edge_direction(p_hc_apply apply)
{
    p_assemble_graph_hash_node vertex_node, tmp1, tmp2;
    p_assemble_vertex vertex, ref_vertex;
    p_assemble_graph_hash_node node;

    node = hc_assemble_utils_get_reference_source_vertex_with_head(
        apply->read_thread_graph->vertics_hash->node_list, apply->read_thread_graph->vertics_hash->node_list, apply->read_thread_graph);
    if (__glibc_unlikely(node == NULL)) {
        return;
    }
    ref_vertex = node->value;
    CDL_FOREACH_SAFE(apply->read_thread_graph->vertics_hash->node_list, vertex_node, tmp1, tmp2)
    {
        vertex = vertex_node->value;
        vertex->color = ASSEMBLE_GRAPH_VERTEX_COLOR_NONE;
    }
    hc_assemble_seq_graph_travel_graph(ref_vertex, apply);

    CDL_FOREACH_SAFE(apply->read_thread_graph->vertics_hash->node_list, vertex_node, tmp1, tmp2)
    {
        vertex = vertex_node->value;
        if (vertex->color == ASSEMBLE_GRAPH_VERTEX_COLOR_NONE) {
            assemble_graph_remove_vertex(apply->read_thread_graph, vertex);
        }
    }
}

/**
 * Simplify this graph, merging vertices together and restructuring the graph in an
 * effort to minimize the number of overall vertices in the graph without changing
 * in any way the sequences implied by a complex enumeration of all paths through the graph.
 */
void hc_assemble_seq_graph_simplify_graph(p_hc_apply apply)
{
    int prev_graph_used = 0;
    int i = 0;

    // start off with one round of zipping of chains for performance reasons
    hc_assemble_seq_graph_zip_linear_chains(apply);

    for (; i < HC_ASSEMBLE_SEQ_GREAPH_MAX_REASONABLE_SIMPLIFICATION_CYCLES; i++) {
        // no simplification algorithm could run, so stop
        if (!hc_assemble_seq_graph_simplify_graph_once(apply)) {
            break;
        }

        // we get five cycles before we start looking for changes in the graph
        // by cloning ourselves and then checking for any changes
        if (i > HC_ASSEMBLE_SEQ_GREAPH_MAX_CHANGE_LOOPS) {
            // the previous graph and this graph have the same structure, so the simplification
            // algorithms are looping endless between states.  Just break and consider ourselves done
            if (prev_graph_used && hc_assemble_seq_graph_equal_graphs(apply->read_thread_graph)) {
                break;
            }

            hc_assemble_seq_graph_clone_graph(apply->read_thread_graph);
            prev_graph_used = 1;
        }
    }
}

/**
 * Run one full cycle of the graph simplification algorithms
 * @return true if any algorithms said they did some simplification
 */
static int hc_assemble_seq_graph_simplify_graph_once(p_hc_apply apply)
{
    // iterate until we haven't don't anything useful
    int didSomeWork = 0;

    hc_assemble_utils_check_hash_work(apply->read_graph.hash_buffer_1);
    apply->read_graph.hash_buffer_1->new_node_memcpy = hc_assemble_utils_hash_pt_new_node_memcpy;

    didSomeWork |= hc_assemble_seq_graph_merge_diamonds(apply);
    // didSomeWork |= hc_assemble_seq_graph_merge_tails(apply);
    didSomeWork |= hc_assemble_seq_graph_split_common_suffices(apply);
    didSomeWork |= hc_assemble_seq_graph_merge_common_suffices(apply);
    didSomeWork |= hc_assemble_seq_graph_zip_linear_chains(apply);

    apply->read_graph.hash_buffer_1->new_node_memcpy = hc_assemble_utils_hash_new_node_memcpy;
    assemble_graph_hash_reset(apply->read_graph.hash_buffer_1);

    return didSomeWork;
}

/**
 * @brief 复制seq graph
 *
 * @param from  from
 * @param apply apply
 */
static void hc_assemble_seq_graph_clone_graph(p_assemble_graph from)
{
    p_hc_shared_lib_list_head vertex_list, edge_list;
    p_assemble_graph_hash_node vertex_st;
    p_assemble_vertex vertex;
    p_assemble_edge from_edge, to_edge;
    uint32_t vertex_sum = 0, edge_sum = 0;

    vertex_list = hc_assemble_utils_check_list_work(from->simplify_graph.vertex);
    edge_list = hc_assemble_utils_check_list_work(from->simplify_graph.edge);
    ring_mempool_auto_enlarge_reset(from->simplify_graph.edge_pool);

    CDL_FOREACH(from->vertics_hash->node_list, vertex_st)
    {
        vertex = vertex_st->value;
        hc_shared_list_insert(vertex_list, vertex);
        vertex_sum += 1;
    }

    CDL_FOREACH2(from->edge_list, from_edge, link.all_st.next)
    {
        to_edge = ring_mempool_auto_enlarge_malloc(from->simplify_graph.edge_pool);
        memcpy(to_edge, from_edge, sizeof(assemble_edge));
        hc_shared_list_insert(edge_list, to_edge);
        edge_sum += 1;
    }
    from->simplify_graph.vertex_sum = vertex_sum;
    from->simplify_graph.edge_sum = edge_sum;
}

/**
 * @brief 比较seq vertex
 *
 * @param a         vertex a
 * @param b         vertex b
 *
 * @return int  1:failed, 0:succeed
 */
static inline int hc_assemble_seq_graph_compare_vertex(p_hc_shared_lib_list_item from_a, p_assemble_vertex b)
{
    p_assemble_vertex a = from_a->data;

    if (a->data_len != b->data_len) {
        return 1;
    }
    return memcmp(a->data, b->data, a->data_len);
}

/**
 * @brief 比较seq edge
 *
 * @param a         edge a
 * @param b         edge b
 *
 * @return int  1:failed, 0:succeed
 */
static inline int hc_assemble_seq_graph_compare_edge(p_hc_shared_lib_list_item from_a, p_assemble_edge b)
{
    p_assemble_edge a = from_a->data;
    p_assemble_vertex vertex_a, vertex_b;

    vertex_a = a->link.from.vertex;
    vertex_b = b->link.from.vertex;
    if (vertex_a->data_len != vertex_b->data_len || (memcmp(vertex_a->data, vertex_b->data, vertex_a->data_len) != 0)) {
        return 1;
    }

    vertex_a = a->link.to.vertex;
    vertex_b = b->link.to.vertex;
    if (vertex_a->data_len != vertex_b->data_len || (memcmp(vertex_a->data, vertex_b->data, vertex_a->data_len) != 0)) {
        return 1;
    }

    return 0;
}

/**
 * @brief graph相等
 *
 * @param graph     graph
 * @param apply     apply
 *
 * @return int      1:相等,0:不等
 */
static int hc_assemble_seq_graph_equal_graphs(p_assemble_graph graph)
{
    p_hc_shared_lib_list_head vertex_list, edge_list;
    p_hc_shared_lib_list_item list_el;
    p_assemble_graph_hash_node vertex_st;
    p_assemble_vertex vertex;
    p_assemble_edge edge;

    vertex_list = graph->simplify_graph.vertex;
    edge_list = graph->simplify_graph.edge;

    // vertex
    CDL_FOREACH(graph->vertics_hash->node_list, vertex_st)
    {
        vertex = vertex_st->value;

        list_el = NULL;
        CDL_SEARCH(vertex_list->nodes, list_el, vertex, hc_assemble_seq_graph_compare_vertex);
        if (!list_el) {
            return 0;
        }
    }

    // edge
    CDL_FOREACH2(graph->edge_list, edge, link.all_st.next)
    {
        list_el = NULL;
        CDL_SEARCH(edge_list->nodes, list_el, edge, hc_assemble_seq_graph_compare_edge);
        if (!list_el) {
            return 0;
        }
    }

    return 1;
}

/**
 * Make sure the reference sequence is properly represented in the provided graph
 *
 * @param graph the graph to check
 * @param refHaplotype the reference haplotype
 */
int hc_assemble_seq_graph_sanity_check_reference_graph(p_hc_apply apply)
{
    p_assemble_vertex ref_source, ref_sink;
    p_assemble_graph_hash_node ref_source_el, ref_sink_el;
    uint32_t ref_len;

    if (__glibc_unlikely(apply->read_thread_graph->vertices_sum == 1)) {
        ref_source = apply->read_thread_graph->vertics_hash->node_list->value;
        return !memcmp(ref_source->data, apply->read_thread_graph->origin_ref, ref_source->data_len);
    }

    ref_source_el = hc_assemble_utils_get_reference_source_vertex(
        apply->read_thread_graph->vertics_hash->node_list, apply->read_thread_graph->vertics_hash->node_list, apply->read_thread_graph);
    ref_sink_el = hc_assemble_utils_get_reference_sink_vertex(apply->read_thread_graph->vertics_hash->node_list,
                                                              apply->read_thread_graph->vertics_hash->node_list, apply->read_thread_graph);

    if (!ref_source_el || !ref_sink_el) {
        return 0;
    }
    ref_source = ref_source_el->value;
    ref_sink = ref_sink_el->value;
    ref_len = hc_assemble_seq_graph_get_reference_bytes(ref_source, ref_sink, apply);

    if (ref_len != apply->read_thread_graph->origin_ref_len ||
        memcmp(apply->read_thread_graph->seq, apply->read_thread_graph->origin_ref, ref_len) != 0) {
        return 0;
    }

    return 1;
}

/**
 * Walk along the reference path in the graph and pull out the corresponding bases
 *
 * @param fromVertex    starting vertex
 * @param toVertex      ending vertex
 *
 * @return              byte[] array holding the reference bases, this can be null if there are no nodes between the starting and ending
 * vertex (insertions for example)
 */
static uint32_t hc_assemble_seq_graph_get_reference_bytes(p_assemble_vertex fromVertex, p_assemble_vertex toVertex, p_hc_apply apply)
{
    uint8_t* bytes = apply->read_thread_graph->seq;
    uint32_t seq_len = 0;
    p_assemble_vertex vertex = fromVertex;

    apply->read_thread_graph->seq_len = 0;
    memcpy(bytes, vertex->data, vertex->data_len);
    bytes += vertex->data_len;
    seq_len += vertex->data_len;

    // advance along the reference path
    vertex = hc_assemble_seq_graph_get_next_reference_vertex(vertex);
    while (vertex != NULL && vertex != toVertex) {
        memcpy(bytes, vertex->data, vertex->data_len);
        bytes += vertex->data_len;
        seq_len += vertex->data_len;
        // advance along the reference path
        vertex = hc_assemble_seq_graph_get_next_reference_vertex(vertex);
    }
    if (vertex != NULL && vertex == toVertex) {
        memcpy(bytes, vertex->data, vertex->data_len);
        bytes += vertex->data_len;
        seq_len += vertex->data_len;
    }

    return seq_len;
}

/**
 * Traverse the graph and get the next reference vertex if it exists
 *
 * @param vertex the current vertex, can be null
 *
 * @return the next vertex (but not necessarily on the reference path if allowNonRefPaths is true) if it exists, otherwise null
 */
static p_assemble_vertex hc_assemble_seq_graph_get_next_reference_vertex(p_assemble_vertex vertex)
{
    p_assemble_edge edge;

    if (vertex == NULL || !vertex->out_edges) {
        return NULL;
    }

    CDL_FOREACH2(vertex->out_edges, edge, link.from.next)
    {
        if (!edge->is_ref || !edge->link.to.vertex) {
            continue;
        }
        return edge->link.to.vertex;
    }

    return NULL;
}
