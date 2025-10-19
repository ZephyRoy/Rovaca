/**
 * @file hc_assemble_base_graph.c
 * @author Cao Ce (caoce@genomics.cn)
 * @brief Graph基础计算函数
 * @version 0.1
 * @date 2021-11-03
 *
 * @copyright Copyright (c) 2021 ROVACA SDK
 *
 */

#include "hc_assemble_base_graph.h"

#include "debug.h"
#include "hc_apply.h"
#include "hc_assemble.h"
#include "hc_func.h"
#include "hc_marco.h"
#include "list.h"

/**
 * @brief 搜索出一条路径后执行的操作
 *
 * @param last_vertex   上一个点
 * @param iter_st       遍历储存
 *
 * @return int 0
 */
static int hc_assemble_base_graph_iter_vertex_one_path_done(p_assemble_vertex last_vertex, p_assemble_graph_iter iter_st)
{
    (void)last_vertex;
    p_hc_shared_lib_list_item path_head = iter_st->run_path;
    p_hc_shared_lib_list_item path_el;

    CDL_FOREACH(path_head, path_el)
    {
        p_assemble_vertex vertex = path_el->data;
        if (!vertex->color) {
            vertex->color = ASSEMBLE_GRAPH_VERTEX_COLOR_RED;
        }
    }

    return 0;
}

static int hc_assemble_base_graph_iter_vertex(p_hc_apply apply)
{
    p_assemble_graph_hash_node vertex_node;
    p_assemble_graph_iter travel_config = &apply->read_thread_graph->graph_iter;
    p_assemble_vertex vertex;
    int ret = 0, done = 0;

    CDL_FOREACH(apply->read_thread_graph->vertics_hash->node_list, vertex_node)
    {
        vertex = vertex_node->value;
        vertex->color = ASSEMBLE_GRAPH_VERTEX_COLOR_NONE;
    }

    CDL_FOREACH(apply->read_thread_graph->vertics_hash->node_list, vertex_node)
    {
        vertex = vertex_node->value;
        if (!vertex->out_degree || vertex->color) {
            continue;
        }

        travel_config->input_layer.from = vertex;
        travel_config->input_layer.to = NULL;
        travel_config->input_layer.apply = apply;
        travel_config->cycle_detected = hc_assemble_utils_default_graph_func_cycle;
        travel_config->one_path_done = hc_assemble_base_graph_iter_vertex_one_path_done;
        travel_config->all_path_done = NULL;

        ret += hc_assemble_base_graph_iter_dfs_non_recursion(travel_config);
        if (__glibc_unlikely(ret > 0)) {
            break;
        }
        done += 1;
    }
    if (__glibc_unlikely(!done && apply->read_thread_graph->vertices_sum > 1)) {
        ret = 1;
    }

    return ret;
}

/**
 * @brief 搜索Graph内Cycle
 *
 * @param apply apply
 *
 * @return int 0:succeed, 1:cycle
 */
int hc_assemble_base_graph_cycle_finder(p_hc_apply apply)
{
    return hc_assemble_base_graph_iter_vertex(apply);
#if 0
    int ret;
    p_assemble_graph_iter travel_config = &apply->read_thread_graph->graph_iter;

    travel_config->input_layer.from         = apply->read_graph.ref_source;
    travel_config->input_layer.to           = NULL;
    travel_config->input_layer.apply        = apply;
    travel_config->cycle_detected           = hc_assemble_utils_default_graph_func_cycle;
    travel_config->one_path_done            = NULL;
    travel_config->all_path_done            = NULL;

    ret = hc_assemble_base_graph_iter_dfs_non_recursion(travel_config);

    return ret;
#endif
}

/**
 * Finds the path upwards in the graph from this vertex to the first diverging node, including that (lowest common ancestor) vertex.
 * Note that nodes are excluded if their pruning weight is less than the pruning factor.
 *
 * @param vertex         the original vertex
 * @param pruneFactor    the prune factor to use in ignoring chain pieces if edge multiplicity is < pruneFactor
 * @param giveUpAtBranch stop trying to find a path if a vertex with multiple incoming or outgoing edge is found
 * @return the path if it can be determined or null if this vertex either doesn't merge onto another path or
 * has an ancestor with multiple incoming edges before hitting the reference path
 */
static p_hc_shared_lib_list_head hc_assemble_base_graph_find_path_upwards_to_lowest_common_ancestor(p_assemble_vertex vertex,
                                                                                                    p_hc_apply apply)
{
    p_assemble_graph_hash_table visited_nodes = hc_assemble_utils_check_hash_work(apply->read_graph.hash_buffer_1);

    p_hc_shared_lib_list_head path = hc_assemble_utils_check_list_work(apply->read_thread_graph->alt_path);
    p_assemble_vertex run_vetex = vertex;
    p_assemble_graph_hash_node hash_node;

    visited_nodes->new_node_memcpy = hc_assemble_utils_hash_pt_new_node_memcpy;

    while (run_vetex && run_vetex->in_degree && !(run_vetex->in_degree != 1 || run_vetex->out_degree >= 2)) {
        p_assemble_edge edge = run_vetex->from_edges;
        uint32_t key = run_vetex->key;

        // if it has too low a weight, don't use it (or previous vertexes) for the path
        if (edge->weight < HC_ASSEMBLE_GRAPH_CHAIN_PRUNE_FACTOR) {
            hc_shared_list_reset(path);
        }
        else {
            hc_shared_list_insert_prev(path, run_vetex);
        }
        run_vetex = edge->link.from.vertex;

        // Check that we aren't stuck in a loop
        hash_node = assemble_graph_hash_search_pt_with_key(visited_nodes, run_vetex, key);
        if (__glibc_unlikely(hash_node && hash_node->value == run_vetex)) {
            hc_shared_list_reset(path);
            assemble_graph_hash_reset(visited_nodes);
            return NULL;
        }
        else {
            assemble_graph_hash_insert_with_key(visited_nodes, &run_vetex, sizeof(void*), run_vetex, key);
        }
    }
    assemble_graph_hash_reset(visited_nodes);
    visited_nodes->new_node_memcpy = hc_assemble_utils_hash_new_node_memcpy;

    if (run_vetex) {
        hc_shared_list_insert_prev(path, run_vetex);
    }
    if (run_vetex && run_vetex->out_degree > 1) {
        return path;
    }
    else {
        hc_shared_list_reset(path);
        return NULL;
    }
}

static p_assemble_edge hc_assemble_base_graph_get_heaviest_incoming_edge(p_assemble_vertex vertex)
{
    p_assemble_edge edge, res = NULL;

    CDL_FOREACH2(vertex->from_edges, edge, link.to.next)
    {
        if (edge->weight == 1) {
            continue;
        }
        if (!res || edge->weight > res->weight) {
            res = edge;
        }
    }
    return res;
}

static p_assemble_edge hc_assemble_base_graph_get_heaviest_outgoing_edge(p_assemble_vertex vertex)
{
    p_assemble_edge edge, res = NULL;

    CDL_FOREACH2(vertex->out_edges, edge, link.from.next)
    {
        if (edge->weight == 1) {
            continue;
        }
        if (!res || edge->weight > res->weight) {
            res = edge;
        }
    }
    return res;
}

/**
 * Traverse the graph and get the next reference vertex if it exists
 * @param vertex the current vertex, can be null
 * @param allowNonRefPaths if true, allow sub-paths that are non-reference if there is only a single outgoing edge
 * @param blacklistedEdge optional edge to ignore in the traversal down; useful to exclude the non-reference dangling paths
 * @return the next vertex (but not necessarily on the reference path if allowNonRefPaths is true) if it exists, otherwise null
 */
p_assemble_vertex hc_assemble_base_graph_get_next_reference_vertex(p_assemble_vertex vertex, int allowNonRefPaths,
                                                                   p_assemble_edge blacklistedEdge)
{
    p_assemble_edge outgoingEdges, edge_el, edge_ret = NULL;
    p_assemble_vertex ret;
    uint32_t edge_sum = 0;

    if (vertex == NULL) {
        return NULL;
    }

    outgoingEdges = vertex->out_edges;
    if (!outgoingEdges) {
        return NULL;
    }

    CDL_FOREACH2(outgoingEdges, edge_el, link.from.next)
    {
        if (edge_el->is_ref) {
            ret = edge_el->link.to.vertex;
            return ret;
        }
    }

    if (!allowNonRefPaths) {
        return NULL;
    }

    // if we got here, then we aren't on a reference path
    CDL_FOREACH2(outgoingEdges, edge_el, link.from.next)
    {
        if (edge_el == blacklistedEdge) {
            continue;
        }
        if (!edge_ret) {
            edge_ret = edge_el;
        }
        edge_sum += 1;
        if (edge_sum >= 2) {
            break;
        }
    }
    if (edge_sum == 1) {
        return edge_ret->link.to.vertex;
    }
    else {
        return NULL;
    }
}

/**
 * Traverse the graph and get the previous reference vertex if it exists
 * @param vertex the current vertex, can be null
 * @return  the previous reference vertex if it exists or null otherwise.
 */
p_assemble_vertex hc_assemble_base_graph_get_prev_reference_vertex(p_assemble_vertex vertex)
{
    p_assemble_edge edge_el;

    if (vertex == NULL) {
        return NULL;
    }

    CDL_FOREACH2(vertex->from_edges, edge_el, link.to.next)
    {
        if (edge_el->is_ref) {
            return edge_el->link.from.vertex;
        }
    }
    return NULL;
}

/**
 * Finds the path in the graph from this vertex to the reference sink, including this vertex
 *
 * @param start           the reference vertex to start from
 * @param direction       describes which direction to move in the graph (i.e. down to the reference sink or up to the source)
 * @param blacklistedEdge edge to ignore in the traversal down; useful to exclude the non-reference dangling paths
 * @return the path (non-null, non-empty)
 */
static int hc_assemble_base_graph_get_reference_path(p_assemble_vertex start, int direction, p_assemble_edge blacklistedEdge,
                                                     p_hc_apply apply)
{
    p_hc_shared_lib_list_head path = hc_assemble_utils_check_list_work(apply->read_thread_graph->ref_path);
    p_assemble_vertex vertex = start;

    while (vertex != NULL) {
        if (__glibc_unlikely(hc_shared_list_search(path, vertex))) {
            hc_shared_list_reset(path);
            return 0;
        }
        hc_shared_list_insert(path, vertex);
        vertex = (direction == ASSEMBLE_GRAPH_DIRECTION_TO ? hc_assemble_base_graph_get_next_reference_vertex(vertex, 1, blacklistedEdge)
                                                           : hc_assemble_base_graph_get_prev_reference_vertex(vertex));
    }

    return path->nodes ? 1 : 0;
}

/**
 * The base sequence for the given path.
 *
 * @param path         the list of vertexes that make up the path
 * @param expandSource if true and if we encounter a source node, then expand (and reverse) the character sequence for that node
 * @return non-null sequence of bases corresponding to the given path
 */
void hc_assemble_base_graph_get_bases_for_path(p_hc_shared_lib_list_head ref, p_hc_shared_lib_list_head alt, int expandSource,
                                               p_hc_apply apply)
{
    p_hc_shared_lib_list_item path_item;
    p_assemble_vertex vertex_el;
    uint32_t i;

    apply->read_thread_graph->seq_len = 0;
    apply->read_thread_graph->ref_len = 0;

    // ref
    CDL_FOREACH(ref->nodes, path_item)
    {
        vertex_el = path_item->data;
        if (__glibc_unlikely(vertex_el->in_degree == 0 && expandSource)) {
            uint8_t* data = vertex_el->data;
            for (i = 0; i < vertex_el->data_len; i++) {
                apply->read_thread_graph->ref[apply->read_thread_graph->ref_len] = data[vertex_el->data_len - 1 - i];
                apply->read_thread_graph->ref_len += 1;
            }
        }
        else {
            apply->read_thread_graph->ref[apply->read_thread_graph->ref_len] = hc_assemble_graph_get_suffix(vertex_el);
            apply->read_thread_graph->ref_len += 1;
        }
    }
    apply->read_thread_graph->ref[apply->read_thread_graph->ref_len] = 0;

    // alt
    CDL_FOREACH(alt->nodes, path_item)
    {
        vertex_el = path_item->data;
        if (__glibc_unlikely(vertex_el->in_degree == 0 && expandSource)) {
            uint8_t* data = vertex_el->data;
            for (i = 0; i < vertex_el->data_len; i++) {
                apply->read_thread_graph->seq[apply->read_thread_graph->seq_len] = data[vertex_el->data_len - 1 - i];
                apply->read_thread_graph->seq_len += 1;
            }
        }
        else {
            apply->read_thread_graph->seq[apply->read_thread_graph->seq_len] = hc_assemble_graph_get_suffix(vertex_el);
            apply->read_thread_graph->seq_len += 1;
        }
    }
    apply->read_thread_graph->seq[apply->read_thread_graph->seq_len] = 0;
}

/**
 * Generates the CIGAR string from the Smith-Waterman alignment of the dangling path (where the
 * provided vertex is the sink) and the reference path.
 *
 * @param aligner
 * @param vertex      the sink of the dangling chain
 * @param pruneFactor the prune factor to use in ignoring chain pieces if edge multiplicity is < pruneFactor
 * @param recoverAll  recover even branches with forks
 * @return a SmithWaterman object which can be null if no proper alignment could be generated
 */
static int hc_assemble_base_graph_generate_cigar_against_downwards_reference_path(p_assemble_vertex vertex, p_hc_apply apply)
{
    p_hc_shared_lib_list_item path_item;
    int alt_size = 0;
    p_hc_shared_lib_list_head path;

    // find the lowest common ancestor path between this vertex and the diverging master path if available
    path = hc_assemble_base_graph_find_path_upwards_to_lowest_common_ancestor(vertex, apply);
    if (path) {
        CDL_COUNT(path->nodes, path_item, alt_size);
    }
    else {
        return 0;
    }
    if (alt_size < HC_ASSEMBLE_BASE_GRAPH_MIN_DANGLING_BRANCH_LENGTH + 1 ||
        hc_assemble_base_graph_is_ref_source(path->nodes->data, apply->read_thread_graph)) {  // add 1 to include the LCA
        apply->gatk_sw.output.cigar_len = 0;
        return 0;
    }

    // now get the reference path from the LCA
    if (!hc_assemble_base_graph_get_reference_path(path->nodes->data, ASSEMBLE_GRAPH_DIRECTION_TO,
                                                   hc_assemble_base_graph_get_heaviest_incoming_edge(path->nodes->next->data), apply)) {
        hc_shared_list_reset(path);
        return 0;
    }

    hc_assemble_base_graph_get_bases_for_path(apply->read_thread_graph->ref_path, path, 0, apply);

    // run Smith-Waterman to determine the best alignment (and remove trailing deletions since they aren't interesting)
    hc_assemble_base_graph_run_sw(apply);
    hc_assemble_base_graph_remove_trailing_deletions(apply);
    return 1;
}

/**
 * @param v the vertex to test
 * @return  true if this vertex is a reference node (meaning that it appears on the reference path in the graph)
 */
static inline int hc_assemble_base_graph_is_reference_node(p_assemble_vertex vertex)
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
 * Finds the path downwards in the graph from this vertex to the reference sequence, including the highest common descendant vertex.
 * However note that the path is reversed so that this vertex ends up at the end of the path.
 * Also note that nodes are excluded if their pruning weight is less than the pruning factor.
 *
 * @param vertex         the original vertex
 * @param pruneFactor    the prune factor to use in ignoring chain pieces
 * @param giveUpAtBranch stop trying to find a path if a vertex with multiple incoming or outgoing edge is found
 * @return the path if it can be determined or null if this vertex either doesn't merge onto the reference path or
 * has a descendant with multiple outgoing edges before hitting the reference path
 */
static int hc_assemble_base_graph_find_path_downwards_to_highest_common_descendant_of_reference(p_assemble_vertex vertex, p_hc_apply apply)
{
    p_assemble_graph_hash_table visited_nodes = hc_assemble_utils_check_hash_work(apply->read_graph.hash_buffer_1);
    p_hc_shared_lib_list_head path = apply->read_thread_graph->alt_path;
    p_assemble_vertex run_vetex = vertex;
    p_assemble_graph_hash_node hash_node;

    visited_nodes->new_node_memcpy = hc_assemble_utils_hash_pt_new_node_memcpy;

    while (run_vetex && run_vetex->out_degree && !(hc_assemble_base_graph_is_reference_node(run_vetex) || run_vetex->out_degree != 1)) {
        uint32_t key = run_vetex->key;
        p_assemble_edge edge = run_vetex->out_edges;

        // if it has too low a weight, don't use it (or previous vertexes) for the path
        if (edge->weight < HC_ASSEMBLE_GRAPH_CHAIN_PRUNE_FACTOR) {
            hc_shared_list_reset(path);
        }
        else {
            hc_shared_list_insert_prev(path, run_vetex);
        }

        run_vetex = edge->link.to.vertex;

        // Check that we aren't stuck in a loop
        hash_node = assemble_graph_hash_search_pt_with_key(visited_nodes, run_vetex, key);
        if (__glibc_unlikely(hash_node && hash_node->value == run_vetex)) {
            hc_shared_list_reset(path);
            assemble_graph_hash_reset(visited_nodes);
            return 0;
        }
        else {
            assemble_graph_hash_insert_with_key(visited_nodes, &run_vetex, sizeof(void*), run_vetex, key);
        }
    }
    assemble_graph_hash_reset(visited_nodes);
    visited_nodes->new_node_memcpy = hc_assemble_utils_hash_new_node_memcpy;

    if (run_vetex) {
        hc_shared_list_insert_prev(path, run_vetex);
    }
    if (run_vetex && hc_assemble_base_graph_is_reference_node(run_vetex)) {
        return 1;
    }
    else {
        hc_shared_list_reset(path);
        return 0;
    }
}

/**
 * Finds the path downwards in the graph from this vertex to the reference sequence, including the highest common descendant vertex.
 * However note that the path is reversed so that this vertex ends up at the end of the path.
 * Also note that nodes are excluded if their pruning weight is less than the pruning factor.
 *
 * @param vertex         the original vertex
 * @param pruneFactor    the prune factor to use in ignoring chain pieces
 * @param giveUpAtBranch stop trying to find a path if a vertex with multiple incoming or outgoing edge is found
 * @return the path if it can be determined or null if this vertex either doesn't merge onto the reference path or
 * has a descendant with multiple outgoing edges before hitting the reference path
 */
static int hc_assemble_base_graph_generate_cigar_against_upwards_reference_path(p_assemble_vertex vertex, p_hc_apply apply)
{
    p_hc_shared_lib_list_item path_item;
    uint32_t alt_size = 0;

    hc_shared_list_reset(apply->read_thread_graph->ref_path);
    hc_shared_list_reset(apply->read_thread_graph->alt_path);

    // find the lowest common ancestor path between this vertex and the diverging master path if available
    if (hc_assemble_base_graph_find_path_downwards_to_highest_common_descendant_of_reference(vertex, apply)) {
        CDL_COUNT(apply->read_thread_graph->alt_path->nodes, path_item, alt_size);
    }
    if (apply->read_thread_graph->alt_path->nodes == NULL ||
        hc_assemble_base_graph_is_ref_sink(apply->read_thread_graph->alt_path->nodes->data, apply->read_thread_graph) ||
        alt_size < HC_ASSEMBLE_BASE_GRAPH_MIN_DANGLING_BRANCH_LENGTH + 1) {  // add 1 to include the LCA
        apply->gatk_sw.output.cigar_len = 0;
        return 0;
    }

    // now get the reference path from the LCA
    if (!hc_assemble_base_graph_get_reference_path(apply->read_thread_graph->alt_path->nodes->data, ASSEMBLE_GRAPH_DIRECTION_FROM, NULL,
                                                   apply)) {
        return 0;
    }

    // create the Smith-Waterman strings to use
    hc_assemble_base_graph_get_bases_for_path(apply->read_thread_graph->ref_path, apply->read_thread_graph->alt_path, 1, apply);

    // run Smith-Waterman to determine the best alignment (and remove trailing deletions since they aren't interesting)
    hc_assemble_base_graph_run_sw(apply);
    hc_assemble_base_graph_remove_trailing_deletions(apply);
    apply->read_thread_graph->alt_size = alt_size;
    alt_size = 0;
    CDL_COUNT(apply->read_thread_graph->ref_path->nodes, path_item, alt_size);
    apply->read_thread_graph->ref_size = alt_size;
    return 1;
}

static void hc_assemble_base_graph_run_sw(p_hc_apply apply)
{
    p_gatk_sw_storage sw = &apply->gatk_sw;
    gatk_sw_parameters param = GATK_SW_STANDARD_NGS;

    hc_assemble_gatk_sw_init(sw);

    sw->input.ref = (void*)apply->read_thread_graph->ref;
    sw->input.ref_len = apply->read_thread_graph->ref_len;
    sw->input.alt = (void*)apply->read_thread_graph->seq;
    sw->input.alt_len = apply->read_thread_graph->seq_len;
    sw->input.overhang_strategy = GATK_SW_OVERHANG_STRATEGY_LEADING_INDEL;
    sw->paramates = param;

    hc_assemble_gatk_sw_align(sw);
}

/**
 * Removing a trailing deletion from the incoming cigar if present
 *
 * @param c the cigar we want to update
 * @return a non-null Cigar
 */
static void hc_assemble_base_graph_remove_trailing_deletions(p_hc_apply apply)
{
    uint32_t cigar = apply->gatk_sw.output.cigar[apply->gatk_sw.output.cigar_len - 1];

    if (__glibc_likely(bam_cigar_op(cigar) != WORKER_SIGAR_STATUS_DELETION)) {
        return;
    }
    apply->gatk_sw.output.cigar_len -= 1;
}

/**
 * Attempt to attach vertex with out-degree == 0 to the graph
 *
 * @param vertex                  the vertex to recover
 * @param pruneFactor             the prune factor to use in ignoring chain pieces if edge multiplicity is < pruneFactor
 * @param minDanglingBranchLength the minimum length of a dangling branch for us to try to merge it
 * @param aligner
 * @return 1 if we successfully recovered the vertex and 0 otherwise
 */
static int hc_assemble_base_graph_recover_dangling_tail(p_assemble_vertex vertex, p_hc_apply apply)
{
    // generate the CIGAR string from Smith-Waterman between the dangling tail and reference paths
    if (!hc_assemble_base_graph_generate_cigar_against_downwards_reference_path(vertex, apply)) {
        return 0;
    }

    // if the CIGAR is too complex (or couldn't be computed) then we do not allow the merge into the reference path
    if (!apply->gatk_sw.output.cigar_len || !hc_assemble_base_graph_cigar_is_okay_to_merge(apply, 0, 1)) {
        return 0;
    }

    // merge
    return hc_assemble_base_graph_merge_dangling_tail(apply);
}

/**
 * Attempt to attach vertex with in-degree == 0, or a vertex on its path, to the graph
 *
 * @param vertex                  the vertex to recover
 * @param pruneFactor             the prune factor to use in ignoring chain pieces if edge multiplicity is < pruneFactor
 * @param minDanglingBranchLength the minimum length of a dangling branch for us to try to merge it
 * @param recoverAll              recover even branches with forks
 * @param aligner
 * @return 1 if we successfully recovered a vertex and 0 otherwise
 */
static int hc_assemble_base_graph_recover_dangling_head(p_assemble_vertex vertex, p_hc_apply apply)
{
    // generate the CIGAR string from Smith-Waterman between the dangling head and reference paths
    if (!hc_assemble_base_graph_generate_cigar_against_upwards_reference_path(vertex, apply)) {
        return 0;
    }

    // if the CIGAR is too complex (or couldn't be computed) then we do not allow the merge into the reference path
    if (!apply->gatk_sw.output.cigar_len || !hc_assemble_base_graph_cigar_is_okay_to_merge(apply, 1, 0)) {
        return 0;
    }

    // merge
    return hc_assemble_base_graph_merge_dangling_head_legacy(apply);
}

/**
 * calculates the longest suffix match between a sequence and a smaller kmer
 *
 * @param seq      the (reference) sequence
 * @param kmer     the smaller kmer sequence
 * @param seqStart the index (inclusive) on seq to start looking backwards from
 * @return the longest matching suffix
 */
static int hc_assemble_base_graph_longest_suffix_match(uint8_t* seq, uint8_t* kmer, int seqStart, int kmer_len)
{
    int len = 1;

    for (; len <= kmer_len; len++) {
        int seqI = seqStart - len + 1;
        int kmerI = kmer_len - len;
        if (seqI < 0 || seq[seqI] != kmer[kmerI]) {
            return len - 1;
        }
    }
    return kmer_len;
}

/**
 * Actually merge the dangling tail if possible
 *
 * @param danglingTailMergeResult the result from generating a Cigar for the dangling tail against the reference
 * @return 1 if merge was successful, 0 otherwise
 */
static int hc_assemble_base_graph_merge_dangling_tail(p_hc_apply apply)
{
    p_assemble_vertex alt, ref;
    p_assemble_edge edge;
    uint32_t* elements = apply->gatk_sw.output.cigar;
    uint32_t lastElement = elements[apply->gatk_sw.output.cigar_len - 1];
    uint32_t lastRefIndex = hc_assemble_utils_get_ref_length(apply->gatk_sw.output.cigar, apply->gatk_sw.output.cigar_len) - 1;
    int matchingSuffix =
        assemble_min(hc_assemble_base_graph_longest_suffix_match((void*)apply->gatk_sw.input.ref, (void*)apply->gatk_sw.input.alt,
                                                                 lastRefIndex, apply->gatk_sw.input.alt_len),
                     (int)bam_cigar_oplen(lastElement));

    if (matchingSuffix == 0) {
        return 0;
    }
    int altIndexToMerge =
        assemble_max(hc_assemble_utils_get_cigar_read_length(elements, apply->gatk_sw.output.cigar_len) - matchingSuffix - 1, 0);

    // there is an important edge condition that we need to handle here: Smith-Waterman correctly calculates that there is a
    // deletion, that deletion is left-aligned such that the LCA node is part of that deletion, and the rest of the dangling
    // tail is a perfect match to the suffix of the reference path.  In this case we need to push the reference index to merge
    // down one position so that we don't incorrectly cut a base off of the deletion.
    int firstElementIsDeletion = bam_cigar_op(elements[0]) == WORKER_SIGAR_STATUS_DELETION;
    int mustHandleLeadingDeletionCase = firstElementIsDeletion && (bam_cigar_oplen(elements[0]) + matchingSuffix == lastRefIndex + 1);
    int refIndexToMerge = lastRefIndex - matchingSuffix + 1 + (mustHandleLeadingDeletionCase ? 1 : 0);

    // another edge condition occurs here: if Smith-Waterman places the whole tail into an insertion then it will try to
    // merge back to the LCA, which results in a cycle in the graph.  So we do not want to merge in such a case.
    if (refIndexToMerge == 0) {
        return 0;
    }

    // it's safe to merge now
    HC_ASSEMBLE_GRAPH_GET_PATH_INDEX_ITEM(apply->read_thread_graph->alt_path->nodes, (uint32_t)altIndexToMerge, alt);
    HC_ASSEMBLE_GRAPH_GET_PATH_INDEX_ITEM(apply->read_thread_graph->ref_path->nodes, (uint32_t)refIndexToMerge, ref);
    if (!ref || !alt) {
        return 0;
    }
    edge = assemble_graph_vertex_get_edge_between_vertex(alt, ref);
    if (__glibc_likely(!edge)) {
        edge = assemble_graph_vertex_add_edge_to_vertex(apply->read_thread_graph, alt, ref, 1);
        edge->is_ref = 0;
    }

    return 1;
}

/**
 * NOTE: this method is only used for dangling heads and not tails.
 *
 * Determine the maximum number of mismatches permitted on the branch.
 * Unless it's preset (e.g. by unit tests) it should be the length of the branch divided by the kmer size.
 *
 * @param lengthOfDanglingBranch the length of the branch itself
 * @return positive integer
 */
static inline int hc_assemble_base_graph_get_max_mismatches_legacy(int lengthOfDanglingBranch, int kmerSize)
{
    return assemble_max(1, (lengthOfDanglingBranch / kmerSize));
}

/**
 * Finds the index of the best extent of the prefix match between the provided paths, for dangling head merging.
 * Assumes that path1.length >= maxIndex and path2.length >= maxIndex.
 *
 * @param path1    the first path
 * @param path2    the second path
 * @param maxIndex the maximum index to traverse (not inclusive)
 * @return the index of the ideal prefix match or -1 if it cannot find one, must be less than maxIndex
 */
static int hc_assemble_base_graph_best_prefix_match_legacy(p_hc_apply apply, int maxIndex)
{
    int maxMismatches = hc_assemble_base_graph_get_max_mismatches_legacy(maxIndex, apply->read_graph.kmer);
    int mismatches = 0;
    int index = 0;
    int lastGoodIndex = -1;
    while (index < maxIndex) {
        if (apply->gatk_sw.input.ref[index] != apply->gatk_sw.input.alt[index]) {
            if (++mismatches > maxMismatches) {
                return -1;
            }
            lastGoodIndex = index;
        }
        index++;
    }
    // if we got here then we hit the max index
    return lastGoodIndex;
}

static int hc_assemble_base_graph_get_offset_for_ref_end_to_dangling_end(uint32_t* cigar, uint32_t cigar_len)
{
    register uint32_t i = 0;
    int ret = 0;
    uint32_t cigar_op, cigar_op_len;

    for (; i < cigar_len; i++) {
        cigar_op = bam_cigar_op(cigar[i]);
        cigar_op_len = bam_cigar_oplen(cigar[i]);

        if (consumes_ref_bases(cigar_op)) {
            ret += cigar_op_len;
        }
        if (consumes_read_bases(cigar_op)) {
            ret -= cigar_op_len;
        }
    }
    return ret;
}

static int hc_assemble_base_graph_extend_dangling_path_against_reference(p_hc_apply apply, int numNodesToExtend, uint32_t* elements)
{
    (void)elements;
    p_assemble_vertex danglingSource, refSourceSequence, newV, prevV;
    p_hc_assemble_graph_kmer kmer;
    p_assemble_edge sourceEdge, remove_edge;
    char* seq;
    int offsetForRefEndToDanglingEnd =
        hc_assemble_base_graph_get_offset_for_ref_end_to_dangling_end(apply->gatk_sw.output.cigar, apply->gatk_sw.output.cigar_len);
    uint32_t indexOfRefNodeToUse = apply->read_thread_graph->alt_size - 1 + offsetForRefEndToDanglingEnd + numNodesToExtend;
    int i;

    if (indexOfRefNodeToUse >= apply->read_thread_graph->ref_size) {
        return 0;
    }

    danglingSource = apply->read_thread_graph->alt_path->nodes->prev->data;
    hc_shared_list_delete(apply->read_thread_graph->alt_path, danglingSource);
    apply->read_thread_graph->alt_size -= 1;
    refSourceSequence = hc_assemble_utils_get_list_index_item(apply->read_thread_graph->ref_path, indexOfRefNodeToUse);

    seq = hc_assemble_thread_mem_malloc(apply->read_graph.cache_pool, numNodesToExtend + apply->gatk_sw.input.alt_len + 1);
    sprintf(seq, "%.*s%.*s", numNodesToExtend, (char*)refSourceSequence->data, danglingSource->data_len, (char*)danglingSource->data);

    // clean up the source and edge
    sourceEdge = hc_assemble_base_graph_get_heaviest_outgoing_edge(danglingSource);
    prevV = sourceEdge->link.to.vertex;
    remove_edge = assemble_graph_vertex_get_edge_between_vertex(danglingSource, prevV);
    if (likely(remove_edge)) {
        assemble_graph_remove_edge(apply->read_thread_graph, remove_edge);
    }

    // extend the path
    for (i = numNodesToExtend; i > 0; i--) {
        // kmer
        kmer = ring_mempool_auto_enlarge_malloc(apply->read_graph.kmer_buffer);
        memset(kmer, 0, sizeof(hc_assemble_graph_kmer));
        kmer->len = apply->read_graph.kmer;
        kmer->seq = (uint8_t*)(seq + i);

        // newV
        newV = assemble_graph_vertex_create(apply->read_thread_graph, kmer->seq, kmer->len, kmer);
        assemble_graph_add_vertex(apply->read_thread_graph, newV);

        // edge
        assemble_graph_vertex_add_edge_to_vertex(apply->read_thread_graph, newV, prevV, sourceEdge->weight);
        hc_shared_list_insert(apply->read_thread_graph->alt_path, newV);
        apply->read_thread_graph->alt_size += 1;
        prevV = newV;
    }

    return 1;
}

/**
 * Actually merge the dangling head if possible, this is the old codepath that does not handle indels
 *
 * @param danglingHeadMergeResult   the result from generating a Cigar for the dangling head against the reference
 * @return 1 if merge was successful, 0 otherwise
 */
static int hc_assemble_base_graph_merge_dangling_head_legacy(p_hc_apply apply)
{
    uint32_t* elements = apply->gatk_sw.output.cigar;
    uint32_t firstElement = elements[0];
    int indexesToMerge;
    p_assemble_vertex from, to;
    p_assemble_edge edge;

    if (bam_cigar_op(firstElement) != WORKER_SIGAR_STATUS_MATCH) {
        return 0;
    }

    indexesToMerge = hc_assemble_base_graph_best_prefix_match_legacy(apply, bam_cigar_oplen(firstElement));
    if (indexesToMerge <= 0) {
        return 0;
    }

    // we can't push back the reference path
    if (indexesToMerge >= (int)(apply->gatk_sw.input.ref_len - 1)) {
        return 0;
    }

    // but we can manipulate the dangling path if we need to
    if (indexesToMerge >= (int)apply->read_thread_graph->alt_size &&
        !hc_assemble_base_graph_extend_dangling_path_against_reference(apply, indexesToMerge - apply->read_thread_graph->alt_size + 2,
                                                                       elements)) {
        return 0;
    }

    from = hc_assemble_utils_get_list_index_item(apply->read_thread_graph->ref_path, indexesToMerge + 1);
    to = hc_assemble_utils_get_list_index_item(apply->read_thread_graph->alt_path, indexesToMerge);
    if (!from || !to) {
        return 0;
    }
    edge = assemble_graph_vertex_get_edge_between_vertex(from, to);
    if (__glibc_likely(!edge)) {
        assemble_graph_vertex_add_edge_to_vertex(apply->read_thread_graph, from, to, 1);
    }

    return 1;
}

/**
 * Determine whether the provided cigar is okay to merge into the reference path
 *
 * @param cigar                the cigar to analyze
 * @param requireFirstElementM if true, require that the first cigar element be an M operator in order for it to be okay
 * @param requireLastElementM  if true, require that the last cigar element be an M operator in order for it to be okay
 * @return true if it's okay to merge, false otherwise
 */
static int hc_assemble_base_graph_cigar_is_okay_to_merge(p_hc_apply apply, int requireFirstElementM, int requireLastElementM)
{
    uint32_t* elements = apply->gatk_sw.output.cigar;
    uint32_t numElements = apply->gatk_sw.output.cigar_len;

    // don't allow more than a couple of different ops
    if (numElements == 0 || numElements > HC_ASSEMBLE_BASE_GRAPH_MAX_CIGAR_COMPLEXITY) {
        return 0;
    }

    // the first element must be an M
    if (requireFirstElementM && bam_cigar_op(elements[0]) != WORKER_SIGAR_STATUS_MATCH) {
        return 0;
    }

    // the last element must be an M
    if (requireLastElementM && bam_cigar_op(elements[numElements - 1]) != WORKER_SIGAR_STATUS_MATCH) {
        return 0;
    }

    // note that there are checks for too many mismatches in the dangling branch later in the process

    return 1;
}

/**
 * Try to recover dangling tails
 *
 * @param pruneFactor             the prune factor to use in ignoring chain pieces if edge multiplicity is < pruneFactor
 * @param minDanglingBranchLength the minimum length of a dangling branch for us to try to merge it
 * @param recoverAll              recover even branches with forks
 * @param aligner
 */
void hc_assemble_base_graph_recover_dangling_tails(p_hc_apply apply)
{
    p_assemble_graph_hash_node vertex_node;
    p_assemble_vertex vertex;
    int attempted = 0;
    int nRecovered = 0;
    CDL_FOREACH(apply->read_thread_graph->vertics_hash->node_list, vertex_node)
    {
        vertex = vertex_node->value;
        if (vertex->out_degree == 0 && !hc_assemble_base_graph_is_ref_sink(vertex, apply->read_thread_graph)) {
            attempted++;
            nRecovered += hc_assemble_base_graph_recover_dangling_tail(vertex, apply);
        }
    }
}

/**
 * Try to recover dangling heads
 *
 * @param pruneFactor             the prune factor to use in ignoring chain pieces if edge multiplicity is < pruneFactor
 * @param minDanglingBranchLength the minimum length of a dangling branch for us to try to merge it
 * @param recoverAll              recover even branches with forks
 * @param aligner
 */
void hc_assemble_base_graph_recover_dangling_heads(p_hc_apply apply)
{
    p_assemble_graph_hash_node vertex_node;
    p_assemble_vertex vertex;

    CDL_FOREACH(apply->read_thread_graph->vertics_hash->node_list, vertex_node)
    {
        vertex = vertex_node->value;
        if (vertex->in_degree == 0 && !hc_assemble_base_graph_is_ref_source(vertex, apply->read_thread_graph)) {
            // singleton orphan
            if (vertex->out_degree == 0) {
                continue;
            }
            hc_assemble_base_graph_recover_dangling_head(vertex, apply);
        }
    }
}

/**
 * @param run_vetex the vertex to test
 * @return  true if this vertex is a reference sink
 */
int hc_assemble_base_graph_is_ref_sink(p_assemble_vertex vertex, p_assemble_graph graph)
{
    p_assemble_edge edge;

    // confirm that no outgoing edges are reference edges
    CDL_FOREACH2(vertex->out_edges, edge, link.from.next)
    {
        if (edge->is_ref) {
            return 0;
        }
    }
    // confirm that there is an incoming reference edge
    CDL_FOREACH2(vertex->from_edges, edge, link.to.next)
    {
        if (edge->is_ref) {
            return 1;
        }
    }
    // edge case: if the graph only has one node then it's a ref sink, otherwise it's not
    return graph->vertices_sum == 1;
}

/**
 * @param vertex the vertex to test
 * @return  true if this vertex is a reference source
 */
int hc_assemble_base_graph_is_ref_source(p_assemble_vertex vertex, p_assemble_graph graph)
{
    p_assemble_edge edge;

    // confirm that no incoming edges are reference edges
    CDL_FOREACH2(vertex->from_edges, edge, link.to.next)
    {
        if (edge->is_ref) {
            return 0;
        }
    }

    // confirm that there is an outgoing reference edge
    CDL_FOREACH2(vertex->out_edges, edge, link.from.next)
    {
        if (edge->is_ref) {
            return 1;
        }
    }

    // edge case: if the graph only has one node then it's a ref source, otherwise it's not
    return graph->vertices_sum == 1;
}