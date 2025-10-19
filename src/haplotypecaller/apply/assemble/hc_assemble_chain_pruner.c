/**
 * @file hc_assemble_chain_pruner.c
 * @author Cao Ce (caoce@genomics.cn)
 * @brief Chain Pruner
 * @version 0.1
 * @date 2022-04-22
 *
 * @copyright Copyright (c) 2022 ROVACA SDK
 *
 */
#include "graph.h"
#include "hc_apply.h"
#include "hc_assemble.h"
#include "hc_func.h"

static void hc_assemble_chain_pruner_find_all_chains(p_assemble_vertex start, p_hc_shared_lib_list_head chain_stars,
                                                     p_hc_shared_lib_list_head all_path);
static p_assemble_vertex hc_assemble_chain_pruner_find_chain(p_assemble_edge start_edge, p_hc_shared_lib_list_head all_path);
static inline void hc_assemble_chain_pruner_insert_edge_to_path(p_hc_shared_lib_list_head all_path, p_assemble_edge edge,
                                                                p_hc_shared_lib_list_item path);
static int hc_assemble_chain_pruner_detect_to_remove(p_hc_shared_lib_list_item all_edge, p_assemble_graph graph);

/**
 * @brief 使用chain做裁剪
 *
 * @param apply apply
 */
void hc_assemble_chain_pruner_prune_low_weight_chains(p_hc_apply apply)
{
    p_assemble_graph_hash_node vertex_node;
    p_hc_shared_lib_list_head all_path = hc_assemble_utils_check_list_work(apply->read_graph.list_buffer_1);
    p_hc_shared_lib_list_head chain_starts = hc_assemble_utils_check_list_work(apply->read_graph.list_buffer_2);
    p_hc_shared_lib_list_item one_path;
    p_hc_shared_lib_list_item chain_start;
    p_assemble_vertex vertex;

    CDL_FOREACH(apply->read_thread_graph->vertics_hash->node_list, vertex_node)
    {
        vertex = vertex_node->value;
        if (__glibc_likely(vertex->in_degree || !vertex->out_degree)) {
            continue;
        }
        hc_shared_list_insert(chain_starts, vertex);
    }

    CDL_FOREACH(chain_starts->nodes, chain_start)
    {
        vertex = chain_start->data;
        hc_assemble_chain_pruner_find_all_chains(vertex, chain_starts, all_path);
    }

    CDL_FOREACH(all_path->nodes, one_path)
    {
        p_hc_shared_lib_list_item path_st = one_path->data;

        hc_assemble_chain_pruner_detect_to_remove(path_st, apply->read_thread_graph);
    }
    hc_shared_list_reset(all_path);
    hc_shared_list_reset(chain_starts);
    hc_assemble_seq_graph_remove_singleton_orphan_vertices(apply);
}

/**
 * @brief 单次路径遍历裁剪
 *
 * @param start start
 * @param apply apply
 */
static void hc_assemble_chain_pruner_find_all_chains(p_assemble_vertex start, p_hc_shared_lib_list_head chain_starts,
                                                     p_hc_shared_lib_list_head all_path)
{
    p_hc_shared_lib_list_item edges = NULL, one_edge;
    p_assemble_edge edge;
    p_assemble_vertex chain_end;

    CDL_FOREACH2(start->out_edges, edge, link.from.next)
    {
        if (__glibc_unlikely(!edge->link.from.vertex || !edge->link.to.vertex)) {
            continue;
        }
        one_edge = ring_mempool_auto_enlarge_malloc(chain_starts->pool);
        one_edge->data = edge;
        CDL_APPEND(edges, one_edge);
    }

    CDL_FOREACH(edges, one_edge)
    {
        edge = one_edge->data;
        chain_end = hc_assemble_chain_pruner_find_chain(edge, all_path);

        if (chain_end && !hc_shared_list_search(chain_starts, chain_end)) {
            hc_shared_list_insert(chain_starts, chain_end);
        }
    }
}

static inline void hc_assemble_chain_pruner_insert_edge_to_path(p_hc_shared_lib_list_head all_path, p_assemble_edge edge,
                                                                p_hc_shared_lib_list_item path)
{
    p_hc_shared_lib_list_item item = ring_mempool_auto_enlarge_malloc(all_path->pool);

    item->data = edge;
    CDL_APPEND(path, item);
}

static p_assemble_vertex hc_assemble_chain_pruner_find_chain(p_assemble_edge start_edge, p_hc_shared_lib_list_head all_path)
{
    p_assemble_vertex first_vertex = start_edge->link.from.vertex;
    p_assemble_vertex last_vertex = start_edge->link.to.vertex;
    p_hc_shared_lib_list_item path = ring_mempool_auto_enlarge_malloc(all_path->pool);

    path->prev = path;
    path->next = path;
    path->data = start_edge;

    while (1) {
        p_assemble_edge out_edge;

        if (last_vertex->out_degree != 1 || last_vertex->in_degree > 1 || last_vertex == first_vertex) {
            break;
        }
        out_edge = last_vertex->out_edges;
        hc_assemble_chain_pruner_insert_edge_to_path(all_path, out_edge, path);
        last_vertex = out_edge->link.to.vertex;
    }
    hc_shared_list_insert(all_path, path);
    return last_vertex;
}

/**
 * @brief 修剪该图中的所有链，其中路径中的所有边都具有多重性 < HC_ASSEMBLE_GRAPH_CHAIN_PRUNE_FACTOR
 *
 * @param chain_start   chain start
 * @param chain_end     chain end
 * @param iter          iter
 *
 * @return int          succeed
 *
 * @note
 *      For A -[1]> B -[1]> C -[1]> D would be removed with pruneFactor 2
 *      but A -[1]> B -[2]> C -[1]> D would not be because the linear chain includes an edge with weight >= 2
 */
static int hc_assemble_chain_pruner_detect_to_remove(p_hc_shared_lib_list_item all_edge, p_assemble_graph graph)
{
    p_hc_shared_lib_list_item node;

    CDL_FOREACH(all_edge, node)
    {
        p_assemble_edge edge = node->data;

        if (edge->weight >= HC_ASSEMBLE_GRAPH_CHAIN_PRUNE_FACTOR || edge->is_ref) {
            return 0;
        }
    }
    CDL_FOREACH(all_edge, node)
    {
        p_assemble_edge edge = node->data;
        assemble_graph_remove_edge(graph, edge);
    }
    return 1;
}