/**
 * @file hc_assemble_base_graph_iter.c
 * @author Cao Ce (caoce@genomics.cn)
 * @brief Graph Iter
 * @version 0.1
 * @date 2021-11-23
 *
 * @copyright Copyright (c) 2021 ROVACA SDK
 *
 */

#include "graph.h"
#include "hc_assemble.h"
#include "hc_assemble_graph.h"
#include "hc_func.h"
#include "list.h"
#define HC_ASSEMBLE_BASE_GRAPH_ITER_BACK_NODE (8192)
#define HC_ASSEMBLE_BASE_GRAPH_PAIR_LIST_NODE (2048)

/**
 * @brief insert prepend
 *
 * @param head  head
 * @param data  data
 *
 * @return int  1:succeed, 0:failed
 */
static int hc_assemble_base_graph_iter_list_insert_prepend(p_assemble_graph_iter head, void* data)
{
    p_hc_shared_lib_list_item one_node;

    one_node = mem_pool_fast_alloc(head->run_mem);
    if (__glibc_unlikely(!one_node)) {
        return 0;
    }

    one_node->data = data;
    CDL_PREPEND(head->run_path, one_node);

    return 1;
}

/**
 * @brief insert prepend
 *
 * @param head  head
 * @param data  data
 *
 * @return int  1:succeed, 0:failed
 */
static int hc_assemble_base_graph_iter_pair_list_insert_prepend(p_assemble_graph_iter iter, void* data1, void* data2)
{
    p_hc_shared_lib_pair_list_item one_node;

    one_node = mem_pool_fast_alloc(iter->node_mem);
    if (__glibc_unlikely(!one_node)) {
        return 0;
    }

    one_node->first_data = data1;
    one_node->second_data = data2;
    CDL_PREPEND(iter->node_stack, one_node);

    return 1;
}

/**
 * @brief 重置此后的vertex内
 *
 * @param edge_head
 */
static inline void hc_assemble_base_graph_iter_reset_visit(p_assemble_graph_iter path, p_assemble_vertex vertex)
{
    p_hc_shared_lib_list_item one_node, tmp1, tmp2;

    CDL_FOREACH_SAFE(path->run_path, one_node, tmp1, tmp2)
    {
        if (one_node->data == vertex) {
            return;
        }
        ((p_assemble_vertex)(one_node->data))->visited = 0;
        CDL_DELETE(path->run_path, one_node);
        mem_pool_fast_free(path->run_mem, one_node);
    }
}

/**
 * @brief vertex list pop
 *
 * @param head list
 *
 * @return void* pop
 */
static inline void hc_assemble_base_graph_iter_list_pop(p_assemble_graph_iter iter, p_assemble_vertex* path_from,
                                                        p_assemble_vertex* path_head)
{
    p_hc_shared_lib_pair_list_item one_node, tmp1, tmp2;

    CDL_FOREACH_SAFE(iter->node_stack, one_node, tmp1, tmp2)
    {
        *path_from = one_node->first_data;
        *path_head = one_node->second_data;
        CDL_DELETE(iter->node_stack, one_node);
        mem_pool_fast_free(iter->node_mem, one_node);
        return;
    }
}

/**
 * @brief vertex list reset
 *
 * @param head list
 *
 * @return void* pop
 */
static inline void hc_assemble_base_graph_iter_list_reset(p_assemble_graph_iter iter)
{
    p_hc_shared_lib_pair_list_item one_node, tmp1, tmp2;

    CDL_FOREACH_SAFE(iter->node_stack, one_node, tmp1, tmp2)
    {
        CDL_DELETE(iter->node_stack, one_node);
        mem_pool_fast_free(iter->node_mem, one_node);
    }
    mem_pool_fast_reset(iter->node_mem);
    mem_pool_fast_reset(iter->run_mem);
}

/**
 * @brief 初始化iter资源
 *
 * @param[in] input     iter
 *
 * @return int 1:succeed, 0:failed
 */
int hc_assemble_base_graph_iter_init(p_assemble_graph_iter input)
{
    input->run_mem = mem_pool_fast_init(sizeof(hc_shared_lib_list_item), HC_ASSEMBLE_BASE_GRAPH_ITER_BACK_NODE);
    input->node_mem = mem_pool_fast_init(sizeof(hc_shared_lib_pair_list_item), HC_ASSEMBLE_BASE_GRAPH_PAIR_LIST_NODE);
    input->safe_nodes = hc_shared_list_init(HC_ASSEMBLE_GRAPH_POOL_PER_REGION);

    if (input->run_mem && input->node_mem && input->safe_nodes) {
        return 1;
    }
    else {
        return 0;
    }
}

/**
 * @brief 释放iter资源
 *
 * @param[in] input     iter
 */
void hc_assemble_base_graph_iter_finit(p_assemble_graph_iter input)
{
    mem_pool_fast_term(input->run_mem);
    mem_pool_fast_term(input->node_mem);
    hc_shared_list_term(input->safe_nodes);
}

/**
 * @brief 深度优先遍历，无迭代
 *
 * @param from  from
 * @param to    to
 * @param path  path
 *
 * @return int  cycle
 */
int hc_assemble_base_graph_iter_dfs_non_recursion(p_assemble_graph_iter input)
{
    p_assemble_vertex to = input->input_layer.to;
    p_assemble_edge edge;
    p_assemble_vertex vertex, path_from_vertex, path_from_vertex_next;
    uint32_t iter_succeed;

    // 无意义图
    if (__glibc_unlikely(!input->input_layer.from)) {
        return 0;
    }

    // 初始化
    vertex = input->input_layer.from;
    path_from_vertex = NULL;
    input->run_path = NULL;
    input->node_stack = NULL;
    path_from_vertex_next = NULL;
    hc_assemble_utils_check_list_work(input->safe_nodes);

    while (1) {
        // Error Graph / Cycle
        if (__glibc_unlikely(!vertex || vertex->visited)) {
            int ret = input->cycle_detected(vertex, input);
            hc_assemble_base_graph_iter_reset_visit(input, NULL);
            hc_assemble_base_graph_iter_list_reset(input);
            return ret;
        }

        hc_assemble_base_graph_iter_list_insert_prepend(input, vertex);
        vertex->visited = 1;
        if (__glibc_unlikely(vertex == to)) {
            goto end_point;
        }
        switch (vertex->out_degree) {
            case 1:  // 常规路径
                vertex = vertex->out_edges->link.to.vertex;
                break;
            case 0:  // 结束点
            end_point:
                if (__glibc_unlikely(!input->node_stack)) {
                    if (input->one_path_done) {
                        input->one_path_done(vertex, input);
                    }
                    if (input->all_path_done) {
                        input->all_path_done(vertex, input);
                    }
                    // 图遍历结束
                    hc_assemble_base_graph_iter_reset_visit(input, NULL);
                    return 0;
                }
                if (input->one_path_done) {
                    input->one_path_done(vertex, input);
                }

                if (path_from_vertex_next && !hc_shared_list_search(input->safe_nodes, path_from_vertex_next)) {
                    hc_shared_list_insert_prev(input->safe_nodes, path_from_vertex_next);
                }

                hc_assemble_base_graph_iter_list_pop(input, &path_from_vertex, &vertex);
                path_from_vertex_next = vertex;
                hc_assemble_base_graph_iter_reset_visit(input, path_from_vertex);
                break;
            default:  // 分歧点
                iter_succeed = 0;

                path_from_vertex = vertex;
                CDL_FOREACH2(vertex->out_edges, edge, link.from.next)
                {
                    if (edge->link.to.vertex && hc_shared_list_search(input->safe_nodes, edge->link.to.vertex)) {
                        continue;
                    }
                    iter_succeed = 1;
                    if (edge->link.to.vertex && edge != vertex->out_edges) {
                        hc_assemble_base_graph_iter_pair_list_insert_prepend(input, vertex, edge->link.to.vertex);
                    }
                }
                if (__glibc_unlikely(!iter_succeed)) {
                    goto end_point;
                }
                vertex = vertex->out_edges->link.to.vertex;
                path_from_vertex_next = vertex;
                break;
        }
    }

    hc_assemble_base_graph_iter_reset_visit(input, NULL);
    hc_assemble_base_graph_iter_list_reset(input);
    return 0;
}

/**
 * @brief 深度优先遍历，无迭代，向上
 *
 * @param from  from
 * @param to    to
 * @param path  path
 *
 * @return int  cycle
 */
int hc_assemble_base_graph_iter_dfs_non_recursion_up(p_assemble_graph_iter input)
{
    p_assemble_vertex to = input->input_layer.to;
    p_assemble_edge edge;
    p_assemble_vertex vertex, path_from_vertex, path_from_vertex_next;
    uint32_t iter_succeed;

    // 无意义图
    if (__glibc_unlikely(!input->input_layer.from)) {
        return 0;
    }

    // 初始化
    vertex = input->input_layer.from;
    input->run_path = NULL;
    input->node_stack = NULL;
    path_from_vertex = NULL;
    path_from_vertex_next = NULL;
    hc_assemble_utils_check_list_work(input->safe_nodes);

    while (1) {
        // Error Graph / Cycle
        if (__glibc_unlikely(!vertex || vertex->visited)) {
            int ret = input->cycle_detected(vertex, input);
            hc_assemble_base_graph_iter_reset_visit(input, NULL);
            hc_assemble_base_graph_iter_list_reset(input);
            return ret;
        }

        hc_assemble_base_graph_iter_list_insert_prepend(input, vertex);
        vertex->visited = 1;
        if (__glibc_unlikely(vertex == to)) {
            goto end_point;
        }
        switch (vertex->in_degree) {
            case 1:  // 常规路径
                vertex = vertex->from_edges->link.from.vertex;
                break;
            case 0:  // 结束点
            end_point:
                if (__glibc_unlikely(!input->node_stack)) {
                    if (input->one_path_done) {
                        input->one_path_done(vertex, input);
                    }
                    if (input->all_path_done) {
                        input->all_path_done(vertex, input);
                    }
                    // 图遍历结束
                    hc_assemble_base_graph_iter_reset_visit(input, NULL);
                    return 0;
                }
                if (input->one_path_done) {
                    input->one_path_done(vertex, input);
                }
                if (path_from_vertex_next && !hc_shared_list_search(input->safe_nodes, path_from_vertex_next)) {
                    hc_shared_list_insert_prev(input->safe_nodes, path_from_vertex_next);
                }

                hc_assemble_base_graph_iter_list_pop(input, &path_from_vertex, &vertex);
                path_from_vertex_next = vertex;
                hc_assemble_base_graph_iter_reset_visit(input, path_from_vertex);
                break;
            default:  // 分歧点
                iter_succeed = 0;
                path_from_vertex = vertex;
                CDL_FOREACH2(vertex->from_edges, edge, link.to.next)
                {
                    if (edge->link.from.vertex && hc_shared_list_search(input->safe_nodes, edge->link.from.vertex)) {
                        continue;
                    }
                    iter_succeed = 1;

                    if (edge->link.from.vertex && edge != vertex->from_edges) {
                        hc_assemble_base_graph_iter_pair_list_insert_prepend(input, vertex, edge->link.from.vertex);
                    }
                }
                // 结束点
                if (__glibc_unlikely(!iter_succeed)) {
                    goto end_point;
                }

                vertex = vertex->from_edges->link.from.vertex;
                path_from_vertex_next = vertex;
                break;
        }
    }

    hc_assemble_base_graph_iter_reset_visit(input, NULL);
    hc_assemble_base_graph_iter_list_reset(input);
    return 0;
}
