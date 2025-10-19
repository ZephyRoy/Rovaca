/**
 * @file hc_assemble_dijkstra_shortest_path.c
 * @author Cao Ce (caoce@genomics.cn)
 * @brief The implementation of Dijkstra algorithm to get the shortest path of a pair of vertices in a graph
 * @version 0.1
 * @date 2021-12-16
 *
 * @copyright Copyright (c) 2021 Rovaca SDK
 *
 */

#include "hc_assemble_dijkstra_shortest_path.h"

#include <math.h>

#include "debug.h"
#include "hc_apply.h"
#include "hc_func.h"
#include "hc_marco.h"
#include "rbtree.h"

#define HC_ASSEMBLE_SEQ_PATH_FINDER_MAX_NUMBER_OF_HAPLOTYPES (128)

static void hc_assemble_dijkstra_reset_tree(p_assemble_dijkstra_path graph);
static void hc_assemble_dijkstra_reset_ref(p_assemble_dijkstra_path graph, p_hc_shared_lib_list_head head);
static void hc_assemble_dijkstra_add_result(p_assemble_dijkstra_one_path_node add, p_assemble_dijkstra_path graph,
                                            p_hc_shared_lib_list_head out);
static void hc_assemble_dijkstra_init_insert_node(p_assemble_vertex initial_vertex, p_assemble_dijkstra_path graph, p_hc_apply apply);
static void hc_assemble_dijkstra_add_path(p_assemble_dijkstra_path graph, p_assemble_dijkstra_one_path_node path, p_assemble_edge edge,
                                          int total_outgoing_multiplicity);
static double hc_assemble_dijkstra_compute_log_penalty_score(int edge_multiplicity, int total_outgoing_multiplicity);
static int hc_assemble_dijkstra_insert_path(struct rb_root* root, p_assemble_dijkstra_one_path_node data, p_assemble_dijkstra_path graph);
static int hc_assemble_dijkstra_compare_edge_path(p_hc_shared_lib_list_item first, p_hc_shared_lib_list_item second,
                                                  p_assemble_dijkstra_path graph);
static inline int hc_assemble_dijkstra_compare_path_node(p_assemble_dijkstra_one_path_node first, p_assemble_dijkstra_one_path_node second,
                                                         p_assemble_dijkstra_path graph);

/**
 * Implement Dijkstra's algorithm as described in https://en.wikipedia.org/wiki/K_shortest_path_routing
 */
void hc_assemble_dijkstra_find_best_haplotypes(p_assemble_vertex source, p_assemble_vertex sink, p_hc_apply apply)
{
    p_assemble_dijkstra_path graph = &apply->read_graph.dijkstra_path_finder;
    p_assemble_dijkstra_one_path_node path_to_extend;
    p_assemble_edge edge;
    struct rb_root* queue = &graph->queue;
    struct rb_node* node;
    p_hc_shared_lib_list_head head = hc_assemble_utils_check_list_work(apply->read_graph.list_buffer_1);

    hc_assemble_dijkstra_init_insert_node(source, graph, apply);
    node = rb_first(queue);

    while (node && graph->result_size <= HC_ASSEMBLE_SEQ_PATH_FINDER_MAX_NUMBER_OF_HAPLOTYPES) {
        p_assemble_vertex vertex_to_extend;
        p_hc_shared_lib_list_item item, tmp1, tmp2;
        int total_outgoing_multiplicity = 0;

        path_to_extend = rb_entry(node, assemble_dijkstra_one_path_node, node);
        vertex_to_extend = path_to_extend->last_vertex;

        if (sink == vertex_to_extend) {
            hc_assemble_dijkstra_add_result(path_to_extend, graph, head);
            goto next_loop;
        }

        // 过多的call
        if (__glibc_unlikely(vertex_to_extend->call + 1 >= HC_ASSEMBLE_SEQ_PATH_FINDER_MAX_NUMBER_OF_HAPLOTYPES)) {
            mem_pool_fast_free(graph->path_pool, path_to_extend);
            goto next_loop;
        }
        vertex_to_extend->call += 1;

        CDL_FOREACH2(vertex_to_extend->out_edges, edge, link.from.next) { total_outgoing_multiplicity += edge->weight; }
        CDL_FOREACH2(vertex_to_extend->out_edges, edge, link.from.next)
        {
            hc_assemble_dijkstra_add_path(graph, path_to_extend, edge, total_outgoing_multiplicity);
        }
        CDL_FOREACH_SAFE(path_to_extend->edges_in_order, item, tmp1, tmp2)
        {
            CDL_DELETE(path_to_extend->edges_in_order, item);
            mem_pool_fast_free(graph->path_node_buffer, item);
        }
        mem_pool_fast_free(graph->path_pool, path_to_extend);
    next_loop:
        rb_erase(node, queue);
        node = rb_first(queue);
    }

    hc_assemble_dijkstra_reset_ref(graph, head);
    hc_assemble_dijkstra_reset_tree(graph);
}

static void hc_assemble_dijkstra_reset_tree(p_assemble_dijkstra_path graph)
{
    p_assemble_dijkstra_one_path_node path;
    p_hc_shared_lib_list_item item, tmp1, tmp2;
    struct rb_root* queue = &graph->queue;
    struct rb_node* node;

    while ((node = rb_first(queue))) {
        path = rb_entry(node, assemble_dijkstra_one_path_node, node);
        CDL_FOREACH_SAFE(path->edges_in_order, item, tmp1, tmp2)
        {
            CDL_DELETE(path->edges_in_order, item);
            mem_pool_fast_free(graph->path_node_buffer, item);
        }
        rb_erase(node, queue);
    }
    mem_pool_fast_reset(graph->path_node_buffer);
}

/**
 * @brief 重置ref质量
 *
 * @param graph     graph
 * @param head      list head
 */
static void hc_assemble_dijkstra_reset_ref(p_assemble_dijkstra_path graph, p_hc_shared_lib_list_head head)
{
    p_hc_shared_lib_list_item item, list_item;
    p_assemble_dijkstra_one_path_node search_item, iter;
    uint32_t all_size = 0;

    CDL_FOREACH(head->nodes, list_item)
    {
        uint32_t get_result = 0;
        search_item = list_item->data;
        CDL_FOREACH(graph->result->nodes, item)
        {
            iter = item->data;
            if (__glibc_unlikely(iter->st->seq_len == search_item->st->seq_len &&
                                 memcmp(search_item->st->seq, iter->st->seq, iter->st->seq_len) == 0)) {
                get_result = 1;
                if (search_item->is_reference) {
                    iter->score = search_item->score;
                }
                break;
            }
        }
        if (__glibc_likely(!get_result)) {
            hc_shared_list_insert(graph->result, search_item);
            all_size += 1;
        }
    }
    hc_shared_list_reset(head);
    graph->all_result_size += all_size;
}

static void hc_assemble_dijkstra_add_result(p_assemble_dijkstra_one_path_node add, p_assemble_dijkstra_path graph,
                                            p_hc_shared_lib_list_head out)
{
    p_hc_shared_lib_list_item item, tmp1, tmp2;
    p_assemble_edge edge;
    p_assemble_vertex vertex;

    add->st = ring_mempool_auto_enlarge_malloc(graph->path_out_st_pool);
    memset(add->st, 0, HC_ASSEMBLE_SW_ST_CLEAN_LEN);

    add->st->seq_len =
        sprintf((char*)add->st->seq, "%s%.*s", HC_ASSEMBLE_SW_PAD_STRING, graph->start_vertex->data_len, (char*)graph->start_vertex->data);

    CDL_FOREACH_SAFE(add->edges_in_order, item, tmp1, tmp2)
    {
        edge = item->data;
        vertex = edge->link.to.vertex;

        if (__glibc_unlikely(!vertex)) {
            return;
        }
        if (__glibc_unlikely(add->st->seq_len + vertex->data_len + strlen(HC_ASSEMBLE_SW_PAD_STRING) >= ASSEMBLE_GRAPH_MAX_SEQ_LEN)) {
            add->st->seq_len = ASSEMBLE_GRAPH_MAX_SEQ_LEN;
            return;
        }

        add->st->seq_len += sprintf((char*)add->st->seq + add->st->seq_len, "%.*s", vertex->data_len, (char*)vertex->data);
        CDL_DELETE(add->edges_in_order, item);
        mem_pool_fast_free(graph->path_node_buffer, item);
    }
    add->st->seq_len += sprintf((char*)add->st->seq + add->st->seq_len, "%s", HC_ASSEMBLE_SW_PAD_STRING);

    add->done = 0;
    hc_shared_list_insert(out, add);
    graph->result_size += 1;
}

/**
 * Create a new Path containing no edges and starting at initial_vertex
 * @param initial_vertex the starting vertex of the path
 * @param graph the graph this path will follow through
 */
static void hc_assemble_dijkstra_init_insert_node(p_assemble_vertex initial_vertex, p_assemble_dijkstra_path graph, p_hc_apply apply)
{
    p_assemble_dijkstra_one_path_node node;
    p_assemble_graph_hash_node vertex_node;
    p_assemble_vertex vertex;

    // reset
    graph->start_vertex = initial_vertex;
    graph->queue = RB_ROOT;
    graph->cached_seq = apply->read_thread_graph->seq;
    graph->cached_ref = apply->read_thread_graph->ref;
    graph->result_size = 0;

    // init node
    node = mem_pool_fast_alloc(graph->path_pool);
    memset(node, 0, sizeof(assemble_dijkstra_one_path_node));
    node->last_vertex = initial_vertex;
    hc_assemble_dijkstra_insert_path(&graph->queue, node, graph);

    // graph node
    CDL_FOREACH(apply->read_thread_graph->vertics_hash->node_list, vertex_node)
    {
        vertex = vertex_node->value;
        vertex->call = 0;
    }

    if (graph->result->nodes) {
        return;
    }

    // ref
    node = mem_pool_fast_alloc(graph->path_pool);
    memset(node, 0, sizeof(assemble_dijkstra_one_path_node));
    node->is_reference = 1;
    node->done = 0;
    node->score = 0;
    node->st = ring_mempool_auto_enlarge_malloc(graph->path_out_st_pool);
    memset(node->st, 0, HC_ASSEMBLE_SW_ST_CLEAN_LEN);

    node->st->seq_len = sprintf((char*)node->st->seq, "%s%.*s%s", HC_ASSEMBLE_SW_PAD_STRING, apply->read_thread_graph->origin_ref_len,
                                (char*)apply->read_thread_graph->origin_ref, HC_ASSEMBLE_SW_PAD_STRING);
    hc_shared_list_insert(graph->result, node);
    graph->result_size += 1;
}

/**
 * @brief reset dijkstra result
 *
 * @param graph dijkstra graph
 */
void hc_assemble_dijkstra_reset_result(p_assemble_dijkstra_path graph)
{
    p_assemble_dijkstra_one_path_node hap_item;
    p_hc_shared_lib_list_item item;

    CDL_FOREACH(graph->result->nodes, item)
    {
        hap_item = item->data;
        if (hap_item->st) {
            memset(hap_item->st->cigar, 0, sizeof(uint32_t) * hap_item->st->cigar_len);
            memset(hap_item->st->cigar_back, 0, sizeof(uint32_t) * hap_item->st->cigar_back_len);
        }
        mem_pool_fast_free(graph->path_pool, hap_item);
    }
    hc_assemble_utils_check_list_work(graph->result);
    ring_mempool_auto_enlarge_reset(graph->path_out_st_pool);
    graph->all_result_size = 0;
}

/**
 * Create a new Path extending p with edge
 *
 * @param p the path to extend.
 * @param edge the edge to extend path with.
 *
 * @throws IllegalArgumentException if {@code p} or {@code edge} are {@code null}, or {@code edge} is
 * not part of {@code p}'s graph, or {@code edge} does not have as a source the last vertex in {@code p}.
 */
static void hc_assemble_dijkstra_add_path(p_assemble_dijkstra_path graph, p_assemble_dijkstra_one_path_node path, p_assemble_edge edge,
                                          int total_outgoing_multiplicity)
{
    p_assemble_dijkstra_one_path_node node;
    p_hc_shared_lib_list_item item, item_to_add, tmp1, tmp2;
    p_assemble_edge edge_before;

    node = mem_pool_fast_alloc(graph->path_pool);
    memset(node, 0, sizeof(assemble_dijkstra_one_path_node));
    node->last_vertex = edge->link.to.vertex;
    node->is_reference = 1;

    CDL_FOREACH_SAFE(path->edges_in_order, item, tmp1, tmp2)
    {
        edge_before = item->data;

        item_to_add = mem_pool_fast_alloc(graph->path_node_buffer);
        item_to_add->data = edge_before;
        CDL_APPEND(node->edges_in_order, item_to_add);
        node->is_reference &= edge_before->is_ref;
    }
    item_to_add = mem_pool_fast_alloc(graph->path_node_buffer);
    item_to_add->data = edge;
    CDL_APPEND(node->edges_in_order, item_to_add);

    node->score = path->score + hc_assemble_dijkstra_compute_log_penalty_score(edge->weight, total_outgoing_multiplicity);
    node->is_reference &= edge->is_ref;

    hc_assemble_dijkstra_insert_path(&graph->queue, node, graph);
}

static double hc_assemble_dijkstra_compute_log_penalty_score(int edge_multiplicity, int total_outgoing_multiplicity)
{
    return log10(edge_multiplicity) - log10(total_outgoing_multiplicity);
}

static int hc_assemble_dijkstra_compare_edge_path(p_hc_shared_lib_list_item first, p_hc_shared_lib_list_item second,
                                                  p_assemble_dijkstra_path graph)
{
    p_hc_shared_lib_list_item nodes;
    uint8_t *first_seq = graph->cached_seq, *second_seq = graph->cached_ref;
    uint32_t first_seq_len = 0, second_seq_len = 0, min_len, i;
    p_assemble_vertex vertex;
    p_assemble_edge edge;

    // first
    if (first) {
        edge = first->data;
        vertex = edge->link.from.vertex;
        memcpy(first_seq, vertex->data, vertex->data_len);
        first_seq_len += vertex->data_len;
        CDL_FOREACH(first, nodes)
        {
            edge = nodes->data;
            vertex = edge->link.to.vertex;
            memcpy(first_seq + first_seq_len, vertex->data, vertex->data_len);
            first_seq_len += vertex->data_len;
        }
    }

    // second
    if (second) {
        edge = second->data;
        vertex = edge->link.from.vertex;
        memcpy(second_seq, vertex->data, vertex->data_len);
        second_seq_len += vertex->data_len;
        CDL_FOREACH(second, nodes)
        {
            edge = nodes->data;
            vertex = edge->link.to.vertex;
            memcpy(second_seq + second_seq_len, vertex->data, vertex->data_len);
            second_seq_len += vertex->data_len;
        }
    }
    for (min_len = assemble_min(first_seq_len, second_seq_len), i = 0; i < min_len; i++) {
        int ret = first_seq[i] - second_seq[i];
        if (ret == 0) {
            continue;
        }
        else if (ret > 0) {
            return -1;
        }
        else {
            return 1;
        }
    }

    if (first_seq_len > second_seq_len) {
        return -1;
    }
    else if (first_seq_len < second_seq_len) {
        return 1;
    }
    else {
        return 0;
    }
}

/**
 * @brief 比较2个路径点大小
 *
 * @param first     第一个
 * @param second    第二个
 *
 * @return int      >0:第一个大; <0:第二个大; =0:相等
 */
static inline int hc_assemble_dijkstra_compare_path_node(p_assemble_dijkstra_one_path_node first, p_assemble_dijkstra_one_path_node second,
                                                         p_assemble_dijkstra_path graph)
{
    if (first->score > second->score) {
        return -1;
    }
    else if (first->score < second->score) {
        return 1;
    }

    return hc_assemble_dijkstra_compare_edge_path(first->edges_in_order, second->edges_in_order, graph);
}

static int hc_assemble_dijkstra_insert_path(struct rb_root* root, p_assemble_dijkstra_one_path_node data, p_assemble_dijkstra_path graph)
{
    struct rb_node** new = &(root->rb_node), *parent = NULL;

    /* Figure out where to put new node */
    while (*new) {
        p_assemble_dijkstra_one_path_node this = container_of(*new, assemble_dijkstra_one_path_node, node);
        int result = hc_assemble_dijkstra_compare_path_node(data, this, graph);

        parent = *new;
        if (result < 0)
            new = &((*new)->rb_left);
        else if (result > 0)
            new = &((*new)->rb_right);
        else
            new = &((*new)->rb_left);
    }

    rb_link_node(&data->node, parent, new);
    rb_insert_color(&data->node, root);

    return 1;
}