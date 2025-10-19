//
//  graph.c
//  adjacencylist
//
//  Created by Sam Goldman on 6/21/11.
//  Copyright 2011 Sam Goldman. All rights reserved.
//

// #include "hc_main.h"
#include "graph.h"

#include <stdlib.h>

#include "hash.h"
#include "hc_assemble_graph.h"
#include "hc_func.h"
#include "hc_marco.h"
#include "utlist.h"

static int hc_assemble_utils_working_graph_vertex_equal_func(void* old, void* new);
static int assemble_graph_hash_equal_func(void* old, void* new);
static void assemble_graph_hash_new_node_memcpy(void* new_node, void* input_pt, int hash_len, p_assemble_graph_hash_table hash_table);
static p_assemble_edge assemble_graph_edge_create(p_assemble_graph graph, p_assemble_vertex from_vertex, p_assemble_vertex to_vertex,
                                                  double weight);
static void assemble_graph_vertex_new_node_memcpy(void* new_node, void* vertex_pt, int hash_len, p_assemble_graph_hash_table hash_table);
static inline void* assemble_graph_add_read_threading_vertex(p_assemble_graph graph, p_assemble_vertex vertex);
static inline void* assemble_graph_add_seq_vertex(p_assemble_graph graph, p_assemble_vertex vertex);

/**
 * @brief 插入相同vertex时Hash Table执行函数
 *
 * @param[in] old   old data
 * @param[in] new   new data
 */
static int hc_assemble_utils_working_graph_vertex_equal_func(void* old, void* new)
{
    (void)old;
    (void)new;
    return 0;
}

/**
 * @brief 图转换，转换为常规Read Thread Graph
 *
 * @param[in] graph 图
 */
void assemble_graph_to_normal_read_thread_graph(p_assemble_graph graph)
{
    graph->graph_type = ASSEMBLE_GRAPH_TYPE_READ_THREADING_GRAPH;
    graph->vertics_hash->equal_func = assemble_graph_vertex_equal_func;
}

/**
 * @brief 图转换，转换为活动Read Thread Graph
 *
 * @param[in] graph 图
 */
void assemble_graph_to_working_read_thread_graph(p_assemble_graph graph)
{
    graph->graph_type = ASSEMBLE_GRAPH_TYPE_READ_THREADING_GRAPH;
    graph->vertics_hash->equal_func = hc_assemble_utils_working_graph_vertex_equal_func;
}

/**
 * @brief 图转换，转换为Seq Graph
 *
 * @param[in] graph 图
 */
void assemble_graph_to_seq_graph(p_assemble_graph graph)
{
    graph->graph_type = ASSEMBLE_GRAPH_TYPE_SEQ_GRAPH;
    graph->vertics_hash->equal_func = (void*)1;
}

/**
 * @brief 在vertex上新建一个edge
 *
 * @param graph     graph
 * @param vertex    vertex
 * @param weight    权
 *
 * @return p_assemble_edge 新edge
 */
static p_assemble_edge assemble_graph_edge_create(p_assemble_graph graph, p_assemble_vertex from_vertex, p_assemble_vertex to_vertex,
                                                  double weight)
{
    p_assemble_edge edge;

    edge = ring_mempool_auto_enlarge_malloc(graph->edge_pool);
    if (__glibc_unlikely(!edge)) {
        return NULL;
    }
    memset(edge, 0, sizeof(assemble_edge));

    edge->weight = weight;
    edge->link.from.vertex = from_vertex;
    edge->link.to.vertex = to_vertex;

    return edge;
}

/**
 * @brief 插入相同vertex时Hash Table执行函数
 *
 * @param old   old data
 * @param new   new data
 */
int assemble_graph_vertex_equal_func(void* old, void* new)
{
    (void)old;
    p_assemble_vertex new_vertex = new;
    p_hc_assemble_graph_kmer new_kmer = new_vertex->storage;

    if (!new_kmer->is_ref) {
        // old_vertex->weight += 1; todo
        return 1;
    }

    return 0;
}

/**
 * @brief Seq Graph插入相同vertex时Hash Table执行函数
 *
 * @param old   old data
 * @param new   new data
 * @return int  0
 */
static int assemble_graph_hash_equal_func(void* old, void* new)
{
    (void)old;
    (void)new;
    return 1;
}

/**
 * @brief 插入vertex时，复制vertex内容的Hash Table执行函数
 *
 * @param new_node      新Hash Nnode
 * @param vertex_pt     插入Vertex
 * @param hash_len      Hash Len
 * @param hash_table    Hash Table
 *
 */
static void assemble_graph_hash_new_node_memcpy(void* new_node, void* input_pt, int hash_len, p_assemble_graph_hash_table hash_table)
{
    (void)hash_table;
    p_assemble_graph_hash_node insert_node = new_node;

    memset(insert_node, 0, sizeof(assemble_graph_hash_node));

    insert_node->key = input_pt;
    insert_node->keylen = hash_len;
    insert_node->value = input_pt;
}

/**
 * @brief 插入vertex时，复制vertex内容的Hash Table执行函数
 *
 * @param new_node      新Hash Nnode
 * @param vertex_pt     插入Vertex
 * @param hash_len      Hash Len
 * @param hash_table    Hash Table
 *
 */
static void assemble_graph_vertex_new_node_memcpy(void* new_node, void* vertex_pt, int hash_len, p_assemble_graph_hash_table hash_table)
{
    (void)hash_len;
    (void)hash_table;
    p_assemble_graph_hash_node insert_node = new_node;
    p_assemble_vertex insert_vertex = vertex_pt;
    uint32_t index = insert_node->index;

    memset(insert_node, 0, sizeof(assemble_graph_hash_node));

    insert_node->key = insert_vertex->data;
    insert_node->keylen = insert_vertex->data_len;
    insert_node->value = insert_vertex;
    insert_vertex->key = index;
}

/**
 * @brief 新建一个graph
 *
 * @return p_assemble_graph 新graph
 */
p_assemble_graph assemble_graph_create(void)
{
    p_assemble_graph graph = NULL;

    graph = malloc(sizeof(assemble_graph));
    if (!graph) {
        return NULL;
    }
    memset(graph, 0, sizeof(assemble_graph));

    if (!(graph->vertex_pool = ring_mempool_auto_enlarge_init(sizeof(assemble_vertex), ASSEMBLE_GRAPH_VERTEX_INIT_NUM)) ||
        !(graph->edge_pool = ring_mempool_auto_enlarge_init(sizeof(assemble_edge), ASSEMBLE_GRAPH_EDGE_INIT_NUM)) ||
        !(graph->alt_path = hc_shared_list_init(ASSEMBLE_GRAPH_EDGE_INIT_NUM)) ||
        !(graph->ref_path = hc_shared_list_init(ASSEMBLE_GRAPH_EDGE_INIT_NUM)) ||
        !(graph->vertics_hash = assemble_graph_hash_init(ASSEMBLE_GRAPH_HASH_NODES, assemble_graph_vertex_equal_func,
                                                         assemble_graph_vertex_new_node_memcpy)) ||
        !(graph->non_unique_kmers =
              assemble_graph_hash_init(ASSEMBLE_GRAPH_HASH_NODES, assemble_graph_hash_equal_func, assemble_graph_hash_new_node_memcpy)) ||
        !(graph->simplify_graph.edge_pool = ring_mempool_auto_enlarge_init(sizeof(assemble_edge), ASSEMBLE_GRAPH_EDGE_INIT_NUM)) ||
        !(graph->simplify_graph.vertex = hc_shared_list_init(ASSEMBLE_GRAPH_EDGE_INIT_NUM)) ||
        !(graph->simplify_graph.edge = hc_shared_list_init(ASSEMBLE_GRAPH_EDGE_INIT_NUM)) ||
        !hc_assemble_base_graph_iter_init(&graph->graph_iter)) {
        if (graph->vertex_pool) {
            ring_mempool_auto_enlarge_term(graph->vertex_pool);
        }
        if (graph->edge_pool) {
            ring_mempool_auto_enlarge_term(graph->edge_pool);
        }
        if (graph->alt_path) {
            hc_shared_list_term(graph->alt_path);
        }
        if (graph->ref_path) {
            hc_shared_list_term(graph->ref_path);
        }
        free(graph);
        return NULL;
    }

    return graph;
}

/**
 * @brief 依据data，新建一个vertex
 *
 * @param graph     graph
 * @param data      data
 *
 * @return p_assemble_vertex    新vertex
 */
p_assemble_vertex assemble_graph_vertex_create(p_assemble_graph graph, void* data, uint32_t data_len, void* storage)
{
    p_assemble_vertex vertex;

    vertex = ring_mempool_auto_enlarge_malloc(graph->vertex_pool);
    if (__glibc_unlikely(!vertex)) {
        return NULL;
    }

    vertex->data = data;
    vertex->data_len = data_len;
    vertex->storage = storage;
    vertex->vertex_len = 1;
    vertex->from_edges = NULL;
    vertex->out_edges = NULL;
    vertex->in_degree = 0;
    vertex->out_degree = 0;
    vertex->visited = 0;
    vertex->call = 0;
    vertex->is_dup = 0;

    return vertex;
}

/**
 * @brief graph添加已存在的vertex
 *
 * @param graph     graph
 * @param vertex    vertex
 */
void* assemble_graph_add_vertex(p_assemble_graph graph, p_assemble_vertex vertex)
{
    if (graph->graph_type == ASSEMBLE_GRAPH_TYPE_READ_THREADING_GRAPH) {
        return assemble_graph_add_read_threading_vertex(graph, vertex);
    }
    else {
        return assemble_graph_add_seq_vertex(graph, vertex);
    }
}

static inline void* assemble_graph_add_read_threading_vertex(p_assemble_graph graph, p_assemble_vertex vertex)
{
    void* res;
    uint32_t key;
    p_assemble_graph_hash_node vertex_hash;

    key = assemble_graph_hash_get_key(graph->non_unique_kmers, vertex->data, vertex->data_len);
    vertex_hash = assemble_graph_hash_search_with_key(graph->non_unique_kmers, vertex->data, vertex->data_len, key);

    if (__glibc_unlikely(vertex_hash != NULL)) {
        void* equal_func = graph->vertics_hash->equal_func;

        graph->vertics_hash->equal_func = hc_assemble_utils_working_graph_vertex_equal_func;
        res = assemble_graph_hash_insert_with_key(graph->vertics_hash, vertex->data, vertex->data_len, vertex, key);
        graph->vertics_hash->equal_func = equal_func;
        if (res) {
            p_assemble_vertex search_vertex = res;
            search_vertex->is_dup = 1;
            graph->vertices_sum += 1;
        }
        return res;
    }

    res = assemble_graph_hash_insert(graph->vertics_hash, vertex->data, vertex->data_len, vertex);
    if (res) {
        graph->vertices_sum += 1;
    }
    return res;
}

/**
 * @brief seq graph添加已存在的vertex
 *
 * @param graph     graph
 * @param vertex    vertex
 *
 * @return void*    vertex
 */
static inline void* assemble_graph_add_seq_vertex(p_assemble_graph graph, p_assemble_vertex vertex)
{
    void* res = assemble_graph_hash_insert(graph->vertics_hash, &vertex, sizeof(void*), vertex);
    if (res) {
        graph->vertices_sum += 1;
    }
    return res;
}

/**
 * @brief graph添加已存在的vertex
 *
 * @param graph     graph
 * @param vertex    vertex
 */
p_assemble_vertex assemble_graph_find_vertex(p_assemble_graph graph, void* data, int data_len)
{
    return assemble_graph_hash_search(graph->vertics_hash, data, data_len);
}

/**
 * @brief graph去除已添加的vertex
 *
 * @param graph     graph
 * @param vertex    vertex
 */
void assemble_graph_remove_vertex(p_assemble_graph graph, p_assemble_vertex vertex)
{
    p_assemble_edge el, tmp1, tmp2;

    if (vertex->out_edges) {
        CDL_FOREACH_SAFE2(vertex->out_edges, el, tmp1, tmp2, link.from.prev, link.from.next) { assemble_graph_remove_edge(graph, el); }
    }
    if (vertex->from_edges) {
        CDL_FOREACH_SAFE2(vertex->from_edges, el, tmp1, tmp2, link.to.prev, link.to.next) { assemble_graph_remove_edge(graph, el); }
    }

    if (graph->graph_type == ASSEMBLE_GRAPH_TYPE_READ_THREADING_GRAPH) {
        // read thread graph, seq graph full vertex
        assemble_graph_hash_remove_item(graph->vertics_hash, vertex->data, vertex->data_len, vertex);
    }
    else {
        // seq graph
        assemble_graph_hash_remove_item(graph->vertics_hash, &vertex, sizeof(void*), vertex);
    }
    graph->vertices_sum -= 1;
}

/**
 * @brief 删除edge
 *
 * @param graph graph
 * @param edge edge
 */
void assemble_graph_remove_edge(p_assemble_graph graph, p_assemble_edge edge)
{
    p_assemble_edge el, tmp1, tmp2;
    p_assemble_vertex vertex_el;

    if (!edge->link.all_st.prev || !edge->link.all_st.next) {
        return;
    }

    if (edge->link.from.vertex) {
        vertex_el = edge->link.from.vertex;
        CDL_FOREACH_SAFE2(vertex_el->out_edges, el, tmp1, tmp2, link.from.prev, link.from.next)
        {
            if (el != edge) {
                continue;
            }
            CDL_DELETE2(vertex_el->out_edges, edge, link.from.prev, link.from.next);
            edge->link.from.vertex->out_degree -= 1;
        }
        edge->link.from.vertex = NULL;
    }

    if (edge->link.to.vertex) {
        vertex_el = edge->link.to.vertex;
        CDL_FOREACH_SAFE2(vertex_el->from_edges, el, tmp1, tmp2, link.to.prev, link.to.next)
        {
            if (el != edge) {
                continue;
            }
            CDL_DELETE2(vertex_el->from_edges, edge, link.to.prev, link.to.next);
            edge->link.to.vertex->in_degree -= 1;
        }
        edge->link.to.vertex = NULL;
    }

    CDL_DELETE2(graph->edge_list, edge, link.all_st.prev, link.all_st.next);
    edge->link.all_st.prev = NULL;
    edge->link.all_st.next = NULL;

    graph->edge_sum -= 1;
}

/**
 * @brief 两个Vertex间的直接edge增加weight
 *
 * @param from      vertex
 * @param to        vertex
 * @param weight    weight
 *
 * @return int      是否成功
 */
int assemble_graph_vertex_add_edge_weight_between_vertex(p_assemble_vertex from, p_assemble_vertex to, double weight)
{
    p_assemble_edge find, to_find;
    int res = 0;

    CDL_FOREACH2(from->out_edges, find, link.from.next)
    {
        CDL_FOREACH2(to->from_edges, to_find, link.to.next)
        {
            if (find == to_find) {
                to_find->weight += weight;
                res = 1;
            }
        }
    }
    return res;
}

/**
 * @brief 查找两个Vertex间的直接edge
 *
 * @param from      vertex
 * @param to        vertex
 *
 * @return p_assemble_edge edge
 */
p_assemble_edge assemble_graph_vertex_get_edge_between_vertex(p_assemble_vertex from, p_assemble_vertex to)
{
    p_assemble_edge find, to_find;

    CDL_FOREACH2(from->out_edges, find, link.from.next)
    {
        CDL_FOREACH2(to->from_edges, to_find, link.to.next)
        {
            if (find == to_find) {
                return find;
            }
        }
    }
    return NULL;
}

/**
 * @brief 在两个vertex间添加edge
 *
 * @param graph     graph
 * @param from      from
 * @param to        to
 * @param weight    weight
 * @return p_assemble_edge  edge
 */
p_assemble_edge assemble_graph_vertex_add_edge_to_vertex(p_assemble_graph graph, p_assemble_vertex from, p_assemble_vertex to,
                                                         double weight)
{
    p_assemble_edge edge;

    edge = assemble_graph_edge_create(graph, from, to, weight);
    if (!edge) {
        return NULL;
    }

    if (from) {
        CDL_APPEND2(from->out_edges, edge, link.from.prev, link.from.next);
        from->out_degree++;
    }
    if (to) {
        CDL_APPEND2(to->from_edges, edge, link.to.prev, link.to.next);
        to->in_degree++;
    }
    CDL_APPEND2(graph->edge_list, edge, link.all_st.prev, link.all_st.next);
    graph->edge_sum += 1;

    if ((from && from->storage) && (to && to->storage)) {
        edge->is_ref = ((p_hc_assemble_graph_kmer)(from->storage))->is_ref && ((p_hc_assemble_graph_kmer)(to->storage))->is_ref;
    }

    return edge;
}

void assemble_graph_reset(p_assemble_graph graph)
{
    assemble_graph_hash_reset(graph->vertics_hash);
    assemble_graph_hash_reset(graph->non_unique_kmers);
    ring_mempool_auto_enlarge_reset(graph->vertex_pool);
    ring_mempool_auto_enlarge_reset(graph->edge_pool);
    graph->edge_list = NULL;
    hc_shared_list_reset(graph->alt_path);
    hc_shared_list_reset(graph->ref_path);
    graph->vertices_sum = 0;
    graph->edge_sum = 0;
}

void assemble_graph_free(p_assemble_graph graph)
{
    assemble_graph_hash_destroy(graph->vertics_hash);
    assemble_graph_hash_destroy(graph->non_unique_kmers);
    ring_mempool_auto_enlarge_term(graph->vertex_pool);
    ring_mempool_auto_enlarge_term(graph->edge_pool);
    hc_shared_list_term(graph->alt_path);
    hc_shared_list_term(graph->ref_path);
    hc_assemble_base_graph_iter_finit(&graph->graph_iter);
    ring_mempool_auto_enlarge_term(graph->simplify_graph.edge_pool);
    hc_shared_list_term(graph->simplify_graph.vertex);
    hc_shared_list_term(graph->simplify_graph.edge);
    free(graph);
}
