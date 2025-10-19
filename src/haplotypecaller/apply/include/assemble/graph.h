//
//  graph.h
//  adjacencylist
//
//  Created by Sam Goldman on 6/21/11.
//  Copyright 2011 Sam Goldman. All rights reserved.
//
#ifndef __ASSEMBLE_LIB_GRAPH_H__
#define __ASSEMBLE_LIB_GRAPH_H__

#include "hash.h"
#include "list.h"
#include "mem_pool_auto_enlarge.h"
#include "pair_list.h"

#define ASSEMBLE_GRAPH_VERTEX_INIT_NUM (8192)
#define ASSEMBLE_GRAPH_EDGE_INIT_NUM   (8192)
#define ASSEMBLE_GRAPH_MAX_SEQ_LEN     (2048)

typedef struct assemble_edge_t assemble_edge, *p_assemble_edge;
typedef struct assemble_vertex_t assemble_vertex, *p_assemble_vertex;  // MultiDeBruijnVertex
typedef struct hc_apply_t hc_apply, *p_hc_apply;

enum ASSEMBLE_GRAPH_TYPE_T { ASSEMBLE_GRAPH_TYPE_READ_THREADING_GRAPH, ASSEMBLE_GRAPH_TYPE_SEQ_GRAPH, ASSEMBLE_GRAPH_TYPE_MAX };

enum ASSEMBLE_GRAPH_DIRECTION_T { ASSEMBLE_GRAPH_DIRECTION_FROM, ASSEMBLE_GRAPH_DIRECTION_TO };

enum ASSEMBLE_GRAPH_VERTEX_COLOR_T {
    ASSEMBLE_GRAPH_VERTEX_COLOR_NONE = 0,
    ASSEMBLE_GRAPH_VERTEX_COLOR_RED = 1,
    ASSEMBLE_GRAPH_VERTEX_COLOR_YELLOW = 2,
    ASSEMBLE_GRAPH_VERTEX_COLOR_GREEN = 4,
    ASSEMBLE_GRAPH_VERTEX_COLOR_ALL = 8
};

/**
 * @brief 图深度优先遍历，无迭代
 */
typedef struct assemble_graph_iter_t
{
    struct
    {
        p_assemble_vertex from;  // from
        p_assemble_vertex to;    // to

        p_hc_apply apply;  // apply
    } input_layer;

    p_mem_pool_fast_main_storage run_mem;         // 正在运行的路径内存
    p_hc_shared_lib_list_item run_path;           // 正在运行的路径
    p_mem_pool_fast_main_storage node_mem;        // 尚未探索的转接节点内存
    p_hc_shared_lib_pair_list_item node_stack;    // 尚未探索的转接节点
    p_hc_shared_lib_list_head safe_nodes;         // 安全的转接节点
    p_assemble_graph_hash_table hash_safe_nodes;  // hash safe

    int (*cycle_detected)(p_assemble_vertex, struct assemble_graph_iter_t*);  // cycle
    int (*one_path_done)(p_assemble_vertex, struct assemble_graph_iter_t*);   // 完成一个Path
    int (*all_path_done)(p_assemble_vertex, struct assemble_graph_iter_t*);   // 所有Path都完成
} assemble_graph_iter, *p_assemble_graph_iter;

typedef struct assemble_edge_t
{
    double weight;
    uint32_t is_ref;

    struct
    {
        // for from vertex
        struct
        {
            p_assemble_vertex vertex;
            p_assemble_edge prev;
            p_assemble_edge next;
        } from;

        // for out vertex
        struct
        {
            p_assemble_vertex vertex;
            p_assemble_edge prev;
            p_assemble_edge next;
        } to;

        // for st
        struct
        {
            p_assemble_edge prev;
            p_assemble_edge next;
        } all_st;

    } link;
} assemble_edge, *p_assemble_edge;

// MultiDeBruijnVertex
typedef struct assemble_vertex_t
{
    uint32_t vertex_len;  // vertex(seq) tag len
    uint32_t data_len;    // data(seq) len
    uint32_t call;        // graph call
    uint32_t is_dup;      // dup vertex
    void* data;           // vertex data
    void* storage;        // vertex st

    p_assemble_edge from_edges;
    p_assemble_edge out_edges;

    int in_degree;
    int out_degree;

    uint32_t visited;  // iter visit
    uint32_t color;    // color
    uint32_t key;
} assemble_vertex, *p_assemble_vertex;

typedef struct assemble_simplify_graph_storage_t
{
    p_hc_shared_lib_list_head vertex;
    p_hc_shared_lib_list_head edge;

    uint32_t vertex_sum;
    uint32_t edge_sum;

    p_auto_enlarge_mempool edge_pool;
} assemble_simplify_graph_storage, *p_assemble_simplify_graph_storage;

typedef struct assemble_graph_t
{
    p_assemble_graph_hash_table vertics_hash;
    p_assemble_edge edge_list;

    uint32_t vertices_sum;
    uint32_t edge_sum;

    p_assemble_graph_hash_table non_unique_kmers;
    p_auto_enlarge_mempool vertex_pool;
    p_auto_enlarge_mempool edge_pool;

    assemble_graph_iter graph_iter;

    uint32_t alt_size;
    uint32_t ref_size;
    p_hc_shared_lib_list_head alt_path;
    p_hc_shared_lib_list_head ref_path;

    uint32_t graph_type;
    uint32_t seq_len;
    uint32_t ref_len;
    uint32_t cigar_len;
    uint32_t origin_ref_start;
    uint32_t origin_ref_len;
    // TODO:pooling
    uint8_t seq[ASSEMBLE_GRAPH_MAX_SEQ_LEN];
    uint8_t ref[ASSEMBLE_GRAPH_MAX_SEQ_LEN];
    uint8_t origin_ref[ASSEMBLE_GRAPH_MAX_SEQ_LEN];
    uint32_t cigar[ASSEMBLE_GRAPH_MAX_SEQ_LEN / 2];

    assemble_simplify_graph_storage simplify_graph;
} assemble_graph, *p_assemble_graph;

p_assemble_graph assemble_graph_create(void);
p_assemble_vertex assemble_graph_find_vertex(p_assemble_graph graph, void* data, int data_len);
p_assemble_vertex assemble_graph_vertex_create(p_assemble_graph graph, void* data, uint32_t data_len, void* storage);
// int assemble_graph_vertex_seq_equal_func(void* old_data, void* new_data);
int assemble_graph_vertex_equal_func(void* old_data, void* new_data);

void* assemble_graph_add_vertex(p_assemble_graph graph, p_assemble_vertex vertex);
void assemble_graph_remove_vertex(p_assemble_graph graph, p_assemble_vertex vertex);
void assemble_graph_remove_edge(p_assemble_graph graph, p_assemble_edge edge);

p_assemble_edge assemble_graph_vertex_add_edge_to_vertex(p_assemble_graph graph, p_assemble_vertex from, p_assemble_vertex to,
                                                         double weight);
int assemble_graph_vertex_add_edge_weight_between_vertex(p_assemble_vertex from, p_assemble_vertex to, double weight);
p_assemble_edge assemble_graph_vertex_get_edge_between_vertex(p_assemble_vertex from, p_assemble_vertex to);

void assemble_graph_to_normal_read_thread_graph(p_assemble_graph graph);
void assemble_graph_to_working_read_thread_graph(p_assemble_graph graph);
void assemble_graph_to_seq_graph(p_assemble_graph graph);

void assemble_graph_reset(p_assemble_graph graph);
void assemble_graph_free(p_assemble_graph graph);

#endif  // !__ASSEMBLE_LIB_GRAPH_H__