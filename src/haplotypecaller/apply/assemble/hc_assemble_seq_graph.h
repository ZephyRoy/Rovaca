#ifndef __PRIVATE_HC_ASSEMBLE_SEQ_GRAPH_H__
#define __PRIVATE_HC_ASSEMBLE_SEQ_GRAPH_H__

#include "graph.h"
#include "hc_apply.h"
#include "list.h"
/**
 * How many cycles of the graph simplifications algorithms will we run before
 * thinking something has gone wrong and throw an exception?
 */
#define HC_ASSEMBLE_SEQ_GREAPH_MAX_REASONABLE_SIMPLIFICATION_CYCLES (100)
#define HC_ASSEMBLE_SEQ_GREAPH_MAX_CHANGE_LOOPS                     (5)

#define HC_ASSEMBLE_SEQ_GRAPH_SIMLIFY_DECLARE(name)                 static int hc_assemble_seq_graph_##name(p_hc_apply apply);

/**
 * Merge until the graph has no vertices that are candidates for merging
 */
#define HC_ASSEMBLE_SEQ_GRAPH_SIMLIFY(name, fun)                                                                                 \
    static int hc_assemble_seq_graph_##name(p_hc_apply apply)                                                                    \
    {                                                                                                                            \
        p_assemble_graph_hash_node vertex_node;                                                                                  \
        p_hc_shared_lib_list_item node;                                                                                          \
        int did_at_least_one_transform = 0;                                                                                      \
        int found_nodes_to_merge = 1;                                                                                            \
                                                                                                                                 \
        while (found_nodes_to_merge) {                                                                                           \
            p_hc_shared_lib_list_head all_head = hc_assemble_utils_check_list_work(apply->read_graph.vertex_spliter.all_vertex); \
                                                                                                                                 \
            found_nodes_to_merge = 0;                                                                                            \
            CDL_FOREACH(apply->read_thread_graph->vertics_hash->node_list, vertex_node)                                          \
            {                                                                                                                    \
                hc_shared_list_insert(all_head, vertex_node->value);                                                             \
            }                                                                                                                    \
                                                                                                                                 \
            CDL_FOREACH(all_head->nodes, node)                                                                                   \
            {                                                                                                                    \
                found_nodes_to_merge = fun(node->data, apply);                                                                   \
                                                                                                                                 \
                if (found_nodes_to_merge) {                                                                                      \
                    did_at_least_one_transform = 1;                                                                              \
                    break;                                                                                                       \
                }                                                                                                                \
            }                                                                                                                    \
        }                                                                                                                        \
        return did_at_least_one_transform;                                                                                       \
    }

static inline p_hc_shared_lib_list_head hc_assemble_seq_graph_check_list_work(p_hc_shared_lib_list_head list);
static int hc_assemble_seq_graph_is_linear_chain_start(p_assemble_vertex vertex);
static inline int hc_assemble_seq_graph_is_reference_node(p_assemble_vertex vertex);
static p_hc_shared_lib_list_head hc_assemble_seq_graph_trace_linear_chain(p_assemble_vertex zipStart, p_hc_apply apply);
static p_assemble_vertex hc_assemble_seq_graph_merge_linear_chain_vertices(p_hc_shared_lib_list_head vertices, p_hc_apply apply);
static int hc_assemble_seq_graph_merge_linear_chain(p_hc_shared_lib_list_head linearChain, p_hc_apply apply);
static int hc_assemble_seq_graph_travel_one_path_done(p_assemble_vertex last_vertex, p_assemble_graph_iter iter_st);
static void hc_assemble_seq_graph_travel_graph(p_assemble_vertex ref_vertex, p_hc_apply apply);
static int hc_assemble_seq_graph_simplify_graph_once(p_hc_apply apply);
static void hc_assemble_seq_graph_clone_graph(p_assemble_graph from);
static inline int hc_assemble_seq_graph_compare_vertex(p_hc_shared_lib_list_item a, p_assemble_vertex b);
static inline int hc_assemble_seq_graph_compare_edge(p_hc_shared_lib_list_item from_a, p_assemble_edge b);
static int hc_assemble_seq_graph_equal_graphs(p_assemble_graph graph);
static uint32_t hc_assemble_seq_graph_get_reference_bytes(p_assemble_vertex fromVertex, p_assemble_vertex toVertex, p_hc_apply apply);
static p_assemble_vertex hc_assemble_seq_graph_get_next_reference_vertex(p_assemble_vertex vertex);

#endif