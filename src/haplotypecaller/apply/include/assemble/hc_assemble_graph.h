#ifndef __SHARED_HC_ASSEMBLE_GRAPH_H__
#define __SHARED_HC_ASSEMBLE_GRAPH_H__

#include "hc_assemble_dijkstra_shortest_path.h"
#include "hc_assemble_simple_thread_mem.h"
#include "list.h"
#include "pair_list.h"
#include "uthash.h"
#define HC_ASSEMBLE_GRAPH_POOL_PER_REGION         (65536)
#define HC_ASSEMBLE_DOWN_SAMPLE_MAX_SLOT_PER_BASE (200)

typedef struct hc_assemble_graph_kmer_t hc_assemble_graph_kmer, *p_hc_assemble_graph_kmer;
typedef struct hc_apply_one_read_t hc_apply_one_read, *p_hc_apply_one_read;

enum HC_KMER_HASH_STATUS {
    HC_KMER_HASH_STATUS_DUP = 0,
    HC_KMER_HASH_STATUS_MATCHED = 1,
};

/**
 * @brief Read->Kmer缓存
 */
typedef struct hc_assemble_graph_read_t
{
    const uint8_t* seq;
    uint8_t* qual;

    const uint8_t* dup_seq;
    uint32_t dup_stop;

    uint32_t len : 29;
    uint32_t is_ref : 1;
    uint32_t is_seq : 1;
    uint32_t is_dup : 1;
    uint32_t start;
    uint32_t end;
    char* name;
    struct hc_assemble_graph_read_t* next;
    struct hc_assemble_graph_read_t* prev;
} hc_assemble_graph_read, *p_hc_assemble_graph_read;

/**
 * @brief Kmer缓存
 */
typedef struct hc_assemble_graph_kmer_t
{
    p_hc_assemble_graph_read item;

    uint8_t* seq;

    uint32_t len : 29;
    uint32_t is_ref : 1;
    uint32_t is_seq : 1;
    uint32_t is_dup : 1;

    // read link
    struct hc_assemble_graph_kmer_t* prev;
    struct hc_assemble_graph_kmer_t* next;
} hc_assemble_graph_kmer, *p_hc_assemble_graph_kmer;

/**
 * @brief Graph计算主储存
 */
typedef struct hc_assemble_read_graph_t
{
    uint32_t graph_reads_count;
    p_hc_assemble_graph_read graph_reads;
    p_assemble_graph_hash_table read_hash;

    uint32_t kmer_non_unique_sum;

    hc_assemble_graph_read ref;
    p_assemble_vertex ref_source;
    p_assemble_vertex ref_end;

    int kmer;

    p_auto_enlarge_mempool item_buffer;
    p_auto_enlarge_mempool kmer_buffer;
    p_hc_shared_lib_list_head list_buffer_1;
    p_hc_shared_lib_list_head list_buffer_2;
    p_assemble_graph_hash_table hash_buffer_1;
    p_assemble_graph_hash_table hash_buffer_2;
    struct
    {
        p_assemble_vertex prefix_vertex;
        p_assemble_vertex suffix_vertex;
        p_assemble_vertex prefix_vertex_main_graph;
        p_assemble_vertex suffix_vertex_main_graph;
        p_hc_shared_lib_list_head edges_to_remove;
        p_hc_shared_lib_list_head all_vertex;

        p_hc_shared_lib_pair_list_head middles;

        p_assemble_graph seq_graph;
    } vertex_spliter;

    assemble_dijkstra_path dijkstra_path_finder;

    p_hc_assemble_thread_mem cache_pool;
} hc_assemble_read_graph, *p_hc_assemble_read_graph;

#endif  // !__SHARED_HC_ASSEMBLE_GRAPH_H__