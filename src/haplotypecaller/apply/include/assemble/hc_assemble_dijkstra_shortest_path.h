#ifndef __SHARED_HC_ASSEMBLE_DIJKSTRA_SHORTEST_PATH_H__
#define __SHARED_HC_ASSEMBLE_DIJKSTRA_SHORTEST_PATH_H__

#include "graph.h"
#include "list.h"
#include "rbtree_shared.h"
#define HC_ASSEMBLE_SW_PAD_STRING   "NNNNNNNNNN"
#define HC_ASSEMBLE_SW_ST_CLEAN_LEN offsetof(hc_assemble_dijkstra_out_path_storage, readIndelRange)

typedef struct hc_assemble_dijkstra_index_range_t
{
    int from;
    int to;
} hc_assemble_dijkstra_index_range, *p_hc_assemble_dijkstra_index_range;

typedef struct hc_assemble_dijkstra_out_path_storage_t
{
    uint32_t seq_len;
    uint32_t cigar_len;
    uint32_t cigar_back_len;

    // cigar builder
    uint32_t section;
    uint32_t lastOperator;
    uint32_t leadingDeletionBasesRemoved;
    uint32_t trailingDeletionBasesRemoved;
    uint32_t trailingDeletionBasesRemovedInMake;

    int start_shift;
    int end_shift;

    uint32_t alignmentStartHapwrtRef;
    hc_assemble_dijkstra_index_range genomeLocation;

    hc_assemble_dijkstra_index_range refIndelRange;
    hc_assemble_dijkstra_index_range readIndelRange;

    uint8_t seq[ASSEMBLE_GRAPH_MAX_SEQ_LEN];
    uint32_t cigar[ASSEMBLE_GRAPH_MAX_SEQ_LEN / 2];
    uint32_t cigar_back[ASSEMBLE_GRAPH_MAX_SEQ_LEN / 2];
} hc_assemble_dijkstra_out_path_storage, *p_hc_assemble_dijkstra_out_path_storage;
/*TODO:A*/
typedef struct assemble_dijkstra_one_path_node_t
{
    uint32_t kmer_size;
    uint32_t is_reference;
    uint32_t done;
    double score;

    p_hc_shared_lib_list_item edges_in_order;
    p_assemble_vertex last_vertex;

    p_hc_assemble_dijkstra_out_path_storage st;

    struct rb_node node;
} assemble_dijkstra_one_path_node, *p_assemble_dijkstra_one_path_node;

typedef struct assemble_dijkstra_one_path_node_t assemble_haplotype_t, assemble_haplotype, *p_assemble_haplotype;

typedef struct assemble_dijkstra_path_t
{
    p_assemble_vertex start_vertex;
    // cacl
    struct rb_root queue;

    // ret
    p_hc_shared_lib_list_head result;
    uint32_t result_size;
    uint32_t all_result_size;
    uint32_t ref_len;

    p_assemble_haplotype ref_haplotype;

    p_mem_pool_fast_main_storage path_pool;   // assemble_dijkstra_one_path_node
    p_auto_enlarge_mempool path_out_st_pool;  // hc_assemble_dijkstra_out_path_storage
    p_mem_pool_fast_main_storage path_node_buffer;

    uint8_t* cached_seq;
    uint8_t* cached_ref;
    uint8_t* ref;
} assemble_dijkstra_path, *p_assemble_dijkstra_path;

#endif