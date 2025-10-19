#ifndef __HC_ASSEMBLE_BASE_GRAPH_H__
#define __HC_ASSEMBLE_BASE_GRAPH_H__
#include "graph.h"

#define HC_ASSEMBLE_BASE_GRAPH_MIN_DANGLING_BRANCH_LENGTH (4)
#define HC_ASSEMBLE_BASE_GRAPH_MAX_CIGAR_COMPLEXITY       (3)

#define HC_ASSEMBLE_GRAPH_GET_PATH_INDEX_ITEM(head, repeat_time, ret) \
    do {                                                              \
        p_hc_shared_lib_list_item tmp;                                \
        uint32_t iter = 0;                                            \
        ret = NULL;                                                   \
        CDL_FOREACH(head, tmp)                                        \
        {                                                             \
            ret = tmp->data;                                          \
            if (iter == repeat_time) {                                \
                break;                                                \
            }                                                         \
            iter += 1;                                                \
        }                                                             \
        if (iter != repeat_time) {                                    \
            ret = NULL;                                               \
        }                                                             \
    } while (0);

static p_hc_shared_lib_list_head hc_assemble_base_graph_find_path_upwards_to_lowest_common_ancestor(p_assemble_vertex vertex,
                                                                                                    p_hc_apply apply);
static p_assemble_edge hc_assemble_base_graph_get_heaviest_incoming_edge(p_assemble_vertex vertex);
static p_assemble_edge hc_assemble_base_graph_get_heaviest_outgoing_edge(p_assemble_vertex vertex);
p_assemble_vertex hc_assemble_base_graph_get_next_reference_vertex(p_assemble_vertex vertex, int allowNonRefPaths,
                                                                   p_assemble_edge blacklistedEdge);
p_assemble_vertex hc_assemble_base_graph_get_prev_reference_vertex(p_assemble_vertex vertex);
static int hc_assemble_base_graph_get_reference_path(p_assemble_vertex start, int direction, p_assemble_edge blacklistedEdge,
                                                     p_hc_apply apply);
void hc_assemble_base_graph_get_bases_for_path(p_hc_shared_lib_list_head ref, p_hc_shared_lib_list_head alt, int expandSource,
                                               p_hc_apply apply);
static void hc_assemble_base_graph_run_sw(p_hc_apply apply);
static void hc_assemble_base_graph_remove_trailing_deletions(p_hc_apply apply);
static int hc_assemble_base_graph_recover_dangling_tail(p_assemble_vertex vertex, p_hc_apply apply);
static int hc_assemble_base_graph_generate_cigar_against_downwards_reference_path(p_assemble_vertex vertex, p_hc_apply apply);
static int hc_assemble_base_graph_cigar_is_okay_to_merge(p_hc_apply apply, int requireFirstElementM, int requireLastElementM);
static int hc_assemble_base_graph_longest_suffix_match(uint8_t* seq, uint8_t* kmer, int seqStart, int kmer_len);
static int hc_assemble_base_graph_merge_dangling_tail(p_hc_apply apply);
static inline int hc_assemble_base_graph_get_max_mismatches_legacy(int lengthOfDanglingBranch, int kmerSize);
static int hc_assemble_base_graph_best_prefix_match_legacy(p_hc_apply apply, int maxIndex);
static int hc_assemble_base_graph_get_offset_for_ref_end_to_dangling_end(uint32_t* cigar, uint32_t cigar_len);
static int hc_assemble_base_graph_extend_dangling_path_against_reference(p_hc_apply apply, int numNodesToExtend, uint32_t* elements);
static int hc_assemble_base_graph_merge_dangling_head_legacy(p_hc_apply apply);

static int hc_assemble_base_graph_recover_dangling_head(p_assemble_vertex vertex, p_hc_apply apply);
static int hc_assemble_base_graph_generate_cigar_against_upwards_reference_path(p_assemble_vertex vertex, p_hc_apply apply);
static int hc_assemble_base_graph_find_path_downwards_to_highest_common_descendant_of_reference(p_assemble_vertex vertex, p_hc_apply apply);
static inline int hc_assemble_base_graph_is_reference_node(p_assemble_vertex vertex);

#endif