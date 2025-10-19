/**
 * @file hc_apply_main.c
 * @author caoce (caoce@genomics.cn)
 * @brief HC 计算Apply计算
 * @version 0.1
 * @date 2021-05-18
 *
 * @copyright Copyright (c) 2021 ROVACA SDK
 *
 */
#include "hc_assemble_main.h"

#include "debug.h"
#include "hc_marco.h"
#include "mem_pool_auto_enlarge.h"

#define HC_APPLY_MAX_FINDER_NODE   (65536)
#define HC_APPLY_MAX_HAPLOTYPE_SUM (258)

pApplyArgument apply_arguments;

/**
 * @brief HC运算Region, 主函数
 *
 * @param region    region
 */
void hc_apply_main(p_hc_apply apply)
{
    hc_assemble_dijkstra_reset_result(&apply->read_graph.dijkstra_path_finder);
    hc_assemble_reads(apply);
}

#define HC_READS_DATA_BUFFER_SIZE 1024 * 1024 * 128
#define AUX_SIZE                  1024 * 1024 * 4
#define HC_BUFF_PADDING           128

/**
 * @brief HC运算Region, 初始化
 *
 * @return p_hc_apply   apply
 */
p_hc_apply hc_apply_init()
{
    p_hc_apply ret = malloc(sizeof(hc_apply));
    if (!ret) {
        return NULL;
    }
    memset(ret, 0, sizeof(hc_apply));

    if (__glibc_unlikely(!hc_assemble_utils_init_basic_read_graph(&ret->read_graph))) {
        free(ret);
        return NULL;
    }
    ret->reads_count = 0;
    ret->read_graph.graph_reads_count = 0;
    ret->read_thread_graph = assemble_graph_create();
    ret->read_graph.dijkstra_path_finder.path_node_buffer = mem_pool_fast_init(sizeof(hc_shared_lib_list_item), HC_APPLY_MAX_FINDER_NODE);
    ret->read_graph.dijkstra_path_finder.path_pool = mem_pool_fast_init(sizeof(assemble_dijkstra_one_path_node), HC_APPLY_MAX_READ_BUFFER);
    ret->read_graph.dijkstra_path_finder.path_out_st_pool =
        ring_mempool_auto_enlarge_init(sizeof(hc_assemble_dijkstra_out_path_storage), HC_APPLY_MAX_HAPLOTYPE_SUM);
    ret->sw_avx = sw_avx_init(1);
    return ret;
}

/**
 * @brief set target region and reference
 * note: ref now refers to the reference on the target region.
 */
int hc_apply_set_target(p_hc_apply apply, p_hc_region_active_storage region, uint8_t* ref, uint32_t char_ref_len)
{
    if (!(apply && region)) {
        return -1;
    }

    apply->ref_start = assemble_max((int32_t)region->paddedSpan.start - HC_APPLY_REFERENCE_PADDING, (int32_t)1);
    apply->ref_end = assemble_min(char_ref_len, region->paddedSpan.end + HC_APPLY_REFERENCE_PADDING);
    apply->ref_start -= 1;
    apply->ref_end -= 1;

    apply->region = region;
    apply->ref = ref + apply->ref_start;
    apply->ref_len = apply->ref_end - apply->ref_start + 1;

    apply->span_start = apply->region->paddedSpan.start;
    apply->span_end = apply->region->paddedSpan.end;
    apply->ref_haplotype = ref + region->paddedSpan.start - 1;
    apply->ref_hap_len = region->paddedSpan.end - region->paddedSpan.start + 1;
    apply->ref_chr_len = char_ref_len;
    apply->reads_count = 0;
    return 0;
}

int hc_apply_set_region_reads(p_hc_apply apply, p_hc_apply_one_read run_read)
{
    if (!run_read) {
        return -1;
    }
    CDL_APPEND(apply->reads, run_read);
    apply->reads_count++;
    return 0;
}

/**
 * @brief HC运算Region, 反初始化
 *
 * @param p_hc_apply    apply
 */
void hc_apply_finit(p_hc_apply apply)
{
    hc_assemble_utils_finit_basic_read_graph(&apply->read_graph);
    mem_pool_fast_term(apply->read_graph.dijkstra_path_finder.path_node_buffer);
    mem_pool_fast_term(apply->read_graph.dijkstra_path_finder.path_pool);
    ring_mempool_auto_enlarge_term(apply->read_graph.dijkstra_path_finder.path_out_st_pool);
    assemble_graph_free(apply->read_thread_graph);
    sw_avx_finit(apply->sw_avx);
    free(apply);
}

void hc_apply_reset(p_hc_apply apply)
{
    p_hc_apply_one_read read, tmp1, tmp2;
    CDL_FOREACH_SAFE(apply->reads, read, tmp1, tmp2)
    {
        if (read->mempolicy) {
            continue;
        }
        free(read->cigar_str);
        free(read->cigar_cache);
        free(read->cigar_cache_back);
        free(read->seq_cache);
        free(read);
    }

    apply->reads = NULL;
    apply->reads_count = 0;
}
