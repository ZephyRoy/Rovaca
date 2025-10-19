/**
 * @file hc_assemble_reads.c
 * @author Cao Ce (caoce@genomics.cn)
 * @brief 在给定区域的读取数据上运行汇编程序的高级功能，返回数据结构，并提供后续HC步骤所需的结果信息
 * @version 0.1
 * @date 2021-05-19
 *
 * @copyright Copyright (c) 2021 ROVACA SDK
 *
 */

#include "assemble_interface.h"
#include "hc_apply.h"
#include "hc_assemble.h"
#include "hc_assemble_cigar_utils.h"
#include "hc_assemble_debug_print.h"
#include "hc_func.h"
#include "hc_marco.h"
#include "rbtree.h"
#include "rbtree_shared.h"

static void hc_assemble_hash_new_node_memcpy(void* new_node, void* input_pt, int hash_len, p_assemble_graph_hash_table hash_table);
static void hc_assemble_correct_overlapping_base_qualities(p_hc_apply apply);
static inline int hc_assemble_reads_compare_int(int a, int b);
static int hc_assemble_reads_compare_read_coodinate(p_hc_apply_one_read first, p_hc_apply_one_read second);
static int hc_assemble_reads_insert_read_by_rbtree(struct rb_root* root, p_hc_apply_one_read read);
static void hc_assemble_reads_sort_by_coordinates(p_hc_apply apply);

#define ASSEMBLE_CHECK_VALID_READ(read) (read->cigar_len && read->read_len)
/**
 * @brief HC 组装
 *
 * @param apply     apply储存
 */
void hc_assemble_reads(p_hc_apply apply)
{
    hc_assemble_reads_sort_by_coordinates(apply);
    if (__glibc_unlikely(apply_arguments->debug.debug_print)) {
        hc_assemble_utils_print_reads_after_finalize_region(apply);
    }
    hc_assemble_correct_overlapping_base_qualities(apply);

    if (!apply->region->active) {
        return;
    }
    int res = 0;
    for (size_t i = 0; i < apply_arguments->kmers->nkmer; ++i) {
        apply->read_graph.kmer = apply_arguments->kmers->kmer[i];
        res += hc_assemble_graph_build(apply);
    }
    int iter = 1;
    while (__glibc_unlikely(!res)) {
        apply->read_graph.kmer += HC_ASSEMBLE_KMER_ITERATION_INCREASE;
        res = hc_assemble_graph_build(apply);
        if (res || iter == HC_ASSEMBLE_MAX_KMER_ITERATIONS_TO_ATTEMPT) {
            break;
        }
        iter += 1;
    }
    if (__glibc_unlikely(apply_arguments->debug.debug_print)) {
        hc_assemble_utils_print_path(apply);
    }
}

/**
 * @brief HC Finalize
 *
 * @param read          单条read
 * @param apply         apply储存
 */
int hc_assemble_finalize_region(p_hc_apply apply, p_hc_apply_one_read read)
{
    if (hc_apply_utils_has_well_defined_fragment_size(read)) {
        hc_assemble_utils_revert_soft_clip(read);
    }
    else {
        hc_assemble_utils_soft_clip_to_hard(read);
    }
    hc_assemble_utils_hc_read_to_cache(read);
    if (!ASSEMBLE_CHECK_VALID_READ(read)) {
        return 0;
    }

    hc_assemble_utils_clip_low_qual_ends(read, apply_arguments->soft_clip_low_quality_ends ? SOFTCLIP : HARDCLIP);
    if (!ASSEMBLE_CHECK_VALID_READ(read)) {
        return 0;
    }

    if (!hc_assemble_utils_read_is_unmapped(read)) {
        hc_apply_utils_hard_clip_adaptor_sequence(read);
        if (!ASSEMBLE_CHECK_VALID_READ(read)) {
            return 0;
        }
    }
    hc_assemble_utils_hard_clip_to_region(read, apply->span_start, apply->span_end);
    if (!ASSEMBLE_CHECK_VALID_READ(read)) {
        return 0;
    }
    merge_consecutive_identical_cigar(read);
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
static void hc_assemble_hash_new_node_memcpy(void* new_node, void* input_pt, int hash_len, p_assemble_graph_hash_table hash_table)
{
    (void)hash_table;
    p_assemble_graph_hash_node insert_node = new_node;
    p_hc_apply_one_read read = input_pt;

    memset(insert_node, 0, sizeof(assemble_graph_hash_node));

    insert_node->key = read->read_data;
    insert_node->keylen = hash_len;
    insert_node->value = input_pt;
}

/**
 * @brief 修正配对重叠的read内碱基质量
 *
 * @param[in] apply apply
 */
static void hc_assemble_correct_overlapping_base_qualities(p_hc_apply apply)
{
    p_hc_apply_one_read read, tmp1, tmp2;
    hts_pos_t mate_start;
    p_hc_apply_one_read same_read;
    p_assemble_graph_hash_node node;
    uint32_t key;
    int len;
    void* memcpy_fun = apply->read_graph.read_hash->equal_func;
    apply->read_graph.read_hash->new_node_memcpy = hc_assemble_hash_new_node_memcpy;
    CDL_FOREACH_SAFE(apply->reads, read, tmp1, tmp2)
    {
        mate_start = read->read.core.mpos + 1;
        if (!PAIRED(&read->read) || MATE_UNMAPPED(&read->read) || mate_start == 0 || mate_start > read->pos_end) {
            continue;
        }
        // 在pair reads cache中寻找qname相同的read.
        len = strlen((const char*)(read->read_data));
        key = assemble_graph_hash_get_key(apply->read_graph.read_hash, read->read_data, len);
        node = assemble_graph_hash_search_with_key(apply->read_graph.read_hash, read->read_data, len, key);
        if (!node) {
            assemble_graph_hash_insert_with_key(apply->read_graph.read_hash, read->read_data, len, read, key);
            continue;
        }
        same_read = node->value;
        hc_assemble_utils_adjust_overlapping_paired_qual(same_read, read);
        assemble_graph_hash_remove_item_with_key(apply->read_graph.read_hash, same_read, key);
    }
    apply->read_graph.read_hash->equal_func = memcpy_fun;
    assemble_graph_hash_reset(apply->read_graph.read_hash);
}

static inline int hc_assemble_reads_compare_int(int a, int b)
{
    if (a > b) {
        return 1;
    }
    else if (a < b) {
        return -1;
    }
    else {
        return 0;
    }
}

static int hc_assemble_reads_compare_read_coodinate(p_hc_apply_one_read first, p_hc_apply_one_read second)
{
    int result = hc_assemble_utils_compare_coordinates(first, second);
    if (result != 0) {
        return result;
    }

    // This is done to mimic SAMRecordCoordinateComparator's behavior
    if (REVERSE_STRAND(&first->read) != REVERSE_STRAND(&second->read)) {
        return REVERSE_STRAND(&first->read) ? 1 : -1;
    }

    if ((result = strcmp((char*)first->read_data, (char*)second->read_data)) != 0) {
        return result;
    }
    result = hc_assemble_reads_compare_int(first->read.core.flag, second->read.core.flag);
    if (result != 0) {
        return result;
    }
    result = hc_assemble_reads_compare_int(first->read.core.qual, second->read.core.qual);
    if (result != 0) {
        return result;
    }
    if (PAIRED(&first->read) && PAIRED(&second->read)) {
        result = hc_assemble_reads_compare_int(first->read.core.mtid, second->read.core.mtid);
        if (result != 0) {
            return result;
        }
        result = hc_assemble_reads_compare_int(first->read.core.mpos, second->read.core.mpos);
        if (result != 0) {
            return result;
        }
    }
    result = hc_assemble_reads_compare_int(first->read.core.isize, second->read.core.isize);

    return result;
}

static int hc_assemble_reads_insert_read_by_rbtree(struct rb_root* root, p_hc_apply_one_read read)
{
    struct rb_node** new = &(root->rb_node), *parent = NULL;

    /* Figure out where to put new node */
    while (*new) {
        p_hc_apply_one_read this = container_of(*new, hc_apply_one_read, node);
        int result = hc_assemble_reads_compare_read_coodinate(read, this);

        parent = *new;
        if (result < 0)
            new = &((*new)->rb_left);
        else if (result > 0)
            new = &((*new)->rb_right);
        else
            new = &((*new)->rb_right);
    }

    rb_link_node(&read->node, parent, new);
    rb_insert_color(&read->node, root);

    return 1;
}

static void hc_assemble_reads_sort_by_coordinates(p_hc_apply apply)
{
    p_hc_apply_one_read read_to_sort, tmp1, tmp2;

    apply->read_sorted = RB_ROOT;

    CDL_FOREACH_SAFE(apply->reads, read_to_sort, tmp1, tmp2) { hc_assemble_reads_insert_read_by_rbtree(&apply->read_sorted, read_to_sort); }
}