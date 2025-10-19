/**
 * @file hc_assemble_utils.c
 * @author Cao Ce (caoce@genomics.cn)
 * @brief HC Assemble实用工具
 * @version 0.1
 * @date 2021-05-19
 *
 * @copyright Copyright (c) 2021 ROVACA SDK
 *
 */

#include "assemble_interface.h"
#include "debug.h"
#include "hash.h"
#include "hc_apply.h"
#include "hc_assemble.h"
#include "hc_assemble_cigar_utils.h"
#include "hc_func.h"
#include "hc_marco.h"

static inline void hc_assemble_utils_repair_first_soft(p_hc_apply_one_read read, uint32_t cigar, uint32_t* cigar_len);
static int hc_assemble_utils_clip_cigar(p_hc_apply_one_read read, int start, int stop, clipping_type clip);
static int hc_assemble_utils_alignment_start_shift(uint32_t* cigar, int cigar_count, int numClipped);

/**
 * @brief 取消SOFT CLIP
 *
 * @param read read
 */
void hc_assemble_utils_revert_soft_clip(p_hc_apply_one_read read)
{
    uint32_t i, this_cigar_len, this_cigar;
    uint32_t* cigar = read->cigar;
    uint32_t origin_cigar_len = read->cigar_len;

    read->cigar_len = 0;
    read->read_len = 0;
    read->ref_len = 0;

    for (i = 0; i < origin_cigar_len; i++) {
        this_cigar_len = bam_cigar_oplen(cigar[i]);
        this_cigar = bam_cigar_op(cigar[i]);

        // 不是SOFT CLIP
        if (this_cigar != WORKER_SIGAR_SOFT_CLIP) {
            read->cigar_cache[read->cigar_len] = cigar[i];
            read->cigar_len += 1;
            if (consumes_read_bases(this_cigar)) {
                read->read_len += this_cigar_len;
            }
            if (consumes_ref_bases(this_cigar)) {
                read->ref_len += this_cigar_len;
            }
            continue;
        }

        // 第一个CIGAR是SOFT，特殊处理
        if (i == 0) {
            hc_assemble_utils_repair_first_soft(read, cigar[i], &this_cigar_len);
        }
        else {
            // 其它CIGAR是SOFT，直接转换
            read->cigar_cache[read->cigar_len] = bam_cigar_gen(this_cigar_len, BAM_CMATCH);
            read->cigar_len += 1;
        }
        read->read_len += this_cigar_len;
        read->ref_len += this_cigar_len;
    }
    return;
    // for debug.
    memset(read->cigar_str, 0, strlen(read->cigar_str));
    for (uint32_t i = 0; i < read->cigar_len; ++i) {
        char op = bam_cigar_opchr(read->cigar_cache[i]);
        int len = bam_cigar_oplen(read->cigar_cache[i]);
        if (i == 0) {
            sprintf(read->cigar_str, "%d%c", len, op);
            continue;
        }
        sprintf(read->cigar_str + strlen(read->cigar_str), "%d%c", len, op);
    }
}

/**
 * @brief HC转换第一个SOFT CLIP
 *
 * @param read              read
 * @param cigar             cigar
 * @param cigar_len         cigar 长度
 */
static inline void hc_assemble_utils_repair_first_soft(p_hc_apply_one_read read, uint32_t cigar, uint32_t* out_cigar_len)
{
    (void)cigar;
    int transform_size, ret_size;
    uint32_t cigar_len = *out_cigar_len;

    // 超出ref
    if (__glibc_unlikely(read->read.core.pos < cigar_len)) {
        transform_size = read->read.core.pos;
        ret_size = cigar_len - transform_size;
        // ref内
    }
    else {
        transform_size = cigar_len;
        ret_size = 0;
    }
    if (__glibc_unlikely(ret_size)) {
        read->seq += ret_size;
        read->qual += ret_size;
        read->cigar_cache[read->cigar_len] = bam_cigar_gen(ret_size, BAM_CHARD_CLIP);
        read->cigar_len += 1;
    }
    read->cigar_cache[read->cigar_len] = bam_cigar_gen(transform_size, BAM_CMATCH);
    read->cigar_len += 1;
    read->pos_start -= transform_size;
    *out_cigar_len = transform_size;
}

/**
 * @brief 获取新的单一cigar
 *
 * @param[in]   cigar_op          cigar_op
 * @param[out]  cigar_len         cigar_len
 *
 * @retval uint32_t               cigar
 **/
inline uint32_t hc_assemble_utils_get_new_cigar(uint8_t cigar_op, uint32_t cigar_len)
{
    uint32_t ret;

    ret = cigar_len << BAM_CIGAR_SHIFT;
    ret |= cigar_op;
    return ret;
}

/**
 * @brief CIGAR，SOFT转为HARD
 *
 * @param read read
 */
void hc_assemble_utils_soft_clip_to_hard(p_hc_apply_one_read read)
{
    uint32_t i, this_cigar_len, this_cigar, read_len = 0;
    uint32_t* cigar = read->cigar;
    uint32_t origin_cigar_len = read->cigar_len;

    read->cigar_len = 0;
    read->read_len = 0;
    read->ref_len = 0;

    for (i = 0; i < origin_cigar_len; i++) {
        this_cigar_len = bam_cigar_oplen(cigar[i]);
        this_cigar = bam_cigar_op(cigar[i]);

        // 不是SOFT CLIP
        if (this_cigar != WORKER_SIGAR_SOFT_CLIP) {
            read->cigar_cache[read->cigar_len] = cigar[i];
            read->cigar_len += 1;
            if (consumes_read_bases(this_cigar)) {
                read->read_len += this_cigar_len;
                read_len += this_cigar_len;
            }
            if (consumes_ref_bases(this_cigar)) {
                read->ref_len += this_cigar_len;
            }
            continue;
        }
        read->cigar_cache[read->cigar_len] = bam_cigar_gen(this_cigar_len, BAM_CHARD_CLIP);
        read->cigar_len += 1;

        if (read->cigar_len == 1) {
            read->seq += this_cigar_len;
            read->qual += this_cigar_len;
        }
        else {
            read_len += this_cigar_len;
        }
    }
    return;
    // for debug.
    memset(read->cigar_str, 0, strlen(read->cigar_str));
    for (uint32_t i = 0; i < read->cigar_len; ++i) {
        char op = bam_cigar_opchr(read->cigar_cache[i]);
        int len = bam_cigar_oplen(read->cigar_cache[i]);
        if (i == 0) {
            sprintf(read->cigar_str, "%d%c", len, op);
            continue;
        }
        sprintf(read->cigar_str + strlen(read->cigar_str), "%d%c", len, op);
    }
}

/**
 * @brief 复制bam序列
 *
 * @param[out]  out             dst
 * @param[in]   in              src
 * @param[in]   start           start
 * @param[in]   len             len
 **/
void hc_assemble_utils_copy_seq(uint8_t* __restrict out, uint8_t* __restrict in, uint32_t start, uint32_t len)
{
    uint32_t i = start;
    register uint8_t str_pos;
    const char* bqsr_utils_seq_str = BQSR_FASTA_SEQ;

    for (; i < len + start; i++) {
        str_pos = bam_seqi(in, i);
        if (likely(str_pos < BQSR_FASTA_SEQ_LEN)) {
            out[i - start] = bqsr_utils_seq_str[str_pos];
        }
    }
}

/**
 * @brief HARD CLIP掉过低的头尾
 *
 * @param read read
 */
void hc_assemble_utils_clip_low_qual_ends(p_hc_apply_one_read read, clipping_type clip_type)
{
    uint32_t clip_start, clip_end;
    uint8_t* qual = read->qual;

    // check how far we can clip both sides
    clip_start = 0;
    clip_end = read->read_len - 1;

    while (qual[clip_end] <= HC_ASSEMBLE_MIN_TAIL_QUALITY_TO_USE) {
        if (clip_end == 0) break;
        clip_end--;
    }
    while (clip_start < read->read_len && qual[clip_start] <= HC_ASSEMBLE_MIN_TAIL_QUALITY_TO_USE) {
        clip_start++;
    }
    if (likely(clip_start == 0 && clip_end == read->read_len - 1)) {
        return;
    }
    else if (__glibc_unlikely(clip_end <= clip_start)) {
        read->cigar_len = 0;
        return;
    }

    if (clip_end != read->read_len - 1) {
        hc_assemble_utils_clip_read(read, clip_end + 1, read->read_len - 1, clip_type);
    }

    if (clip_start != 0) {
        hc_assemble_utils_clip_read(read, 0, clip_start - 1, clip_type);
    }
}

/**
 * Given a cigar string, soft clip up to leftClipEnd and soft clip starting at rightClipBegin
 * @param start initial index to clip within read bases, inclusive
 * @param stop final index to clip within read bases exclusive
 * @param clippingOperator      type of clipping -- must be either hard clip or soft clip
 * @Modification instructions:
 *  使用htslib原生方法bam_cigar_gen(), consumes_read_bases()等代替原来函数调用。
 *  去除自定义的未知cigar operation，用具体clip操作代替，避免二次赋值。
 */
static int hc_assemble_utils_clip_cigar(p_hc_apply_one_read read, int start, int stop, clipping_type clip_type)
{
    uint32_t* cigar = read->cigar;
    uint32_t cigar_len = read->cigar_len;
    uint32_t* new_cigar = read->cigar_back;
    uint32_t operator, operator_len;
    uint32_t i = 0, j = 0;
    int clipLeft = start == 0;
    int elementStart = 0;

    for (; i < cigar_len; i++) {
        operator= bam_cigar_op(cigar[i]);
        operator_len = bam_cigar_oplen(cigar[i]);
        bool consume_read = bam_cigar_type(operator) & 0x1;

        if (operator== BAM_CHARD_CLIP) {
            new_cigar[j] = bam_cigar_gen(operator_len, operator);
            j += 1;
            continue;
        }
        int elementEnd = elementStart + (consume_read ? operator_len : 0);
        if (elementEnd <= start || elementStart >= stop) {
            if (consume_read || (elementStart != start && elementStart != stop)) {
                new_cigar[j] = bam_cigar_gen(operator_len, operator);
                j += 1;
            }
        }
        else {  // otherwise, some or all of the element is soft-clipped
            int unclippedLength = clipLeft ? elementEnd - stop : start - elementStart;
            int clippedLength = operator_len - unclippedLength;

            if (unclippedLength <= 0) {
                if (consume_read) {
                    new_cigar[j] = bam_cigar_gen(operator_len, clip_type);
                    j += 1;
                }
            }
            else if (clipLeft) {
                new_cigar[j] = bam_cigar_gen(clippedLength, clip_type);
                new_cigar[++j] = bam_cigar_gen(unclippedLength, operator);
                j += 1;
            }
            else {
                new_cigar[j] = bam_cigar_gen(unclippedLength, operator);
                new_cigar[++j] = bam_cigar_gen(clippedLength, clip_type);
                j += 1;
            }
        }
        elementStart = elementEnd;
    }
    return (int)j;
}

/**
 * How many bases to the right does a read's alignment start shift given its cigar and the number of left soft clips
 */
static int hc_assemble_utils_alignment_start_shift(uint32_t* cigar, int cigar_count, int numClipped)
{
    int refBasesClipped = 0;
    register int i = 0;
    uint32_t operator, operator_len;

    int elementStart = 0;  // this and elementEnd are indices in the read's bases

    for (; i < cigar_count; i++) {
        // final CigarElement element : cigar.getCigarElements();
        // final CigarOperator operator = element.getOperator();
        operator= bam_cigar_op(cigar[i]);
        operator_len = bam_cigar_oplen(cigar[i]);

        if (operator== WORKER_SIGAR_HARD_CLIP) {
            continue;
        }
        int elementEnd = elementStart + (consumes_read_bases(operator) ? operator_len : 0);

        if (elementEnd <= numClipped) {  // totally within clipped span -- this includes deletions immediately following clipping
            refBasesClipped += consumes_ref_bases(operator) ? operator_len : 0;
        }
        else if (elementStart < numClipped) {  // clip in middle of element, which means the element necessarily consumes read bases
            int clippedLength = numClipped - elementStart;
            refBasesClipped += consumes_ref_bases(operator) ? clippedLength : 0;
            break;
        }
        elementStart = elementEnd;
    }
    return refBasesClipped;
}

/**
 * @brief Hard Clip掉区域
 *
 * @param read          read
 * @param clip_start    头位置
 * @param clip_end      尾位置
 * @Refactoring instructions:
 *  clip cigar 函数接口修改，支持根据clip type修改cigar.
 *  seq & qual内存reset方法修改，根据偏移一次性修改，无需循环.
 */
void hc_assemble_utils_clip_read(p_hc_apply_one_read read, uint32_t clip_start, uint32_t clip_end, clipping_type clip_type)
{
    uint32_t i, this_cigar_len, this_cigar;
    uint8_t* qual = read->qual;
    uint8_t* seq = read->seq;
    uint32_t new_len = read->read_len - (clip_end - clip_start + 1);

    // clip cigar.
    uint32_t new_cigar_len = hc_assemble_utils_clip_cigar(read, clip_start, clip_end + 1, clip_type);
    int copy_start = (clip_start == 0) ? clip_end + 1 : 0;

    // set position.这里计算start shift需要未经过clip的cigar.
    if (clip_start == 0 && !UNMAPPED(&read->read)) {
        read->pos_start += hc_assemble_utils_alignment_start_shift(read->cigar, read->cigar_len, clip_end + 1);
    }
    read->cigar = read->cigar_back;
    read->cigar_back = read->cigar_back == read->cigar_cache ? read->cigar_cache_back : read->cigar_cache;
    read->cigar_len = new_cigar_len;
    uint32_t new_ref_len = 0;
    read->read_len = 0;

    for (i = 0; i < read->cigar_len; i++) {
        this_cigar_len = bam_cigar_oplen(read->cigar[i]);
        this_cigar = bam_cigar_op(read->cigar[i]);
        read->read_len += bam_cigar_type(this_cigar) & 0x1 ? this_cigar_len : 0;
        new_ref_len += bam_cigar_type(this_cigar) & 0x2 ? this_cigar_len : 0;
    }
    read->ref_len = new_ref_len;
    read->pos_end = new_ref_len + read->pos_start - 1;  // read->pos_start已经是1base了.

    memmove(read->seq, seq + copy_start, new_len);
    memmove(read->qual, qual + copy_start, new_len);
    // 此处省略设置indel qual逻辑。
}

/**
 * @brief Read是否为Unmapped
 *
 * @param read read
 * @return int 成功
 */
int hc_assemble_utils_read_is_unmapped(p_hc_apply_one_read read)
{
    return UNMAPPED(&read->read) || read->read.core.tid > HC_R_NAME_MAX || read->read.core.tid < 0;
}

/**
 * @brief 初始化Read结构
 *
 * @param read  read
 */
void hc_assemble_utils_init_hc_read(p_hc_apply_one_read read)
{
    read->cigar = bam_get_cigar(&read->read);
    read->qual = bam_get_qual(&read->read);
    read->cigar_len = read->read.core.n_cigar;
    read->pos_start = read->read.core.pos;
    read->read_len = read->read.core.l_qseq;
    read->insert_size = read->read.core.isize;
    read->ref_len = 0;
}

/**
 * @brief 使用Cache初始化Read结构
 *
 * @param read  read
 */
void hc_assemble_utils_hc_read_to_cache(p_hc_apply_one_read read)
{
    read->cigar = read->cigar_cache;
    read->cigar_back = read->cigar_cache_back;
    read->insert_size = read->read.core.isize;
}

/**
 * Hard clip the read to the variable region (from refStart to refStop)
 *
 * @param read     the read to be clipped
 * @param refStart the beginning of the variant region (inclusive)
 * @param refStop  the end of the variant region (inclusive)
 * @return the read hard clipped to the variant region (Could return an empty, unmapped read)
 */
void hc_assemble_utils_hard_clip_to_region(p_hc_apply_one_read read, int refStart, int refStop)
{
    int start = read->pos_start;
    int stop = hc_assemble_utils_read_get_end(read);

    if (start <= refStop && stop >= refStart) {
        if (start < refStart && stop > refStop) {
            hc_apply_utils_clip_by_reference_coordinates(refStart - 1, refStop + 1, read);
        }

        if (start < refStart) {
            hc_apply_utils_clip_by_reference_coordinates(-1, refStart - 1, read);
        }

        if (stop > refStop) {
            hc_apply_utils_clip_by_reference_coordinates(refStop + 1, -1, read);
        }
    }
    else {
        read->pos_start = 0;
        read->cigar_len = 0;
        read->read_len = 0;
        read->cigar_len = 0;
    }
}

/**
 * @brief 通过read及kmer start生产kmer
 *
 * @param read          read
 * @param kmerStart     kmer start
 * @param apply         apply
 *
 * @return p_hc_assemble_graph_kmer kmer(always)
 */
p_hc_assemble_graph_kmer hc_assemble_utils_fill_kmer_by_read(p_hc_assemble_graph_read read, int kmerStart, p_hc_apply apply)
{
    p_hc_assemble_graph_kmer kmer_item;

    kmer_item = ring_mempool_auto_enlarge_malloc(apply->read_graph.kmer_buffer);
    if (__glibc_unlikely(!kmer_item)) {
        statics_print_error("No Mem");
        exit(-1);
    }
    kmer_item->item = read;
    kmer_item->seq = (uint8_t*)(read->seq + kmerStart);
    kmer_item->len = apply->read_graph.kmer;
    kmer_item->is_ref = read->is_ref;
    kmer_item->is_seq = read->is_seq;
    kmer_item->is_dup = read->is_dup;
    kmer_item->prev = NULL;
    kmer_item->next = NULL;

    return kmer_item;
}

/**
 * @brief 通过Kmer建立Vertex并加入graph
 *
 * @param kmer_item         kmer
 * @param apply             apply
 *
 * @return p_assemble_vertex vertex(always)
 */
p_assemble_vertex hc_assemble_utils_add_graph_vertex_by_kmer(p_hc_assemble_graph_kmer kmer_item, p_hc_apply apply)
{
    p_assemble_vertex vertex;

    vertex = assemble_graph_vertex_create(apply->read_thread_graph, kmer_item->seq, kmer_item->len, kmer_item);
    if (__glibc_unlikely(!vertex)) {
        statics_print_error("No Mem");
        exit(-1);
    }
    vertex = assemble_graph_add_vertex(apply->read_thread_graph, vertex);
    if (__glibc_unlikely(!vertex)) {
        statics_print_error("No Mem");
        exit(-1);
    }
    return vertex;
}

void hc_assemble_utils_simple_print_all_edge(p_hc_apply apply)
{
    p_assemble_edge edge_el;

    CDL_FOREACH2(apply->read_thread_graph->edge_list, edge_el, link.all_st.next)
    {
        int i;

        printf("edge weight %lf from seq ", edge_el->weight);
        if (edge_el->link.from.vertex) {
            for (i = 0; i < ((p_hc_assemble_graph_kmer)(edge_el->link.from.vertex->storage))->len; i++) {
                printf("%c", ((p_hc_assemble_graph_kmer)(edge_el->link.from.vertex->storage))->seq[i]);
            }
        }
        else {
            printf("null");
        }
        printf(" to seq ");
        if (edge_el->link.to.vertex) {
            for (i = 0; i < ((p_hc_assemble_graph_kmer)(edge_el->link.to.vertex->storage))->len; i++) {
                printf("%c", ((p_hc_assemble_graph_kmer)(edge_el->link.to.vertex->storage))->seq[i]);
            }
        }
        else {
            printf("null");
        }
        printf("\n");
    }
    puts("===========================================================");
}

void hc_assemble_utils_simple_print_all_vertex_edge(p_hc_apply apply)
{
    p_assemble_graph_hash_node vertex_st;
    p_assemble_vertex vertex;
    int i;

    CDL_FOREACH(apply->read_thread_graph->vertics_hash->node_list, vertex_st)
    {
        vertex = vertex_st->value;

        printf("seq: \"");
        for (i = 0; i < ((p_hc_assemble_graph_kmer)(vertex->storage))->len; i++) {
            printf("%c", ((p_hc_assemble_graph_kmer)(vertex->storage))->seq[i]);
        }
        printf("\", len %d, \n", vertex->data_len);
    }
    puts("===========================================================");
}

/**
 * @brief 通过cigar获取read length
 *
 * @param cigar         cigar
 * @param cigar_len     cigar长度
 *
 * @return int  read length
 */
int hc_assemble_utils_get_cigar_read_length(uint32_t* cigar, uint32_t cigar_len)
{
    int read_length = 0;
    uint32_t i = 0, op, op_len;

    for (; i < cigar_len; i++) {
        op = bam_cigar_op(cigar[i]);
        op_len = bam_cigar_oplen(cigar[i]);

        switch (op) {
            case WORKER_SIGAR_STATUS_MATCH:
            case WORKER_SIGAR_STATUS_INSERTION:
            case WORKER_SIGAR_SOFT_CLIP:
            case WORKER_SIGAR_S_MATCH:
            case WORKER_SIGAR_S_MISMATCH: read_length += op_len; break;
            default: break;
        }
    }
    return read_length;
}

/**
 * @return the reference source vertex pulled from the graph, can be null if it doesn't exist in the graph
 */
p_assemble_graph_hash_node hc_assemble_utils_get_reference_source_vertex_with_head(p_assemble_graph_hash_node head,
                                                                                   p_assemble_graph_hash_node end, p_assemble_graph graph)
{
    p_assemble_graph_hash_node vertex_node;
    p_assemble_vertex vertex;

    for (vertex_node = head; vertex_node; vertex_node = (vertex_node->next == end ? NULL : vertex_node->next)) {
        vertex = vertex_node->value;
        if (hc_assemble_base_graph_is_ref_source(vertex, graph)) {
            return vertex_node;
        }
    }

    return NULL;
}

/**
 * @return the reference source vertex pulled from the graph, can be null if it doesn't exist in the graph
 */
p_assemble_graph_hash_node hc_assemble_utils_get_reference_source_vertex(p_assemble_graph_hash_node head, p_assemble_graph_hash_node end,
                                                                         p_assemble_graph graph)
{
    p_assemble_graph_hash_node vertex_node;
    p_assemble_vertex vertex;

    if (!head->next || head->next == end) {
        return NULL;
    }
    vertex_node = head->next;

    for (; vertex_node; vertex_node = (vertex_node->next == end ? NULL : vertex_node->next)) {
        vertex = vertex_node->value;
        if (hc_assemble_base_graph_is_ref_source(vertex, graph)) {
            return vertex_node;
        }
    }

    return NULL;
}

p_assemble_graph_hash_node hc_assemble_utils_get_reference_sink_vertex_with_head(p_assemble_graph_hash_node head,
                                                                                 p_assemble_graph_hash_node end, p_assemble_graph graph)
{
    p_assemble_graph_hash_node vertex_node;
    p_assemble_vertex vertex;

    for (vertex_node = head; vertex_node; vertex_node = (vertex_node->next == end ? NULL : vertex_node->next)) {
        vertex = vertex_node->value;
        if (hc_assemble_base_graph_is_ref_sink(vertex, graph)) {
            return vertex_node;
        }
    }

    return NULL;
}

/**
 * @return the reference sink vertex pulled from the graph, can be null if it doesn't exist in the graph
 */
p_assemble_graph_hash_node hc_assemble_utils_get_reference_sink_vertex(p_assemble_graph_hash_node head, p_assemble_graph_hash_node end,
                                                                       p_assemble_graph graph)
{
    p_assemble_graph_hash_node vertex_node;
    p_assemble_vertex vertex;

    if (!head->next || head->next == end) {
        return NULL;
    }
    vertex_node = head->next;

    for (; vertex_node; vertex_node = (vertex_node->next == end ? NULL : vertex_node->next)) {
        vertex = vertex_node->value;
        if (hc_assemble_base_graph_is_ref_sink(vertex, graph)) {
            return vertex_node;
        }
    }

    return NULL;
}

/**
 * @brief 遍历Graph，Cycle时默认操作
 *
 * @param vertex    发生cycle时最后一个vertex
 * @param input     遍历储存
 *
 * @return int 1
 */
int hc_assemble_utils_default_graph_func_cycle(p_assemble_vertex vertex, p_assemble_graph_iter input)
{
    (void)vertex;
    (void)input;
    return 1;
}

/**
 * @brief 遍历Graph，一条路径完成时默认操作
 *
 * @param vertex    路径完成时最后一个vertex
 * @param input     遍历储存
 *
 * @return int 0
 */
int hc_assemble_utils_default_graph_func_one_path_done(p_assemble_vertex vertex, p_assemble_graph_iter input)
{
    (void)vertex;
    (void)input;
    return 0;
}

/**
 * @brief 遍历Graph，所有路径完成时默认操作
 *
 * @param vertex    路径完成时最后一个vertex
 * @param input     遍历储存
 *
 * @return int 0
 */
int hc_assemble_utils_default_graph_func_all_path_done(p_assemble_vertex vertex, p_assemble_graph_iter input)
{
    (void)vertex;
    (void)input;
    return 0;
}

/**
 * @brief 检查list可用
 *
 * @param list list
 *
 * @return p_hc_shared_lib_list_head list
 */
inline p_hc_shared_lib_list_head hc_assemble_utils_check_list_work(p_hc_shared_lib_list_head list)
{
    if (list->nodes) {
        hc_shared_list_reset(list);
    }
    return list;
}

/**
 * @brief 检查list可用
 *
 * @param list list
 *
 * @return p_hc_shared_lib_list_head list
 */
inline p_hc_shared_lib_pair_list_head hc_assemble_utils_check_pair_list_work(p_hc_shared_lib_pair_list_head list)
{
    if (list->nodes) {
        hc_shared_pair_list_reset(list);
    }
    return list;
}
/**
 * @brief 检查hash可用
 *
 * @param hash hash
 *
 * @return p_assemble_graph_hash_table hash
 */
inline p_assemble_graph_hash_table hc_assemble_utils_check_hash_work(p_assemble_graph_hash_table hash)
{
    if (hash->node_list) {
        assemble_graph_hash_reset(hash);
    }
    return hash;
}

/**
 * @brief Seq Graph插入相同vertex时Hash Table执行函数
 *
 * @param old   old data
 * @param new   new data
 * @return int  0
 */
static int hc_assemble_utils_hash_equal_func(void* old, void* new)
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
void hc_assemble_utils_hash_new_node_memcpy(void* new_node, void* input_pt, int hash_len, p_assemble_graph_hash_table hash_table)
{
    (void)hash_table;
    p_assemble_graph_hash_node insert_node = new_node;

    memset(insert_node, 0, sizeof(assemble_graph_hash_node));

    insert_node->key = input_pt;
    insert_node->keylen = hash_len;
    insert_node->value = input_pt;
}

void hc_assemble_utils_hash_pt_new_node_memcpy(void* new_node, void* input_pt, int hash_len, p_assemble_graph_hash_table hash_table)
{
    (void)hash_table;
    p_assemble_graph_hash_node insert_node = new_node;
    void* key = insert_node->key;

    memset(insert_node, 0, sizeof(assemble_graph_hash_node));

    insert_node->key = key;
    insert_node->keylen = hash_len;
    insert_node->value = input_pt;
}

#define mem_size_t unsigned long long
#define MB         (mem_size_t)(1 << 20)

/**
 * @brief 初始化read graph
 *
 * @param graph graph
 *
 * @return int 1:succeed, 0:failed
 */
int hc_assemble_utils_init_basic_read_graph(p_hc_assemble_read_graph graph)
{
    graph->kmer_buffer = ring_mempool_auto_enlarge_init(sizeof(hc_assemble_graph_kmer), HC_ASSEMBLE_GRAPH_POOL_PER_REGION);
    graph->item_buffer = ring_mempool_auto_enlarge_init(sizeof(hc_assemble_graph_read), HC_ASSEMBLE_GRAPH_POOL_PER_REGION);
    graph->list_buffer_1 = hc_shared_list_init(HC_ASSEMBLE_GRAPH_POOL_PER_REGION);
    graph->dijkstra_path_finder.result = hc_shared_list_init(HC_ASSEMBLE_GRAPH_POOL_PER_REGION);
    graph->read_hash =
        assemble_graph_hash_init(ASSEMBLE_GRAPH_HASH_NODES, hc_assemble_utils_hash_equal_func, hc_assemble_utils_hash_new_node_memcpy);
    graph->hash_buffer_1 =
        assemble_graph_hash_init(ASSEMBLE_GRAPH_HASH_NODES, hc_assemble_utils_hash_equal_func, hc_assemble_utils_hash_new_node_memcpy);
    graph->hash_buffer_2 =
        assemble_graph_hash_init(ASSEMBLE_GRAPH_HASH_NODES, hc_assemble_utils_hash_equal_func, hc_assemble_utils_hash_new_node_memcpy);
    graph->list_buffer_2 = hc_shared_list_init(HC_ASSEMBLE_GRAPH_POOL_PER_REGION);
    graph->cache_pool = hc_assemble_thread_mem_init(8 * MB);
    graph->vertex_spliter.seq_graph = assemble_graph_create();
    graph->vertex_spliter.middles = hc_shared_pair_list_init(HC_ASSEMBLE_GRAPH_POOL_PER_REGION);
    graph->vertex_spliter.edges_to_remove = hc_shared_list_init(HC_ASSEMBLE_GRAPH_POOL_PER_REGION);
    graph->vertex_spliter.all_vertex = hc_shared_list_init(HC_ASSEMBLE_GRAPH_POOL_PER_REGION);

    if (__glibc_likely(graph->read_hash && graph->kmer_buffer && graph->item_buffer && graph->list_buffer_1 && graph->hash_buffer_1 &&
                       graph->hash_buffer_2 && graph->list_buffer_2 && graph->cache_pool && graph->vertex_spliter.seq_graph &&
                       graph->vertex_spliter.middles && graph->vertex_spliter.edges_to_remove && graph->vertex_spliter.all_vertex)) {
        assemble_graph_to_seq_graph(graph->vertex_spliter.seq_graph);
        return 1;
    }
    else {
        return 0;
    }
}

/**
 * @brief 释放read graph
 *
 * @param graph graph
 */
void hc_assemble_utils_finit_basic_read_graph(p_hc_assemble_read_graph graph)
{
    ring_mempool_auto_enlarge_term(graph->kmer_buffer);
    ring_mempool_auto_enlarge_term(graph->item_buffer);
    hc_shared_list_term(graph->list_buffer_1);
    hc_shared_list_term(graph->dijkstra_path_finder.result);
    assemble_graph_hash_destroy(graph->read_hash);
    assemble_graph_hash_destroy(graph->hash_buffer_1);
    assemble_graph_hash_destroy(graph->hash_buffer_2);
    hc_shared_list_term(graph->list_buffer_2);
    hc_assemble_thread_mem_term(graph->cache_pool);
    assemble_graph_free(graph->vertex_spliter.seq_graph);
    hc_shared_pair_list_term(graph->vertex_spliter.middles);
    hc_shared_list_term(graph->vertex_spliter.edges_to_remove);
    hc_shared_list_term(graph->vertex_spliter.all_vertex);
}

/**
 * @brief Seq Graph创建一个新的vertex，基于输入vertex复制
 *
 * @param vertex        原有vertex
 * @param data_len      len
 * @param apply         apply
 *
 * @return p_assemble_vertex vertex
 */
p_assemble_vertex hc_assemble_utils_create_seq_vertex_copy_from_seq_vertex(p_assemble_vertex vertex, int data_len, p_hc_apply apply)
{
    p_hc_assemble_graph_kmer kmer;
    uint8_t* seq;

    kmer = ring_mempool_auto_enlarge_malloc(apply->read_graph.kmer_buffer);
    memset(kmer, 0, sizeof(hc_assemble_graph_kmer));
    seq = hc_assemble_thread_mem_malloc(apply->read_graph.cache_pool, data_len + 1);
    memcpy(seq, vertex->data, data_len);
    seq[data_len] = 0;

    kmer->seq = seq;
    kmer->len = data_len;
    kmer->item = ((p_hc_assemble_graph_kmer)(vertex->storage))->item;

    vertex = assemble_graph_vertex_create(apply->read_thread_graph, kmer->seq, kmer->len, kmer);

    return vertex;
}

/**
 * @brief Seq Graph创建一个新的vertex，基于输入data复制
 *
 * @param data          data
 * @param data_len      data len
 * @param apply         apply
 *
 * @return p_assemble_vertex    vertex
 */
p_assemble_vertex hc_assemble_utils_create_seq_vertex_from_data(uint8_t* data, int data_len, p_hc_apply apply)
{
    p_hc_assemble_graph_kmer kmer;
    p_assemble_vertex vertex;
    uint8_t* seq;

    kmer = ring_mempool_auto_enlarge_malloc(apply->read_graph.kmer_buffer);
    memset(kmer, 0, sizeof(hc_assemble_graph_kmer));
    seq = hc_assemble_thread_mem_malloc(apply->read_graph.cache_pool, data_len + 1);
    memcpy(seq, data, data_len);
    seq[data_len] = 0;

    kmer->seq = seq;
    kmer->len = data_len;
    kmer->item = NULL;

    vertex = assemble_graph_vertex_create(apply->read_thread_graph, kmer->seq, kmer->len, kmer);

    return vertex;
}

/**
 * @brief Seq Graph创建一个新的vertex，基于输入两个vertex复制
 *
 * @param vertex_a  vertex a
 * @param vertex_b  vertex b
 * @param apply     apply
 *
 * @return p_assemble_vertex vertex
 */
p_assemble_vertex hc_assemble_utils_merge_seq_vertex_from_two_seq_vertex(p_assemble_vertex vertex_a, p_assemble_vertex vertex_b,
                                                                         p_hc_apply apply)
{
    p_hc_assemble_graph_kmer kmer;
    p_assemble_vertex vertex;
    uint8_t* seq;

    // kmer
    kmer = ring_mempool_auto_enlarge_malloc(apply->read_graph.kmer_buffer);
    memset(kmer, 0, sizeof(hc_assemble_graph_kmer));
    kmer->len = vertex_a->data_len + vertex_b->data_len;
    kmer->item = NULL;

    // seq
    seq = hc_assemble_thread_mem_malloc(apply->read_graph.cache_pool, kmer->len);
    if (__glibc_unlikely(!seq)) {
        return NULL;
    }
    kmer->seq = seq;
    memcpy(seq, vertex_a->data, vertex_a->data_len);
    seq += vertex_a->data_len;
    memcpy(seq, vertex_b->data, vertex_b->data_len);

    vertex = assemble_graph_vertex_create(apply->read_thread_graph, kmer->seq, kmer->len, kmer);

    return vertex;
}

/**
 * @brief 倒置一个stack
 *
 * @param stack         stack
 * @param stack_len     stack长度
 */
void hc_assemble_utils_reverse_stack(uint32_t* stack, uint32_t stack_len)
{
    uint32_t i = 0, middle = stack_len >> 1, tmp;

    for (; i < middle; i++) {
        tmp = stack[i];
        stack[i] = stack[stack_len - 1 - i];
        stack[stack_len - 1 - i] = tmp;
    }
}

/**
 * @brief list 获取index值的item
 *
 * @param list      list
 * @param index     index
 *
 * @return void*    item
 */
void* hc_assemble_utils_get_list_index_item(p_hc_shared_lib_list_head list, unsigned int index)
{
    register unsigned int i = 0;
    p_hc_shared_lib_list_item item = list->nodes;

    CDL_FOREACH(list->nodes, item)
    {
        if (i == index) {
            return item->data;
        }
        i += 1;
    }
    return NULL;
}

#include "rbtree.h"

void hc_assemble_utils_print_read_graph(p_hc_apply apply, const char* path)
{
    struct rb_node* node;
    p_hc_apply_one_read el;
    uint32_t i;
    FILE* fp = fopen(path, "w");

    if (!fp) {
        return;
    }

    for (node = rb_first(&apply->read_sorted); node; node = rb_next(node)) {
        el = rb_entry(node, hc_apply_one_read, node);
        fprintf(fp, "%s %ld %.*s ", el->read_data, el->pos_start + 1, el->read_len, el->seq);
        for (i = 0; i < el->read_len; i++) {
            fprintf(fp, "%u ", el->qual[i]);
        }
        fputs("\n", fp);
    }
    fclose(fp);
}

static inline int hc_assemble_utils_compare_int(int a, int b)
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

int hc_assemble_utils_compare_coordinates(p_hc_apply_one_read first, p_hc_apply_one_read second)
{
    int firstRefIndex = first->read.core.tid;
    int secondRefIndex = second->read.core.tid;

    if (firstRefIndex == -1) {
        return (secondRefIndex == -1 ? 0 : 1);
    }
    else if (secondRefIndex == -1) {
        return -1;
    }

    int refIndexDifference = firstRefIndex - secondRefIndex;
    if (refIndexDifference != 0) {
        return refIndexDifference;
    }

    return hc_assemble_utils_compare_int(first->pos_start, second->pos_start);
}