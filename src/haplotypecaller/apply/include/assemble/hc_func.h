#ifndef __SHARED_HC_FUNC_H__
#define __SHARED_HC_FUNC_H__

#include "assemble_interface.h"
#include "graph.h"
#include "hash.h"
#include "hc_apply.h"
#include "hc_assemble_graph.h"

// p_hc_main_storage hc_region_init_main(void);
// void hc_region_finit_main(p_hc_main_storage main_st);
// void hc_set_target_args(p_hc_region region,int force_bed , hts_pos_t *start , hts_pos_t *end , int regcount);
// void hc_set_target_len(p_hc_region region, int target_len);
// int hc_region_main(p_hc_main_storage main_storage, const char* ref);
// int hc_region_run_bases_active(p_hc_region region_storage);
// // interface for m2 module
// void flush_active_region(p_hc_main_storage main_storage, p_hc_region_base_active active_base, p_hc_region_active out_region);
// int hc_region_flush_storage(p_hc_main_storage main_storage,hts_pos_t ref_len);
// // interface for m2 module
// void hc_region_active_base_to_active_region(p_hc_region_base_active active_base, p_hc_region_active out_region,uint32_t items,uint32_t
// current_pos,uint32_t force); void flush_active_region_core(p_hc_region_base_active active_base,p_hc_region_active out_region);

/**
 * hc_region_utils.c
 *
 * hc_consumes_read_bases                       Bases消耗Seq长度
 * consumes_ref_bases                  Bases消耗Ref长度
 * hc_utils_count_high_quality_soft_clips       计算Read内高质量Soft总量
 * hc_region_utils_qual_to_prob                 将按比例缩放的质量得分转换为其真实的log10概率
 * hc_region_utils_approximate_log10_sum_log10  近似Log10和Log10的总和
 * hc_region_utils_add_hq_soft_clips            增加高质量Soft Clip计数
 * hc_region_utils_log10_binomial_coefficient   log10应用于结果
 * hc_region_utils_max_element_index            数组最大数位置
 * hc_region_utils_log_to_log10                 log->log10
 * hc_region_utils_normalize_log10              归一化log10数组
 * hc_region_utils_init_gaussian_kernel         缓存各数值下高斯核结果
 */

// unsigned int hc_utils_count_high_quality_soft_clips(bam1_t* bam_info);
// double hc_region_utils_qual_to_prob(uint16_t qual);
// double hc_region_utils_qual_to_error_prob_log10(uint8_t qual);
// double hc_region_utils_approximate_log10_sum_log10(double a, double b, p_hc_region_math_cache math_st);
// void hc_region_utils_add_hq_soft_clips(p_hc_region_pos_cacl_cache result, double obs);
// double hc_region_utils_log10_binomial_coefficient(int n, int k);
// int hc_region_utils_max_element_index(double* array, int sum);
// double hc_region_utils_log_to_log10(double ln);
// void hc_region_utils_normalize_log10(double* array, int sum);
// void hc_region_utils_init_gaussian_kernel(p_hc_region_math_cache cache);
// void hc_region_utils_init_math(p_hc_region_math_cache cache);

/**
 * hc_region_reference_confidence.c
 *
 * hc_region_rc_calc_genotype_likelihoods       HC计算Genotype Likelihoods
 */
// void hc_region_rc_calc_genotype_likelihoods(p_hc_region_pos_cacl_cache item,
//                                             p_ring_mempool out_pool,
//                                             p_hc_region_math_cache math_st);

// void hc_region_his_rc_calc_genotype_likelihoods(p_hc_region_his_pos_cacl_cache item,
//                                                 p_hc_region_math_cache math_st);
/**
 * hc_assemble_reads.c
 *
 * hc_assemble_reads                            HC 组装
 */
void hc_assemble_reads(p_hc_apply apply);
int hc_assemble_finalize_region(p_hc_apply apply, p_hc_apply_one_read read);
/**
 * hc_assemble_utils.c
 *
 * hc_assemble_utils_get_new_cigar                  获取cigar
 * hc_assemble_utils_revert_soft_clip               取消SOFT CLIP
 * hc_assemble_utils_copy_seq                       复制bam序列
 * hc_assemble_utils_soft_clip_to_hard              CIGAR，SOFT转为HARD
 * hc_assemble_utils_clip_low_qual_ends        HARD CLIP掉过低的头尾
 * hc_assemble_utils_read_is_unmapped               Read是否为Unmapped
 * hc_assemble_utils_init_hc_read                   初始化Read结构
 * hc_assemble_utils_hc_read_to_cache               使用Cache初始化Read结构
 * hc_assemble_utils_clip_read                 按照头尾Hard Clip掉头尾
 * hc_assemble_utils_hard_clip_to_region            依据region进行hard clip
 * hc_assemble_utils_fill_kmer_by_read              通过read及kmer start生产kmer
 * hc_assemble_utils_add_graph_vertex_by_kmer       通过Kmer建立Vertex并加入graph
 * hc_assemble_utils_get_cigar_read_length          通过cigar获取read length
 * hc_assemble_utils_get_reference_source_vertex    reference source vertex pulled from the graph
 * hc_assemble_utils_get_reference_sink_vertex      the reference sink vertex pulled from the graph
 * hc_assemble_utils_default_graph_func_cycle       遍历Graph，Cycle时默认操作
 * hc_assemble_utils_default_graph_func_one_path_done   遍历Graph，一条路径完成时默认操作
 * hc_assemble_utils_default_graph_func_all_path_done   遍历Graph，所有路径完成时默认操作
 * hc_assemble_utils_check_list_work                检查list可用
 * hc_assemble_utils_check_pair_list_work           检查list可用
 * hc_assemble_utils_check_hash_work                检查hash可用
 * hc_assemble_utils_init_basic_read_graph          初始化read graph
 * hc_assemble_utils_finit_basic_read_graph         释放read graph
 * hc_assemble_utils_create_seq_vertex_copy_from_seq_vertex     Seq Graph创建一个新的vertex，基于输入vertex复制
 * hc_assemble_utils_create_seq_vertex_from_data                Seq Graph创建一个新的vertex，基于输入data复制
 * hc_assemble_utils_merge_seq_vertex_from_two_seq_vertex       Seq Graph创建一个新的vertex，基于输入两个vertex复制
 * hc_assemble_utils_reverse_stack                  倒置一个stack
 * hc_assemble_graph_layout                         输出图
 * hc_assemble_utils_get_list_index_item            list 获取index值的item
 * hc_assemble_utils_get_ref_length                             cigar获取ref长度
 * hc_assemble_utils_compare_coordinates
 *
 */
inline uint32_t hc_assemble_utils_get_new_cigar(uint8_t cigar_op, uint32_t cigar_len);
void hc_assemble_utils_revert_soft_clip(p_hc_apply_one_read read);
void hc_assemble_utils_copy_seq(uint8_t* __restrict out, uint8_t* __restrict in, uint32_t start, uint32_t len);
void hc_assemble_utils_soft_clip_to_hard(p_hc_apply_one_read read);
void hc_assemble_utils_clip_low_qual_ends(p_hc_apply_one_read read, clipping_type clip_type);
int hc_assemble_utils_read_is_unmapped(p_hc_apply_one_read read);
void hc_assemble_utils_init_hc_read(p_hc_apply_one_read read);
void hc_assemble_utils_hc_read_to_cache(p_hc_apply_one_read read);
void hc_assemble_utils_clip_read(p_hc_apply_one_read read, uint32_t clip_start, uint32_t clip_end, clipping_type clip_type);
void hc_assemble_utils_hard_clip_to_region(p_hc_apply_one_read read, int refStart, int refStop);
p_hc_assemble_graph_kmer hc_assemble_utils_fill_kmer_by_read(p_hc_assemble_graph_read read, int kmerStart, p_hc_apply apply);
p_assemble_vertex hc_assemble_utils_add_graph_vertex_by_kmer(p_hc_assemble_graph_kmer kmer_item, p_hc_apply apply);
void hc_assemble_utils_simple_print_all_edge(p_hc_apply apply);
void hc_assemble_utils_simple_print_all_vertex_edge(p_hc_apply apply);
int hc_assemble_utils_get_cigar_read_length(uint32_t* cigar, uint32_t cigar_len);
p_assemble_graph_hash_node hc_assemble_utils_get_reference_source_vertex_with_head(p_assemble_graph_hash_node head,
                                                                                   p_assemble_graph_hash_node end, p_assemble_graph graph);
p_assemble_graph_hash_node hc_assemble_utils_get_reference_source_vertex(p_assemble_graph_hash_node head, p_assemble_graph_hash_node end,
                                                                         p_assemble_graph graph);
p_assemble_graph_hash_node hc_assemble_utils_get_reference_sink_vertex_with_head(p_assemble_graph_hash_node head,
                                                                                 p_assemble_graph_hash_node end, p_assemble_graph graph);
p_assemble_graph_hash_node hc_assemble_utils_get_reference_sink_vertex(p_assemble_graph_hash_node head, p_assemble_graph_hash_node end,
                                                                       p_assemble_graph graph);
int hc_assemble_utils_default_graph_func_cycle(p_assemble_vertex vertex, p_assemble_graph_iter input);
int hc_assemble_utils_default_graph_func_one_path_done(p_assemble_vertex vertex, p_assemble_graph_iter input);
int hc_assemble_utils_default_graph_func_all_path_done(p_assemble_vertex vertex, p_assemble_graph_iter input);
inline p_hc_shared_lib_list_head hc_assemble_utils_check_list_work(p_hc_shared_lib_list_head list);
inline p_hc_shared_lib_pair_list_head hc_assemble_utils_check_pair_list_work(p_hc_shared_lib_pair_list_head list);
inline p_assemble_graph_hash_table hc_assemble_utils_check_hash_work(p_assemble_graph_hash_table hash);
int hc_assemble_utils_init_basic_read_graph(p_hc_assemble_read_graph graph);
void hc_assemble_utils_finit_basic_read_graph(p_hc_assemble_read_graph graph);
p_assemble_vertex hc_assemble_utils_create_seq_vertex_copy_from_seq_vertex(p_assemble_vertex vertex, int data_len, p_hc_apply apply);
p_assemble_vertex hc_assemble_utils_create_seq_vertex_from_data(uint8_t* data, int data_len, p_hc_apply apply);
p_assemble_vertex hc_assemble_utils_merge_seq_vertex_from_two_seq_vertex(p_assemble_vertex vertex_a, p_assemble_vertex vertex_b,
                                                                         p_hc_apply apply);

void hc_assemble_utils_reverse_stack(uint32_t* stack, uint32_t stack_len);
void* hc_assemble_utils_get_list_index_item(p_hc_shared_lib_list_head list, unsigned int index);
int hc_assemble_utils_get_ref_length(uint32_t* cigar, uint32_t cigar_len);
int hc_assemble_utils_compare_coordinates(p_hc_apply_one_read first, p_hc_apply_one_read second);
void hc_assemble_utils_hash_new_node_memcpy(void* new_node, void* input_pt, int hash_len, p_assemble_graph_hash_table hash_table);
void hc_assemble_utils_hash_pt_new_node_memcpy(void* new_node, void* input_pt, int hash_len, p_assemble_graph_hash_table hash_table);

/**
 * hc_apply_main.c
 *
 * hc_apply_main                                    HC运算Region, 主函数
 * hc_apply_reset                                   HC运算Region，重置
 * hc_apply_init                                    HC运算Region, 初始化
 * hc_apply_finit                                   HC运算Region, 反初始化
 */
void hc_apply_main(p_hc_apply apply);
void hc_apply_reset(p_hc_apply apply);
p_hc_apply hc_apply_init();
int hc_apply_set_target(p_hc_apply apply, p_hc_region_active_storage region, uint8_t* ref, uint32_t ref_len);
int hc_apply_set_region_reads(p_hc_apply apply, p_hc_apply_one_read run_read);
void hc_apply_finit(p_hc_apply apply);

/**
 * hc_apply_utils.c
 *
 * hc_apply_utils_mate_read_is_unmapped             mate read 是否unmapped
 * hc_apply_utils_has_well_defined_fragment_size    是否可以基于读取及其配对的对齐方式，从读取中可靠地除去读取的衔接子序列
 * hc_apply_utils_hard_clip_adaptor_sequence        HC Assemble 重建缓存read
 * hc_apply_utils_get_soft_start                    获取soft start
 * hc_assemble_utils_adjust_overlapping_paired_qual 如果可能，通过调整碱基质量来修复同一片段的两个重叠读段
 * hc_apply_alloc_read                              申请一条read空间
 * hc_apply_utils_clip_by_reference_coordinates     依据ref裁剪read
 */
inline int hc_apply_utils_mate_read_is_unmapped(p_hc_apply_one_read read);
void hc_apply_utils_hard_clip_adaptor_sequence(p_hc_apply_one_read read);
int hc_apply_utils_has_well_defined_fragment_size(p_hc_apply_one_read read);
int hc_apply_utils_get_soft_start(p_hc_apply_one_read read);
void hc_assemble_utils_adjust_overlapping_paired_qual(p_hc_apply_one_read read1, p_hc_apply_one_read read2);
void hc_apply_utils_clip_by_reference_coordinates(int refStart, int refStop, p_hc_apply_one_read read);
void hc_apply_utils_init_region(p_hc_apply apply, p_hc_region_active_storage region, const uint8_t* chr_ref, uint32_t chr_ref_length);

/**
 * hc_assemble_read_threading_assembler.c
 *
 * hc_assemble_graph_build                          Build Graph
 */
int hc_assemble_graph_build(p_hc_apply apply);

/**
 * hc_assemble_read_threading_graph.c
 *
 * hc_assemble_graph_determine_ref_non_unique_kmers     检测Ref Unique Kmer
 * hc_assemble_graph_read_determine_non_unique_kmers    检测Read Unique Kmer
 * hc_assemble_graph_thread_sequence                    Thread sequence seqForKmers through the current graph, updating the graph as
 * appropriate hc_assemble_graph_get_suffix                         Get the suffix byte of this DeBruijnVertex
 */
int hc_assemble_graph_determine_ref_non_unique_kmers(p_hc_apply apply);
void hc_assemble_graph_read_determine_non_unique_kmers(p_hc_assemble_graph_read read, p_hc_apply apply);
int hc_assemble_graph_thread_sequence(p_hc_apply apply);
uint8_t hc_assemble_graph_get_suffix(p_assemble_vertex vertex);

/**
 * hc_assemble_base_graph.c
 *
 * hc_assemble_base_graph_cycle_finder                  搜索Graph内Cycle
 * hc_assemble_base_graph_recover_dangling_tails        Try to recover dangling tails
 * hc_assemble_base_graph_recover_dangling_haids        Try to recover dangling heads
 * hc_assemble_base_graph_is_ref_sink                   true if this vertex is a reference sink
 * hc_assemble_base_graph_is_ref_source                 true if this vertex is a reference source
 */
int hc_assemble_base_graph_cycle_finder(p_hc_apply apply);
void hc_assemble_base_graph_recover_dangling_tails(p_hc_apply apply);
void hc_assemble_base_graph_recover_dangling_heads(p_hc_apply apply);
int hc_assemble_base_graph_is_ref_sink(p_assemble_vertex vertex, p_assemble_graph graph);
int hc_assemble_base_graph_is_ref_source(p_assemble_vertex vertex, p_assemble_graph graph);

/**
 * hc_assemble_seq_graph.c
 *
 * hc_assemble_seq_graph_to_sequence_graph              Convert this kmer graph to a simple sequence graph
 * hc_assemble_seq_graph_zip_linear_chains              Zip up all of the simple linear chains present in this graph
 * hc_assemble_seq_graph_remove_singleton_orphan_vertices   Remove all vertices in the graph that have in and out degree of 0
 * hc_assemble_seq_graph_simplify_graph                 Simplify this graph
 * hc_assemble_seq_graph_sanity_check_reference_graph    Make sure the reference sequence is properly represented in the provided graph
 */
void hc_assemble_seq_graph_to_sequence_graph(p_assemble_graph graph);
int hc_assemble_seq_graph_zip_linear_chains(p_hc_apply apply);
void hc_assemble_seq_graph_remove_singleton_orphan_vertices(p_hc_apply apply);
void hc_assemble_seq_graph_simplify_graph(p_hc_apply apply);
int hc_assemble_seq_graph_sanity_check_reference_graph(p_hc_apply apply);

/**
 * hc_assemble_base_graph_iter.c
 *
 * hc_assemble_base_graph_iter_init                     初始化iter资源
 * hc_assemble_base_graph_iter_finit                    释放iter资源
 * hc_assemble_base_graph_iter_dfs_non_recursion        深度优先遍历，无迭代
 * hc_assemble_base_graph_iter_dfs_non_recursion_up     深度优先遍历，无迭代，向上
 * hc_assemble_seq_graph_remove_vertices_not_connected_to_ref_regardless_of_edge_direction  Remove all vertices on the graph that cannot be
 * accessed by following any edge, regardless of its direction
 */
int hc_assemble_base_graph_iter_init(p_assemble_graph_iter input);
void hc_assemble_base_graph_iter_finit(p_assemble_graph_iter input);
int hc_assemble_base_graph_iter_dfs_non_recursion(p_assemble_graph_iter input);
int hc_assemble_base_graph_iter_dfs_non_recursion_up(p_assemble_graph_iter input);
void hc_assemble_seq_graph_remove_vertices_not_connected_to_ref_regardless_of_edge_direction(p_hc_apply apply);

/**
 * hc_assemble_vertex_sequence_spliter.c
 *
 * hc_assemble_vertex_sequence_spliter_merge_diamonds_try_to_transform          Merge diamond configurations
 * hc_assemble_vertex_sequence_spliter_merge_tails_try_to_transform             Merge tail configurations
 * hc_assemble_vertex_sequence_spliter_split_common_suffices_try_to_transform   Performs the transformation
 * hc_assemble_vertex_sequence_spliter_split_merge_common_suffices             Merge headless configurations
 */
int hc_assemble_vertex_sequence_spliter_merge_diamonds_try_to_transform(p_assemble_vertex top, p_hc_apply apply);
int hc_assemble_vertex_sequence_spliter_merge_tails_try_to_transform(p_assemble_vertex top, p_hc_apply apply);
int hc_assemble_vertex_sequence_spliter_split_common_suffices_try_to_transform(p_assemble_vertex bottom, p_hc_apply apply);
int hc_assemble_vertex_sequence_spliter_split_merge_common_suffices(p_assemble_vertex vertex, p_hc_apply apply);

/**
 * hc_assemble_dijkstra_shortest_path.c
 *
 * hc_assemble_dijkstra_find_best_haplotypes                                    Implement Dijkstra's algorithm as described
 * hc_assemble_dijkstra_reset_result                                            reset dijkstra result
 */
void hc_assemble_dijkstra_find_best_haplotypes(p_assemble_vertex source, p_assemble_vertex sink, p_hc_apply apply);
void hc_assemble_dijkstra_reset_result(p_assemble_dijkstra_path graph);

/**
 * hc_assemble_cigar_cacl.c
 *
 * hc_assemle_cigar_cacl_calculate_cigar                                        Calculate the cigar elements for this path against the
 * reference sequence.
 */
void hc_assemle_cigar_cacl_calculate_cigar(p_hc_apply apply, p_assemble_dijkstra_one_path_node altSeq);

/**
 * hc_assemble_seq_path_finder.c
 *
 * hc_assemble_seq_finder_find_best_paths                                       Find discover paths by using KBestHaplotypeFinder over each
 * graph.
 */
void hc_assemble_seq_finder_find_best_paths(p_hc_apply apply);

/**
 * hc_assemble_gatk_sw.c
 *
 * hc_assemble_gatk_sw_align                                                    Aligns the alternate sequence to the reference sequence
 */
int hc_assemble_gatk_sw_align(p_gatk_sw_storage align);
void hc_assemble_gatk_sw_init(p_gatk_sw_storage align);

/**
 * hc_assemble_chain_pruner.c
 *
 * void hc_assemble_chain_pruner_prune_low_weight_chains(p_hc_apply apply)      使用chain做裁剪
 */
void hc_assemble_chain_pruner_prune_low_weight_chains(p_hc_apply apply);

/**
 * hc_assemble_random.c
 *
 * hc_assemble_random_next_int                                                  random int
 */
// int32_t hc_assemble_random_next_int(int32_t bound, uint64_t* seed);

/**
 * hc_down_sample.c
 *
 * int hc_down_sample(p_hc_apply apply)                                         down sample
 */

// int hc_down_sample(p_hc_apply apply);

#endif  // !__SHARED_HC_FUNC_H__