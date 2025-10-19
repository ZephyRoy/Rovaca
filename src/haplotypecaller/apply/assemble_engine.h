#ifndef ASSEMBLE_ENGINE_H
#define ASSEMBLE_ENGINE_H

#include <thread>
#include <vector>

#include "assemble_argument.h"
#include "assemble_result.h"
#include "hc_assemble_main.h"

using OneNodeData = p_assemble_dijkstra_one_path_node;
static const int k_assemble_cigar_string_mem = 128;
static const int k_assemble_cigar_cache_mem = 128;
static const int k_reads_buffer_capacity = 1024 * 200;  // 每条read需要约616b, 约300条reads.
class AssembleEngine
{
public:
    /**
     * @brief Assemble主接口
     * @param region active span.
     * @param ref active span padded with REF_PADDING(500).
     * @param ref_len ref + REF_PADDING*2
     * @param region_reads 根据active span所取的reads.
     * @param mem memory resource for all necessary results.
     * @return AssembleResult*
     */
    static AssembleResult* local_assemble(p_hc_apply assembler, p_hc_region_active_storage region, const uint8_t* chr_ref,
                                          uint32_t chr_ref_len, const std::vector<bam1_t*>& region_reads, pMemoryPool mem,
                                          AssembleReadsBuffer* reads_mem);
    /**
     * @brief 根据配置信息初始化Assemble模块参数。
     * @param assemble_argument
     */
    static void init_assemble_argument(const AssembleArgument* assemble_argument);

    /**
     * @brief 反初始化assemble参数。
     */
    static void finit_assemble_argument();

    /**
     * @brief Create a assembler instance.
     * @return p_hc_apply
     */
    static p_hc_apply creat_assembler();

    /**
     * @brief destory a assembler instance.
     * @param assembler
     */
    static void destory_assembler(p_hc_apply assembler);

private:
    /**
     * @brief 从Apply结构体中提取reads，可直接供后边模块调用。
     * @param assembler
     * @param reads
     * @param mem
     * @return true
     * @return false
     */
    static bool extract_reads(p_hc_apply assembler, std::pmr::vector<pReadRecord>& reads, pMemoryPool mem);

    /**
     * @brief 从Apply结构体中提取单倍体，可直接供后边模块调用。
     * @param assembler
     * @param haplotypes
     * @param mem
     * @return true
     * @return false
     */
    static bool extract_haplotypes(p_hc_apply assembler, std::pmr::vector<pHaplotype>& haplotypes, pSimpleInterval region, pMemoryPool mem);

    /**
     * @brief Create a one assemble result object.
     * 生命周期与reads，haplotypes相同，从一个memory pool中创建。
     * @param resource
     * @return AssembleResult*
     */
    static AssembleResult* create_one_assemble_result(pMemoryPool resource);
    static void destory_one_assemble_result(AssembleResult* result);

    /**
     * @brief Create a one assemble read object.
     * 为apply read中的cigar提供缓存空间，每个region需要约200k空间，若空间不够采用malloc分配。
     * @param read
     * @param reads_mem
     * @return p_hc_apply_one_read
     */
    static p_hc_apply_one_read create_one_assemble_read(bam1_t* read, pAssembleReadsBuffer reads_mem);

    /** 根据mempolicy销毁内存。
     * @brief
     * @param read
     * @param reads_mem
     */
    static void destory_one_assemble_read(p_hc_apply_one_read read);

    /**
     * @brief 压缩实际使用内存
     * @param read
     * @param reads_mem
     */
    static void adjust_assemble_read_mem(p_hc_apply_one_read read, pAssembleReadsBuffer reads_mem);

    /**
     * @brief 从Apply结构体中提取cigar.。
     * @param apply_node
     * @param mem
     * @return pCigar
     */
    static pCigar build_cigar_from_apply(OneNodeData apply_node, pMemoryPool mem);

    /**
     * @brief 从Apply结构体中提取bases.
     * @param apply_node
     * @param mem
     * @return pBases
     */
    static pBases build_bases_from_apply(OneNodeData apply_node, pMemoryPool mem);
};

#endif  // ASSEMBLE_ENGINE_H