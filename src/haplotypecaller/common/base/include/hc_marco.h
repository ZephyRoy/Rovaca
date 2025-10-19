#ifndef __SHARED_HC_MAIN_H__
#define __SHARED_HC_MAIN_H__

#ifdef __cplusplus
extern "C"
{
#endif

#include "htslib/sam.h"

#define HC_APPLY_SUM_MAX (64)
#define BAM1_INIT_LEN    (2048)
#define HC_R_NAME_MAX    (4096)
#define BQSR_FASTA_SEQ   "=ACMGRSVTWYHKDBN"
// 固定长度len自动优化
#define BQSR_FASTA_SEQ_LEN (strlen(BQSR_FASTA_SEQ))

#define HTSLIB_SAM_H

    // sam内cigar状态
    enum WORKER_CIGAR_PARSE_STATUS {
        WORKER_SIGAR_STATUS_MATCH = 0,      // M 匹配上了
        WORKER_SIGAR_STATUS_INSERTION = 1,  // I 插入点，loop这次循环
        WORKER_SIGAR_STATUS_DELETION = 2,   // D 缺少点，跳到下一个slot写数据
        WORKER_SIGAR_STATUS_NSTATUS = 3,    // N 屏蔽点，不增加 = D
        WORKER_SIGAR_SOFT_CLIP = 4,         // S 屏蔽点，不增加
        WORKER_SIGAR_HARD_CLIP = 5,         // H 屏蔽点，不增加
        WORKER_SIGAR_PADDING = 6,           // P 屏蔽点，不增加
        WORKER_SIGAR_S_MATCH = 7,           // = = M
        WORKER_SIGAR_S_MISMATCH = 8,        // X = M
        WORKER_SIGAR_MAX = 9,
        WORKER_SIGAR_SPECIAL_READ = 10,
        WORKER_SIGAR_SPECIAL_RF = 11
    };

    enum WGS_SAM_FALGS_T {
        WGS_SAM_FLAG_READ_PAIRED,
        WGS_SAM_FLAG_READ_MAPPED_PROPER_PAIR,
        WGS_SAM_FLAG_READ_UNMAPPED,
        WGS_SAM_FLAG_MATE_UNMAPPED,
        WGS_SAM_FLAG_READ_REV_STRAND,
        WGS_SAM_FLAG_MATE_REV_STRAND,
        WGS_SAM_FLAG_FIRST_IN_PAIR,
        WGS_SAM_FLAG_SECOND_IN_PAIR,
        WGS_SAM_FLAG_N_PRIM_ALIGN,
        WGS_SAM_FLAG_READ_FAILES,
        WGS_SAM_FLAG_READ_IS_PCR,
        WGS_SAM_FLAG_READ_SUPPLE_ALIGN,
        WGS_SAM_FLAG_MAX
    };

// Read flag methods.
#define PAIRED(b)                     (((b)->core.flag & BAM_FPAIRED) != 0)
#define PROPERLY_PAIRED(b)            (((b)->core.flag & BAM_FPROPER_PAIR) != 0)
#define UNMAPPED(b)                   (((b)->core.flag & BAM_FUNMAP) != 0)
#define MATE_UNMAPPED(b)              (((b)->core.flag & BAM_FMUNMAP) != 0)
#define REVERSE_STRAND(b)             (((b)->core.flag & BAM_FREVERSE) != 0)
#define MATE_REVERSE_STRAND(b)        (((b)->core.flag & BAM_FMREVERSE) != 0)
#define FIRST_OF_PAIR(b)              (((b)->core.flag & BAM_FREAD1) != 0)
#define SECOND_OF_PAIR(b)             (((b)->core.flag & BAM_FREAD2) != 0)
#define SECONDARY_ALIGNMENT(b)        (((b)->core.flag & BAM_FSECONDARY) != 0)
#define FAILS_VENDOR_QUALITY_CHECK(b) (((b)->core.flag & BAM_FQCFAIL) != 0)
#define DUPLICATE(b)                  (((b)->core.flag & BAM_FDUP) != 0)
#define SUPPLEMENTARY_ALIGNMENT(b)    (((b)->core.flag & BAM_FSUPPLEMENTARY) != 0)

// read cigar methods.
#define consumes_read_bases(op) (bam_cigar_type(op) & 0x1)
#define consumes_ref_bases(op)  (bam_cigar_type(op) & 0x2)
#define consumes_all_bases(op)  (bam_cigar_type(op) & 0x3)

#define getBit(b, p)            (((b) >> (p)) & 1)
#define setBit(b, p)            ((b) |= (1 << (p)))
#define clearBit(b, p)          ((b) &= ~(1 << (p)))
#define assemble_min(A, B)      (((A) < (B)) ? (A) : (B))
#define assemble_max(A, B)      (((A) > (B)) ? (A) : (B))
#define likely(x)               __builtin_expect(!!(x), 1)
#define unlikely(x)             __builtin_expect(!!(x), 0)

#ifdef __cplusplus
}
#endif
#endif  // !__SHARED_HC_MAIN_H__