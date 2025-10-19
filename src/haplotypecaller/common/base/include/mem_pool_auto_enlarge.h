/**
 * @file mem_pool_auto_alloc.h
 * @author Cao Ce (caoce@genomics.cn)
 * @brief 自动扩充内存的内存池，不可重入，不可同时申请、释放
 * @version 0.1
 * @date 2021-05-24
 *
 * @copyright Copyright (c) 2021
 *
 */

#ifndef __LIB_MEM_POOL_AUTO_ENLARGE__
#define __LIB_MEM_POOL_AUTO_ENLARGE__

#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "utlist.h"

#ifndef MEM_POOL_AUTO_ENLARGE_MALLOC
#define MEM_POOL_AUTO_ENLARGE_MALLOC(size) malloc(size)
#endif
#ifndef MEM_POOL_AUTO_ENLARGE_FREE
#define MEM_POOL_AUTO_ENLARGE_FREE(pt) free(pt)
#endif

#define MEM_POOL_AUTO_ENLARGE_ENLARGE (2)

typedef struct mem_pool_auto_enlarge_mem_t
{
    void* mem;
    void* next_data;
    unsigned int count;

    struct mem_pool_auto_enlarge_mem_t* prev;
    struct mem_pool_auto_enlarge_mem_t* next;
} mem_pool_auto_enlarge_mem, *p_mem_pool_auto_enlarge_mem;

typedef struct auto_enlarge_mempool_t
{
    p_mem_pool_auto_enlarge_mem memory;
    p_mem_pool_auto_enlarge_mem run_pool;

    // 初始化数据
    unsigned int size;
    unsigned int count;
} auto_enlarge_mempool, *p_auto_enlarge_mempool;

p_auto_enlarge_mempool ring_mempool_auto_enlarge_init(unsigned int size, unsigned int count);
void ring_mempool_auto_enlarge_term(p_auto_enlarge_mempool pool);
void ring_mempool_enalrge_pool(p_auto_enlarge_mempool pool);
void* ring_mempool_auto_enlarge_malloc(p_auto_enlarge_mempool pool);
void ring_mempool_auto_enlarge_reset(p_auto_enlarge_mempool pool);
void ring_mempool_auto_enlarge_reset_max(p_auto_enlarge_mempool pool, uint32_t max);

typedef struct mem_pool_fast_main_storage_t mem_pool_fast_main_storage, *p_mem_pool_fast_main_storage;
p_mem_pool_fast_main_storage mem_pool_fast_init(uint32_t size, uint32_t items);
void* mem_pool_fast_alloc(p_mem_pool_fast_main_storage in);
void mem_pool_fast_free(p_mem_pool_fast_main_storage pool, void* block);
void mem_pool_fast_reset(p_mem_pool_fast_main_storage pool);
void mem_pool_fast_term(p_mem_pool_fast_main_storage pool);

#endif  // !__LIB_MEM_POOL_AUTO_ENLARGE__