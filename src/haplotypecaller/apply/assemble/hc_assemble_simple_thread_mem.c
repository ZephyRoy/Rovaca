/**
 * @file hc_assemble_simple_thread_mem.c
 * @author Cao Ce (caoce@genomics.cn)
 * @brief 简易的线程内存池
 * @version 0.1
 * @date 2022-02-09
 *
 * @copyright Copyright (c) 2022 ROVACA SDK
 *
 */

#include "hc_assemble_simple_thread_mem.h"

#include "debug.h"
#include "hc_assemble.h"
#include "mem_pool_auto_enlarge.h"

static void hc_assemble_thread_mem_enalrge_pool(p_hc_assemble_thread_mem pool);

/**
 * @brief 初始化一个内存池
 *
 * @param size          单个item大小
 * @param count         总item数量
 *
 * @return p_hc_assemble_thread_mem       成功
 *         NULL                         失败
 */
p_hc_assemble_thread_mem hc_assemble_thread_mem_init(unsigned int size)
{
    p_hc_assemble_thread_mem_item mem = NULL;
    p_hc_assemble_thread_mem pool = NULL;

    if (!size) {
        return NULL;
    }
    pool = (p_hc_assemble_thread_mem)malloc(sizeof(hc_assemble_thread_mem));
    mem = (p_hc_assemble_thread_mem_item)malloc(sizeof(hc_assemble_thread_mem_item) + size);

    if (!mem || !pool) {
        free(pool);
        free(mem);
        return NULL;
    }
    memset(pool, 0, sizeof(auto_enlarge_mempool));

    mem->used = 0;
    pool->run_pool = mem;
    pool->size = size;

    CDL_APPEND(pool->memory, mem);

    return pool;
}

/**
 * @brief 取消一个内存池
 *
 * @param pool          内存池
 */
void hc_assemble_thread_mem_term(p_hc_assemble_thread_mem pool)
{
    p_hc_assemble_thread_mem_item cache, tmp1, tmp2;

    CDL_FOREACH_SAFE(pool->memory, cache, tmp1, tmp2)
    {
        CDL_DELETE(pool->memory, cache);
        free(cache);
    }

    free(pool);
}

/**
 * @brief 自动扩增内存池
 *
 * @param pool          内存池
 */
static void hc_assemble_thread_mem_enalrge_pool(p_hc_assemble_thread_mem pool)
{
    p_hc_assemble_thread_mem_item mem = NULL;

    if (pool->run_pool->next) {
        pool->run_pool = pool->run_pool->next;
        return;
    }

    mem = (p_hc_assemble_thread_mem_item)malloc(sizeof(hc_assemble_thread_mem_item) + pool->size);

    if (!mem) {
        statics_print_error("No Memory");
        exit(2);
    }

    CDL_APPEND(pool->memory, mem);
    mem->used = 0;
    pool->run_pool = mem;
}

/**
 * @brief 从内存池Alloc一个item
 *
 * @param pool          内存池
 *
 * @return void*        item
 */
void* hc_assemble_thread_mem_malloc(p_hc_assemble_thread_mem pool, uint32_t wantsize)
{
    void* res = NULL;

    while (__glibc_unlikely(pool->run_pool->used + wantsize >= pool->size)) {
        hc_assemble_thread_mem_enalrge_pool(pool);
    }

    res = pool->run_pool->mem + pool->run_pool->used;
    pool->run_pool->used += wantsize;
    return res;
}

/**
 * @brief 重置内存池
 *
 * @param pool          内存池
 */
void hc_assemble_thread_mem_reset(p_hc_assemble_thread_mem pool)
{
    p_hc_assemble_thread_mem_item cache, tmp1, tmp2;

    CDL_FOREACH_SAFE(pool->memory, cache, tmp1, tmp2)
    {
        cache->used = 0;
        if (cache == pool->memory) {
            continue;
        }
        CDL_DELETE(pool->memory, cache);
        free(cache);
    }
    pool->run_pool = pool->memory;
}
