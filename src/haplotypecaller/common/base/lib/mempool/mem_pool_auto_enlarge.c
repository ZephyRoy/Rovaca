/**
 * @file mem_pool_auto_enlarge.c
 * @author Cao Ce (caoce@genomics.cn)
 * @brief 自动扩充内存的内存池，不可重入，不可同时申请、释放
 * @version 0.1
 * @date 2021-11-18
 *
 * @copyright Copyright (c) 2021 ROVACA SDK
 *
 */

#include "mem_pool_auto_enlarge.h"

/**
 * @brief 初始化一个内存池
 *
 * @param size          单个item大小
 * @param count         总item数量
 *
 * @return p_auto_enlarge_mempool       成功
 *         NULL                         失败
 */
p_auto_enlarge_mempool ring_mempool_auto_enlarge_init(unsigned int size, unsigned int count)
{
    unsigned long long wantsize;
    p_mem_pool_auto_enlarge_mem mem = NULL;
    p_auto_enlarge_mempool pool = NULL;
    char *memory;

    if (!size || !count) {
        return NULL;
    }
    wantsize = (unsigned long long)size * (unsigned long long)count;
    pool = (p_auto_enlarge_mempool)MEM_POOL_AUTO_ENLARGE_MALLOC(sizeof(auto_enlarge_mempool));
    memory = (char *)MEM_POOL_AUTO_ENLARGE_MALLOC(wantsize);
    mem = (p_mem_pool_auto_enlarge_mem)MEM_POOL_AUTO_ENLARGE_MALLOC(sizeof(mem_pool_auto_enlarge_mem));

    if (!memory || !mem || !pool) {
        MEM_POOL_AUTO_ENLARGE_FREE(pool);
        MEM_POOL_AUTO_ENLARGE_FREE(mem);
        MEM_POOL_AUTO_ENLARGE_FREE(memory);
        return NULL;
    }
    memset(pool, 0, sizeof(auto_enlarge_mempool));

    mem->mem = memory;
    mem->count = count;
    mem->next_data = memory;
    pool->run_pool = mem;
    pool->count = count;
    pool->size = size;

    DL_APPEND(pool->memory, mem);

    return pool;
}

/**
 * @brief 取消一个内存池
 *
 * @param pool          内存池
 */
void ring_mempool_auto_enlarge_term(p_auto_enlarge_mempool pool)
{
    p_mem_pool_auto_enlarge_mem cache, tmp;

    DL_FOREACH_SAFE(pool->memory, cache, tmp)
    {
        DL_DELETE(pool->memory, cache);
        MEM_POOL_AUTO_ENLARGE_FREE(cache->mem);
        MEM_POOL_AUTO_ENLARGE_FREE(cache);
    }

    MEM_POOL_AUTO_ENLARGE_FREE(pool);
}

/**
 * @brief 自动扩增内存池
 *
 * @param pool          内存池
 */
void ring_mempool_enalrge_pool(p_auto_enlarge_mempool pool)
{
    unsigned long long wantsize;
    p_mem_pool_auto_enlarge_mem mem = NULL;
    char *memory = NULL;

    if (pool->run_pool->next && pool->run_pool->next->count) {
        pool->run_pool = pool->run_pool->next;
        return;
    }

    wantsize = ((unsigned long long)pool->size * (unsigned long long)pool->count);
    mem = (p_mem_pool_auto_enlarge_mem)MEM_POOL_AUTO_ENLARGE_MALLOC(sizeof(mem_pool_auto_enlarge_mem));
    memory = (char *)MEM_POOL_AUTO_ENLARGE_MALLOC(wantsize);

    if (!memory || !mem) {
        MEM_POOL_AUTO_ENLARGE_FREE(mem);
        MEM_POOL_AUTO_ENLARGE_FREE(memory);
        fprintf(stderr, "[Lib Auto Enlarge Pool] No Memory");
        exit(0);
    }

    DL_APPEND(pool->memory, mem);
    mem->mem = memory;
    mem->count = pool->count;
    mem->next_data = memory;
    pool->run_pool = mem;
}

/**
 * @brief 从内存池Alloc一个item
 *
 * @param pool          内存池
 *
 * @return void*        item
 */
void *ring_mempool_auto_enlarge_malloc(p_auto_enlarge_mempool pool)
{
    void *res = NULL;

    if (pool->run_pool->count == 0) {
        ring_mempool_enalrge_pool(pool);
    }

    res = pool->run_pool->next_data;
    pool->run_pool->count -= 1;
    pool->run_pool->next_data = (uint8_t *)(pool->run_pool->next_data) + pool->size;
    return res;
}

/**
 * @brief 重置内存池
 *
 * @param pool          内存池
 */
void ring_mempool_auto_enlarge_reset(p_auto_enlarge_mempool pool)
{
    p_mem_pool_auto_enlarge_mem cache, tmp;

    DL_FOREACH_SAFE(pool->memory, cache, tmp)
    {
        if (cache == pool->memory) {
            continue;
        }
        DL_DELETE(pool->memory, cache);
        MEM_POOL_AUTO_ENLARGE_FREE(cache->mem);
        MEM_POOL_AUTO_ENLARGE_FREE(cache);
    }
    pool->run_pool = pool->memory;
    cache = pool->run_pool;
    cache->count = pool->count;
    cache->next_data = cache->mem;
}

/**
 * @brief 重置内存池,保留最大max数量的元素,防止内存过大
 *
 * @param pool          内存池
 */
void ring_mempool_auto_enlarge_reset_max(p_auto_enlarge_mempool pool, uint32_t max)
{
    p_mem_pool_auto_enlarge_mem cache, tmp;
    uint32_t now_index = 0;
    DL_FOREACH_SAFE(pool->memory, cache, tmp)
    {
        if (now_index >= max) {
            DL_DELETE(pool->memory, cache);
            MEM_POOL_AUTO_ENLARGE_FREE(cache->mem);
            MEM_POOL_AUTO_ENLARGE_FREE(cache);
            now_index++;
            continue;
        }
        now_index++;
        cache->count = pool->count;
        cache->next_data = cache->mem;
    }

    pool->run_pool = pool->memory;
}
