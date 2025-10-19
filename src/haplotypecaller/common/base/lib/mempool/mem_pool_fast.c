/**
 * @file mem_pool_fast.c
 * @author Cao Ce (caoce@genomics.cn)
 * @brief 快速内存池，带释放
 * @version 0.1
 * @date 2022-07-13
 *
 * @copyright Copyright (c) 2022 ROVACA SDK
 *
 */

#include <stdint.h>
#include <stdlib.h>
#include <utlist.h>

typedef struct mem_pool_fast_one_mem_t
{
    struct mem_pool_fast_one_mem_t* next;

    uint8_t data[];
} mem_pool_fast_one_mem, *p_mem_pool_fast_one_mem;
typedef struct mem_pool_fast_one_mem_t mem_pool_fast_one_page, *p_mem_pool_fast_one_page;

typedef struct mem_pool_fast_main_storage_t
{
    uint32_t item_size;
    uint32_t item_sum;
    uint32_t extended;

    p_mem_pool_fast_one_mem mem_list;
    p_mem_pool_fast_one_page all_page;
} mem_pool_fast_main_storage, *p_mem_pool_fast_main_storage;

/**
 * @brief 内存池新增一个页
 *
 * @param[in] in 内存池
 *
 * @return p_mem_pool_fast_main_storage
 */
static inline p_mem_pool_fast_main_storage mem_pool_fast_add_one_page(p_mem_pool_fast_main_storage in)
{
    p_mem_pool_fast_one_page one_page;
    p_mem_pool_fast_one_mem mem;
    uint32_t page_size;
    register uint32_t i;

    page_size = (in->item_size + sizeof(mem_pool_fast_one_mem)) * in->item_sum + sizeof(mem_pool_fast_one_page);
    one_page = malloc(page_size);
    if (__glibc_unlikely(!one_page)) {
        return NULL;
    }
    LL_PREPEND(in->all_page, one_page);
    for (i = 0; i < in->item_sum; i++) {
        mem = (p_mem_pool_fast_one_mem)(&(one_page->data[i * (in->item_size + sizeof(mem_pool_fast_one_mem))]));
        LL_PREPEND(in->mem_list, mem);
    }
    in->extended = 1;
    return in;
}

/**
 * @brief 内存池初始化
 *
 * @param[in] size          item 大小
 * @param[in] items         item 总量
 * @param[in] alloc_pool    外部 pool
 *
 * @return p_mem_pool_fast_main_storage
 */
p_mem_pool_fast_main_storage mem_pool_fast_init(uint32_t size, uint32_t items)
{
    p_mem_pool_fast_main_storage ret;

    ret = malloc(sizeof(mem_pool_fast_main_storage));
    if (__glibc_unlikely(!ret)) {
        return NULL;
    }

    ret->item_size = size;
    ret->item_sum = items;
    ret->mem_list = NULL;
    ret->all_page = NULL;

    if (__glibc_unlikely(!mem_pool_fast_add_one_page(ret))) {
        free(ret);
        return NULL;
    }

    ret->extended = 0;
    return ret;
}

/**
 * @brief 内存池alloc
 *
 * @param[in] in        内存池
 *
 * @return void*
 */
void* mem_pool_fast_alloc(p_mem_pool_fast_main_storage in)
{
    p_mem_pool_fast_one_mem mem;

    if (__glibc_unlikely(!in->mem_list && !mem_pool_fast_add_one_page(in))) {
        return NULL;
    }
    mem = in->mem_list;
    LL_DELETE(in->mem_list, in->mem_list);
    return mem;
}

/**
 * @brief 内存池free
 *
 * @param[in] pool      内存池
 * @param[in] block     free
 */
void mem_pool_fast_free(p_mem_pool_fast_main_storage pool, void* block)
{
    p_mem_pool_fast_one_mem mem = block;

    LL_PREPEND(pool->mem_list, mem);
}

/**
 * @brief 内存池重置
 *
 * @param[in] pool      内存池
 */
void mem_pool_fast_reset(p_mem_pool_fast_main_storage pool)
{
    p_mem_pool_fast_one_page one_page, tmp;
    p_mem_pool_fast_one_mem mem;
    register uint32_t i;

    if (__glibc_likely(!pool->extended)) {
        return;
    }

    pool->mem_list = NULL;
    LL_FOREACH_SAFE(pool->all_page, one_page, tmp)
    {
        if (one_page->next) {
            LL_DELETE(pool->all_page, one_page);
            free(one_page);
            continue;
        }

        for (i = 0; i < pool->item_sum; i++) {
            mem = (p_mem_pool_fast_one_mem)(one_page->data + i * (pool->item_size + sizeof(mem_pool_fast_one_mem)));
            LL_PREPEND(pool->mem_list, mem);
        }
    }
    pool->extended = 0;
}

unsigned int mem_pool_fast_count(p_mem_pool_fast_main_storage pool)
{
    unsigned int ret = 0;
    p_mem_pool_fast_one_mem mem;

    LL_COUNT(pool->mem_list, mem, ret);
    return ret;
}

/**
 * @brief 内存池释放
 *
 * @param[in] pool      内存池
 */
void mem_pool_fast_term(p_mem_pool_fast_main_storage pool)
{
    p_mem_pool_fast_one_page one_page, tmp;

    LL_FOREACH_SAFE(pool->all_page, one_page, tmp)
    {
        LL_DELETE(pool->all_page, one_page);
        free(one_page);
    }
    free(pool);
}
