/**
 * @file
 * BQSR PR MemPool
 */

#include "ring_mem_pool.h"

#include <sys/mman.h>

#define PROTECTION (PROT_READ | PROT_WRITE)

#ifndef MAP_HUGETLB
#define MAP_HUGETLB 0x40000 /* arch specific */
#endif

/* Only ia64 requires this */
#ifdef __ia64__
#define ADDR  (void *)(0x8000000000000000UL)
#define FLAGS (MAP_PRIVATE | MAP_ANONYMOUS | MAP_HUGETLB | MAP_FIXED)
#else
#define ADDR  (void *)(0x0UL)
#define FLAGS (MAP_PRIVATE | MAP_ANONYMOUS | MAP_HUGETLB)
#endif

/**
 * Initialize a memory pool
 *
 * @param pool the pool
 * @param size the size of a memory block
 * @param count the count of memory blocks
 * @return 1 on success, 0 otherwise
 */
int ring_mempool_init(p_ring_mempool pool, unsigned int size, unsigned int count)
{
    unsigned long long wantsize;

    memset(pool, 0, sizeof(ring_mempool));

    if (!size || !count) {
        return 0;
    }
    wantsize = (unsigned long long)size * (unsigned long long)count;
    char *memory = (char *)RING_MEMPOOL_MALLOC(wantsize);
    wantsize = sizeof(void *) * count;
    pool->blocks = (void **)RING_MEMPOOL_MALLOC(wantsize);

    if (!memory || !pool->blocks) {
        RING_MEMPOOL_FREE(memory);
        RING_MEMPOOL_FREE(pool->blocks);
        return 0;
    }

    pool->memory = memory;
    pool->count = count;
    pool->size = size;
    pool->mask = count - 1;
    pool->head_index = count - 1;
    pool->tail_index = 0;

    if (pthread_cond_init(&pool->cnd_mem, NULL) != 0) {
        return 0;
    }
    if (pthread_mutex_init(&pool->mutex, NULL) != 0) {
        pthread_cond_destroy(&(pool->cnd_mem));
        return 0;
    }

    memory -= size;
    while (count--) {
        pool->blocks[count] = memory + size;
        memory += size;
    }

    return 1;
}

/**
 * Finalize a memory pool
 *
 * @param pool the pool
 */
void ring_mempool_term(p_ring_mempool pool)
{
    RING_MEMPOOL_FREE(pool->blocks);
    RING_MEMPOOL_FREE(pool->memory);
    pthread_cond_destroy(&(pool->cnd_mem));
    pthread_mutex_destroy(&(pool->mutex));
    memset(pool, 0, sizeof(ring_mempool));
}

/**
 * Get a memory block from the pool
 *
 * @param pool the pool
 * @return a pointer to memory block on success, NULL on failure
 */
void *ring_mempool_alloc(p_ring_mempool pool)
{
    volatile void *res = NULL;

    // 尝试等待获取
    if (pool->tail_index == (pool->head_index & pool->mask)) {
        return NULL;
    }

    res = pool->blocks[pool->tail_index];
    // pool->blocks[pool->tail_index] = NULL;
    pool->tail_index = ((pool->tail_index + 1) & pool->mask);
    return (void *)res;
}

/**
 * Free a memory block
 *
 * @param pool the pool
 * @param block the block
 */
void ring_mempool_free(p_ring_mempool pool, void *block)
{
    uint32_t free_block;

    free_block = (__sync_fetch_and_add(&pool->head_index, 1) & pool->mask);
    pool->blocks[free_block] = block;
}

/**
 * Restore a memory block
 *
 * @param pool the pool
 */
void ring_mempool_restore(p_ring_mempool pool)
{
    char *memory = pool->memory;
    unsigned int count = pool->count;

    pool->head_index = pool->count - 1;
    pool->tail_index = 0;

    memory -= pool->size;
    while (count--) {
        pool->blocks[count] = memory + pool->size;
        memory += pool->size;
    }
}

int ring_mempool_init_with_huge(p_ring_mempool pool, unsigned int size, unsigned int count)
{
    unsigned long long wantsize;

    memset(pool, 0, sizeof(ring_mempool));

    if (!size || !count) {
        return 0;
    }
    wantsize = (unsigned long long)size * (unsigned long long)count;
    char *memory = mmap(ADDR, wantsize, PROTECTION, FLAGS, 0, 0);
    wantsize = sizeof(void *) * count;
    pool->blocks = mmap(ADDR, wantsize, PROTECTION, FLAGS, 0, 0);

    if (memory == MAP_FAILED || !pool->blocks) {
        munmap(memory, (unsigned long long)size * (unsigned long long)count);
        if (pool->blocks) {
            munmap(pool->blocks, wantsize);
        }
        return 0;
    }

    pool->memory = memory;
    pool->count = count;
    pool->size = size;
    pool->mask = count - 1;
    pool->head_index = count - 1;
    pool->tail_index = 0;

    if (pthread_cond_init(&pool->cnd_mem, NULL) != 0) {
        return 0;
    }
    if (pthread_mutex_init(&pool->mutex, NULL) != 0) {
        pthread_cond_destroy(&(pool->cnd_mem));
        return 0;
    }

    memory -= size;
    while (count--) {
        pool->blocks[count] = memory + size;
        memory += size;
    }

    return 1;
}
void ring_mempool_term_with_huge(p_ring_mempool pool)
{
    munmap(pool->blocks, pool->count * sizeof(void *));
    munmap(pool->memory, pool->count * pool->size);
    pthread_cond_destroy(&(pool->cnd_mem));
    pthread_mutex_destroy(&(pool->mutex));
    memset(pool, 0, sizeof(ring_mempool));
}
