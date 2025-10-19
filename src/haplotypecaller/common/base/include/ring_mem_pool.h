#ifndef __LIB_RING_MEM_POOL__
#define __LIB_RING_MEM_POOL__

#include <pthread.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>

#ifndef RING_MEMPOOL_MALLOC
#define RING_MEMPOOL_MALLOC(size) malloc(size)
#endif
#ifndef RING_MEMPOOL_FREE
#define RING_MEMPOOL_FREE(pt) free(pt)
#endif

typedef struct ring_mempool_t
{
    volatile uint32_t tail_index;
    volatile uint32_t head_index;
    void** blocks;
    uint32_t mask;
    char* memory;
    pthread_cond_t cnd_mem;
    pthread_mutex_t mutex;
    unsigned int count;
    uint32_t size;
} ring_mempool, *p_ring_mempool;

int ring_mempool_init(p_ring_mempool pool, unsigned int size, unsigned int count);
void ring_mempool_term(p_ring_mempool pool);
void* ring_mempool_alloc(p_ring_mempool pool);
void ring_mempool_free(p_ring_mempool pool, void* block);
void ring_mempool_restore(p_ring_mempool pool);

int ring_mempool_init_with_huge(p_ring_mempool pool, unsigned int size, unsigned int count);
void ring_mempool_term_with_huge(p_ring_mempool pool);

#endif  // !__LIB_RING_MEM_POOL__