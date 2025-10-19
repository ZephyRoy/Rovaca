#ifndef __SHARED_HC_ASSEMBLE_SIMPLE_THREAD_MEM_H__
#define __SHARED_HC_ASSEMBLE_SIMPLE_THREAD_MEM_H__

#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "utlist.h"

typedef struct hc_assemble_thread_mem_item_t
{
    unsigned int used;

    struct hc_assemble_thread_mem_item_t* prev;
    struct hc_assemble_thread_mem_item_t* next;

    uint8_t mem[];
} hc_assemble_thread_mem_item, *p_hc_assemble_thread_mem_item;

typedef struct hc_assemble_thread_mem_t
{
    p_hc_assemble_thread_mem_item memory;
    p_hc_assemble_thread_mem_item run_pool;

    // 初始化数据
    unsigned int size;
} hc_assemble_thread_mem, *p_hc_assemble_thread_mem;

p_hc_assemble_thread_mem hc_assemble_thread_mem_init(unsigned int size);
void hc_assemble_thread_mem_term(p_hc_assemble_thread_mem pool);
void* hc_assemble_thread_mem_malloc(p_hc_assemble_thread_mem pool, uint32_t wantsize);
void hc_assemble_thread_mem_reset(p_hc_assemble_thread_mem pool);

#endif  // !__SHARED_HC_ASSEMBLE_SIMPLE_THREAD_MEM_H__