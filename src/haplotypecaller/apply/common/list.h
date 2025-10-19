#ifndef __SHARED_LIB_LIST_H__
#define __SHARED_LIB_LIST_H__

#include "mem_pool_auto_enlarge.h"
#include "utlist.h"

typedef struct hc_shared_lib_list_item_t hc_shared_lib_list_item, *p_hc_shared_lib_list_item;

typedef struct hc_shared_lib_list_item_t
{
    void* data;

    p_hc_shared_lib_list_item prev;
    p_hc_shared_lib_list_item next;
} hc_shared_lib_list_item, *p_hc_shared_lib_list_item;

typedef struct hc_shared_lib_list_head_t
{
    p_hc_shared_lib_list_item nodes;

    p_auto_enlarge_mempool pool;
} hc_shared_lib_list_head, *p_hc_shared_lib_list_head;

p_hc_shared_lib_list_head hc_shared_list_init(int list_block);
int hc_shared_list_insert(p_hc_shared_lib_list_head head, void* data);
int hc_shared_list_insert_prev(p_hc_shared_lib_list_head head, void* data);
int hc_shared_list_delete(p_hc_shared_lib_list_head head, void* data);
void hc_shared_list_reset(p_hc_shared_lib_list_head head);
void hc_shared_list_term(p_hc_shared_lib_list_head head);
int hc_shared_list_search(p_hc_shared_lib_list_head head, void* data);

#endif