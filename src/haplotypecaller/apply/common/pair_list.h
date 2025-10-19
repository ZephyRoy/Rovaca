#ifndef __SHARED_LIB_PAIR_LIST_H__
#define __SHARED_LIB_PAIR_LIST_H__

#include "mem_pool_auto_enlarge.h"
#include "utlist.h"

typedef struct hc_shared_lib_pair_list_item_t hc_shared_lib_pair_list_item, *p_hc_shared_lib_pair_list_item;

typedef struct hc_shared_lib_pair_list_item_t
{
    void* first_data;
    union
    {
        void* second_data;
        struct
        {
            uint32_t second_data_uint32_1;
            uint32_t second_data_uint32_2;
        } two_uint_st;
    };

    // void*                               second_data;

    p_hc_shared_lib_pair_list_item prev;
    p_hc_shared_lib_pair_list_item next;
} hc_shared_lib_pair_list_item, *p_hc_shared_lib_pair_list_item;

typedef struct hc_shared_lib_pair_list_head_t
{
    p_hc_shared_lib_pair_list_item nodes;

    p_auto_enlarge_mempool pool;
} hc_shared_lib_pair_list_head, *p_hc_shared_lib_pair_list_head;

p_hc_shared_lib_pair_list_head hc_shared_pair_list_init(int list_block);
int hc_shared_pair_list_insert(p_hc_shared_lib_pair_list_head head, void* first_data, void* second_data);
int hc_shared_pair_list_delete(p_hc_shared_lib_pair_list_head head, void* first_data, void* second_data);
int hc_shared_pair_list_search(p_hc_shared_lib_pair_list_head head, void* first_data, void* second_data);
void hc_shared_pair_list_reset(p_hc_shared_lib_pair_list_head head);
void hc_shared_pair_list_term(p_hc_shared_lib_pair_list_head head);

#endif