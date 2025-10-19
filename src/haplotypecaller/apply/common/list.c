/**
 * @file list.c
 * @author Caoce (caoce@genomics.cn)
 * @brief list，外部内存
 * @version 0.1
 * @date 2021-11-12
 *
 * @copyright Copyright (c) 2021 ROVACA SDK
 *
 */

#include "list.h"

#include "debug.h"

/**
 * @brief init
 *
 * @param list_block
 * @return p_hc_shared_lib_list_head
 */
p_hc_shared_lib_list_head hc_shared_list_init(int list_block)
{
    p_hc_shared_lib_list_head ret = malloc(sizeof(hc_shared_lib_list_head));

    if (__glibc_unlikely(!ret)) {
        statics_print_error("No mem");
        return NULL;
    }

    ret->nodes = 0;
    ret->pool = ring_mempool_auto_enlarge_init(sizeof(hc_shared_lib_list_item), list_block);
    if (__glibc_unlikely(!ret->pool)) {
        free(ret);
        statics_print_error("No mem");
        return NULL;
    }

    return ret;
}

/**
 * @brief insert
 *
 * @param head
 * @param data
 * @return int
 */
int hc_shared_list_insert(p_hc_shared_lib_list_head head, void* data)
{
    p_hc_shared_lib_list_item one_node;

    one_node = ring_mempool_auto_enlarge_malloc(head->pool);
    if (__glibc_unlikely(!one_node)) {
        return 0;
    }

    one_node->data = data;
    CDL_APPEND(head->nodes, one_node);

    return 1;
}

/**
 * @brief insert prev
 *
 * @param head
 * @param data
 * @return int
 */
int hc_shared_list_insert_prev(p_hc_shared_lib_list_head head, void* data)
{
    p_hc_shared_lib_list_item one_node;

    one_node = ring_mempool_auto_enlarge_malloc(head->pool);
    if (__glibc_unlikely(!one_node)) {
        return 0;
    }

    one_node->data = data;
    CDL_PREPEND(head->nodes, one_node);

    return 1;
}

/**
 * @brief search
 *
 * @param head
 * @param data
 * @return int
 */
int hc_shared_list_search(p_hc_shared_lib_list_head head, void* data)
{
    p_hc_shared_lib_list_item el;

    CDL_FOREACH(head->nodes, el)
    {
        if (el->data == data) {
            return 1;
        }
    }
    return 0;
}

/**
 * @brief delete
 *
 * @param head
 * @param data
 * @return int
 */
int hc_shared_list_delete(p_hc_shared_lib_list_head head, void* data)
{
    p_hc_shared_lib_list_item el, tmp1, tmp2;

    CDL_FOREACH_SAFE(head->nodes, el, tmp1, tmp2)
    {
        if (el->data == data) {
            CDL_DELETE(head->nodes, el);
            return 1;
        }
    }
    return 0;
}

/**
 * @brief reset
 *
 * @param head
 */
void hc_shared_list_reset(p_hc_shared_lib_list_head head)
{
    head->nodes = 0;
    ring_mempool_auto_enlarge_reset(head->pool);
}

/**
 * @brief term
 *
 * @param head
 */
void hc_shared_list_term(p_hc_shared_lib_list_head head)
{
    ring_mempool_auto_enlarge_term(head->pool);
    free(head);
}
