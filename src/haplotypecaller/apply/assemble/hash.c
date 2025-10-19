/**
 * @file
 * BQSR Hash
 *
 * 原始代码来自 https://github.com/skeeto/hashtab
 */

#include "hash.h"

#include "debug.h"
#include "utlist.h"

/* Default hashtable hash function. */
static inline uint32_t assemble_graph_hash_function(void* key, int key_size, uint32_t hashtab_size);

/**
 * @brief 新建hash table
 *
 * @param[in] size                  node总数
 * @param[in] equal_func            hash相等时执行函数
 * @param[in] new_node_memcpy       新建node时复制函数
 *
 * @return p_assemble_graph_hash_table hash table
 */
p_assemble_graph_hash_table assemble_graph_hash_init(int size, int (*equal_func)(void*, void*),
                                                     void (*new_node_memcpy)(void*, void*, int, p_assemble_graph_hash_table))
{
    p_assemble_graph_hash_table new_ht;

    new_ht = (p_assemble_graph_hash_table)malloc(sizeof(assemble_graph_hash_table));
    if (!new_ht) {
        statics_print_error("No hash memory");
        return NULL;
    }
    memset(new_ht, 0, sizeof(assemble_graph_hash_table));
    new_ht->recal_table_mem = ring_mempool_auto_enlarge_init(sizeof(assemble_graph_hash_node), ASSEMBLE_GRAPH_HASH_NODES);
    new_ht->arr = (p_assemble_graph_hash_node*)malloc(sizeof(void*) * size);
    if (!new_ht->arr || !new_ht->recal_table_mem) {
        statics_print_error("No hash memory");
        return NULL;
    }
    memset(new_ht->arr, 0, sizeof(void*) * size);

    new_ht->size = size;
    new_ht->equal_func = equal_func;
    new_ht->new_node_memcpy = new_node_memcpy;

    return new_ht;
}

/**
 * @brief 使用key keylen执行hash search
 *
 * @param[in] hashtable     hash table
 * @param[in] key           key
 * @param[in] keylen        key len
 *
 * @return void*    item
 */
void* assemble_graph_hash_search(p_assemble_graph_hash_table hashtable, void* key, int keylen)
{
    p_assemble_graph_hash_node last_node;
    uint32_t index = assemble_graph_hash_function(key, keylen, hashtable->size);

    if (__glibc_unlikely(hashtable->arr[index] == NULL)) {
        return NULL;
    }

    LL_FOREACH2(hashtable->arr[index], last_node, next_node)
    {
        if (memcmp(key, last_node->key, keylen) == 0) {
            return last_node->value;
        }
    }

    return NULL;
}

/**
 * @brief 获取hash key
 *
 * @param[in] key       value
 * @param[in] keylen    value len
 *
 * @return uint32_t     key
 */
inline uint32_t assemble_graph_hash_get_key(p_assemble_graph_hash_table hashtable, void* key, int keylen)
{
    return assemble_graph_hash_function(key, keylen, hashtable->size);
}

/**
 * @brief 使用key keylen index执行hash search
 *
 * @param[in] hashtable     hash table
 * @param[in] key           key
 * @param[in] keylen        key len
 *
 * @return p_assemble_graph_hash_node    node
 */
p_assemble_graph_hash_node assemble_graph_hash_search_with_key(p_assemble_graph_hash_table hashtable, void* key, int keylen, uint32_t index)
{
    p_assemble_graph_hash_node last_node;

    if (__glibc_unlikely(hashtable->arr[index] == NULL)) {
        return NULL;
    }

    LL_FOREACH2(hashtable->arr[index], last_node, next_node)
    {
        if (memcmp(key, last_node->key, keylen) == 0) {
            return last_node;
        }
    }

    return NULL;
}

/**
 * @brief 使用key keylen index执行hash search
 *
 * @param[in] hashtable     hash table
 * @param[in] key           key
 * @param[in] keylen        key len
 *
 * @return p_assemble_graph_hash_node    node
 */
p_assemble_graph_hash_node assemble_graph_hash_search_pt_with_key(p_assemble_graph_hash_table hashtable, void* value, uint32_t index)
{
    p_assemble_graph_hash_node last_node;

    if (__glibc_unlikely(hashtable->arr[index] == NULL)) {
        return NULL;
    }

    LL_FOREACH2(hashtable->arr[index], last_node, next_node)
    {
        if (last_node->value == value) {
            return last_node;
        }
    }

    return NULL;
}

/**
 * @brief 使用key keylen index执行hash insert
 *
 * @param[in] hashtable     hash table
 * @param[in] key           key
 * @param[in] keylen        key len
 * @param[in] value         value
 * @param[in] index         index
 *
 * @return void*            new node value
 */
void* assemble_graph_hash_insert_with_key(p_assemble_graph_hash_table hashtable, void* key, int keylen, void* value, uint32_t index)
{
    p_assemble_graph_hash_node last_node;
    p_assemble_graph_hash_node new_node;

    /* Search for an existing key. */
    LL_FOREACH2(hashtable->arr[index], last_node, next_node)
    {
        if (memcmp(key, last_node->key, keylen) == 0) {
            if (!hashtable->equal_func(last_node->value, value)) {
                break;
            }
            return last_node->value;
        }
    }

    /* create a new node */
    new_node = (p_assemble_graph_hash_node)ring_mempool_auto_enlarge_malloc(hashtable->recal_table_mem);
    if (new_node == NULL) {
        return NULL;
    }

    new_node->key = key;
    new_node->index = index;
    hashtable->new_node_memcpy(new_node, value, keylen, hashtable);
    new_node->index = index;
    new_node->next_node = NULL;
    new_node->flags = 0;

    /* Tack the new node on the end or right on the table. */
    LL_APPEND2(hashtable->arr[index], new_node, next_node);
    CDL_APPEND(hashtable->node_list, new_node);

    return new_node->value;
}

void* assemble_graph_hash_insert(p_assemble_graph_hash_table hashtable, void* key, int keylen, void* value)
{
    uint32_t index = assemble_graph_hash_function(key, keylen, hashtable->size);

    p_assemble_graph_hash_node last_node;
    p_assemble_graph_hash_node new_node;

    /* Search for an existing key. */
    LL_FOREACH2(hashtable->arr[index], last_node, next_node)
    {
        if (memcmp(key, last_node->key, keylen) == 0) {
            if (!hashtable->equal_func(last_node->value, value)) {
                break;
            }
            return last_node->value;
        }
    }

    /* create a new node */
    new_node = (p_assemble_graph_hash_node)ring_mempool_auto_enlarge_malloc(hashtable->recal_table_mem);
    if (new_node == NULL) {
        return NULL;
    }

    new_node->key = key;
    new_node->index = index;
    hashtable->new_node_memcpy(new_node, value, keylen, hashtable);
    new_node->index = index;
    new_node->next_node = NULL;

    /* Tack the new node on the end or right on the table. */
    LL_APPEND2(hashtable->arr[index], new_node, next_node);
    CDL_APPEND(hashtable->node_list, new_node);

    return new_node->value;
}

/* delete the given key from the hashtable */
void assemble_graph_hash_remove(p_assemble_graph_hash_table hashtable, void* key, int keylen)
{
    p_assemble_graph_hash_node next_node, tmp;
    int index = assemble_graph_hash_function(key, keylen, hashtable->size);
    next_node = hashtable->arr[index];

    LL_FOREACH_SAFE2(hashtable->arr[index], next_node, tmp, next_node)
    {
        if (memcmp(key, next_node->key, keylen) == 0) {
            LL_DELETE2(hashtable->arr[index], next_node, next_node);
            CDL_DELETE(hashtable->node_list, next_node);
        }
    }
}

/* delete the given key from the hashtable */
void assemble_graph_hash_remove_item(p_assemble_graph_hash_table hashtable, void* key, int keylen, void* data)
{
    p_assemble_graph_hash_node next_node, tmp;
    int index = assemble_graph_hash_function(key, keylen, hashtable->size);
    next_node = hashtable->arr[index];

    LL_FOREACH_SAFE2(hashtable->arr[index], next_node, tmp, next_node)
    {
        if (next_node->value == data) {
            LL_DELETE2(hashtable->arr[index], next_node, next_node);
            CDL_DELETE(hashtable->node_list, next_node);
        }
    }
}

void assemble_graph_hash_remove_item_with_key(p_assemble_graph_hash_table hashtable, void* data, uint32_t index)
{
    p_assemble_graph_hash_node next_node, tmp;
    next_node = hashtable->arr[index];

    LL_FOREACH_SAFE2(hashtable->arr[index], next_node, tmp, next_node)
    {
        if (next_node->value == data) {
            LL_DELETE2(hashtable->arr[index], next_node, next_node);
            CDL_DELETE(hashtable->node_list, next_node);
        }
    }
}

/**
 * @brief 释放hash table上所有item
 *
 * @param hashtable hashtable
 *
 * @note 不处理运行时释放
 */
void assemble_graph_hash_reset(p_assemble_graph_hash_table hashtable)
{
    p_assemble_graph_hash_node el, tmp1, tmp2;

    CDL_FOREACH_SAFE(hashtable->node_list, el, tmp1, tmp2) { hashtable->arr[el->index] = NULL; }
    hashtable->node_list = NULL;
    ring_mempool_auto_enlarge_reset(hashtable->recal_table_mem);
}

/* free all resources used by the hashtable */
void assemble_graph_hash_destroy(p_assemble_graph_hash_table hashtable)
{
    ring_mempool_auto_enlarge_term(hashtable->recal_table_mem);
    free(hashtable->arr);
    free(hashtable);
}

#include "uthash.h"
/**
 *  Compute the hashcode for a KMer.
 *  Equivalent to <code>new String(bases, start, length).hashCode()</code>
 */
static uint32_t assemble_graph_hash_function(void* key, int length, uint32_t hashtab_size)
{
    if (length == 0) {
        return 0;
    }

    uint32_t res = 0;
    HASH_SFH(key, length, res);
    return res & (hashtab_size - 1);
}
