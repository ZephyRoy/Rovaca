#ifndef __ASSEMBLE_GRAPH_HASH_H__
#define __ASSEMBLE_GRAPH_HASH_H__
#include <stdint.h>

#include "mem_pool_auto_enlarge.h"

#define ASSEMBLE_GRAPH_HASH_NODES  (32768)
#define HC_ASSEMBLE_GRAPH_HASH_MUL (31)

typedef struct assemble_graph_hash_node_t
{
    void* key;      /* key for the node */
    int keylen;     /* length of the key */
    void* value;    /* value for this node */
    uint32_t index; /* index */
    uint32_t flags; /* flags */

    struct assemble_graph_hash_node_t* next_node; /* next node (open hashtable) */

    struct assemble_graph_hash_node_t* prev; /* hash iter CDL */
    struct assemble_graph_hash_node_t* next; /* hash iter CDL */
} assemble_graph_hash_node, *p_assemble_graph_hash_node;

typedef struct assemble_graph_hash_table_t
{
    p_assemble_graph_hash_node* arr;
    p_assemble_graph_hash_node node_list;

    p_auto_enlarge_mempool recal_table_mem;  // recal table线程内存
    int size;                                /* size of the hash */

    int (*equal_func)(void*, void*);
    void (*new_node_memcpy)(void*, void*, int, struct assemble_graph_hash_table_t*);
} assemble_graph_hash_table, *p_assemble_graph_hash_table;

/* Initialize a new hashtable (set bookingkeeping data) and return a
 * pointer to the hashtable. A hash function may be provided. If no
 * function pointer is given (a NULL pointer), then the built in hash
 * function is used. A NULL pointer returned if the creation of the
 * hashtable failed due to a failed malloc(). */
p_assemble_graph_hash_table assemble_graph_hash_init(int size, int (*equal_func)(void*, void*),
                                                     void (*new_node_memcpy)(void*, void*, int, p_assemble_graph_hash_table));

/* Fetch a value from table matching the key. Returns a pointer to
 * the value matching the given key. */
void* assemble_graph_hash_search(p_assemble_graph_hash_table hashtable, void* key, int keylen);

inline uint32_t assemble_graph_hash_get_key(p_assemble_graph_hash_table hashtable, void* key, int keylen);
p_assemble_graph_hash_node assemble_graph_hash_search_with_key(p_assemble_graph_hash_table hashtable, void* key, int keylen,
                                                               uint32_t index);
p_assemble_graph_hash_node assemble_graph_hash_search_pt_with_key(p_assemble_graph_hash_table hashtable, void* value, uint32_t index);
void* assemble_graph_hash_insert_with_key(p_assemble_graph_hash_table hashtable, void* key, int keylen, void* value, uint32_t index);

/* Put a value into the table with the given key. Returns NULL if
 * malloc() fails to allocate memory for the new node. */
void* assemble_graph_hash_insert(p_assemble_graph_hash_table hashtable, void* key, int keylen, void* value);

/* Delete the given key and value pair from the hashtable. If the key
 * does not exist, no error is given. */
void assemble_graph_hash_remove(p_assemble_graph_hash_table hashtable, void* key, int keylen);
void assemble_graph_hash_remove_item(p_assemble_graph_hash_table hashtable, void* key, int keylen, void* data);
void assemble_graph_hash_remove_item_with_key(p_assemble_graph_hash_table hashtable, void* data, uint32_t index);

/* Free all resources used by the hashtable. */
void assemble_graph_hash_destroy(p_assemble_graph_hash_table hashtable);
void assemble_graph_hash_reset(p_assemble_graph_hash_table hashtable);

#endif  // !__ASSEMBLE_GRAPH_HASH_H__