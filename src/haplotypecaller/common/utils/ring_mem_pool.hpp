#ifndef _RING_MEM_POOL_H_
#define _RING_MEM_POOL_H_

#include <cstdlib>
#include <cstring>

#include "htslib/sam.h"
static uint32_t get_mask(uint32_t count)
{
    uint32_t mask = 1;
    while (mask < count) {
        mask = (mask << 1) | 1;
    }
    return mask;
}
template <typename T>
class RingMemPool
{
public:
    RingMemPool()
        : memory_(nullptr)
        , tail_index_(0)
        , head_index_(0)
        , blocks_(nullptr)
        , mask_(0)
        , size_()
        , count_(0)
    {}

    ~RingMemPool() {}

    int init(size_t size, uint32_t count)
    {
        if (count == 0) {
            return 0;
        }

        unsigned long long wantsize = (unsigned long long)size * (unsigned long long)count;
        uint8_t* memory = new uint8_t[wantsize];

        wantsize = sizeof(void*) * count;
        blocks_ = new void*[wantsize];

        memory_ = memory;
        count_ = count;
        size_ = size;
        mask_ = get_mask(count);
        head_index_ = count - 1;
        tail_index_ = 0;

        memory -= size;
        while (count--) {
            blocks_[count] = memory + size;
            memory += size;
        }
        return 1;
    }

    void term()
    {
        delete[] blocks_;
        delete[] memory_;
    }

    T* alloc()
    {
        if (((tail_index_ + 1) & mask_) == head_index_) {
            return nullptr;  // Pool is full
        }

        T* ret = reinterpret_cast<T*>(blocks_[tail_index_]);
        blocks_[tail_index_] = nullptr;
        tail_index_ = (tail_index_ + 1) & mask_;
        return ret;
    }

    void free(T* block)
    {
        uint32_t free_block = head_index_ & mask_;
        blocks_[free_block] = reinterpret_cast<void*>(block);

        head_index_ = (head_index_ + 1) & mask_;
    }

private:
    uint8_t* memory_;
    uint32_t tail_index_;
    uint32_t head_index_;
    void** blocks_;
    uint32_t mask_;
    uint32_t size_;
    uint32_t count_;
};

/*******************************************************************************
 * @brief 特化bam1_t*
 *
 * @tparam
 * *****************************************************************************
 */

static const uint32_t BAM_DATA_SIZE = 2048;
static const uint32_t BAM_DATA_CAPACITY = sizeof(bam1_t) + BAM_DATA_SIZE;

template <>
bam1_t* RingMemPool<bam1_t>::alloc()
{
    if (((tail_index_ + 1) & mask_) == head_index_) {
        return nullptr;  // Pool is full
    }

    void* res = blocks_[tail_index_];
    blocks_[tail_index_] = nullptr;
    tail_index_ = (tail_index_ + 1) & mask_;

    bam1_t* ret = reinterpret_cast<bam1_t*>(res);
    ret->data = reinterpret_cast<uint8_t*>((uint8_t*)res + sizeof(bam1_t));
    ret->m_data = BAM_DATA_SIZE;
    ret->mempolicy = BAM_USER_OWNS_DATA | BAM_USER_OWNS_STRUCT;
    return ret;
}

#endif  // _RING_MEM_POOL_H_
