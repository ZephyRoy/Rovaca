#include "rovaca_memory_pool.h"

#define likely(x)   __builtin_expect(!!(x), 1)
#define unlikely(x) __builtin_expect(!!(x), 0)

namespace rovaca
{

void* RovacaMemoryPool::do_allocate(size_t _bytes, size_t _alignment)
{
    if (_alignment < 2 || (_alignment & (_alignment - 1)) != 0) {
        return nullptr;
    }

    void* ret;
    if (std::align(_alignment, _bytes, cache_, capacity_)) {
        ret = cache_;
        cache_ = (uint8_t*)cache_ + _bytes;
        capacity_ -= _bytes;
    }
    else {
        if (unlikely(!additional_mem_)) {
            additional_mem_ = new std::pmr::monotonic_buffer_resource{std::pmr::new_delete_resource()};
        }
        ret = additional_mem_->allocate(_bytes, _alignment);
    }

    return ret;
}

}  // namespace rovaca