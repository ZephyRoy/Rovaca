#ifndef DEV_CODE_ROVACA_MEMORY_POOL_H_
#define DEV_CODE_ROVACA_MEMORY_POOL_H_
#include <memory_resource>

namespace rovaca
{

class RovacaMemoryPool : public std::pmr::memory_resource
{
    friend class MemoryPoolGuard;

    /* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
private:
    void *cache_;
    size_t capacity_;
    std::pmr::memory_resource *additional_mem_;

    /* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
public:
    RovacaMemoryPool(uint8_t *cache, size_t capacity)
        : std::pmr::memory_resource{}
        , cache_{cache}
        , capacity_{capacity}
        , additional_mem_{nullptr}
    {}

    ~RovacaMemoryPool() override { delete additional_mem_; }

    /* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
private:
    void *do_allocate(size_t _bytes, size_t _alignment) override;
    void do_deallocate([[maybe_unused]] void *_p, [[maybe_unused]] size_t _bytes, [[maybe_unused]] size_t _alignment) override {}
    bool do_is_equal([[maybe_unused]] const memory_resource &_other) const noexcept override { return false; }

    /*!
     * @brief 以下两个函数配合 MemoryPoolGuard 使用, 对内存的精细化控制
     */
    std::pair<void *, std::size_t> mark() { return {cache_, capacity_}; }
    void reset(const std::pair<void *, std::size_t> &old)
    {
        cache_ = old.first;
        capacity_ = old.second;
    }
};

/*!
 * @brief 记录 RovacaMemoryPool 当前使用的 size，离开作用域的时候重置，精细化使用内存
 */
class MemoryPoolGuard
{
    /* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
private:
    std::pair<void *, std::size_t> old_;
    RovacaMemoryPool *pool_;

    /* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
public:
    MemoryPoolGuard(RovacaMemoryPool *pool)
        : old_(pool->mark())
        , pool_(pool)
    {}
    ~MemoryPoolGuard() { pool_->reset(old_); }
    MemoryPoolGuard(const RovacaMemoryPool &) = delete;
    void operator=(const RovacaMemoryPool &) = delete;
};

}  // namespace rovaca

#endif  // DEV_CODE_ROVACA_MEMORY_POOL_H_
