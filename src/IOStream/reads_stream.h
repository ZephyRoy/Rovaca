#ifndef READS_ITERATOR_H
#define READS_ITERATOR_H

#include <htslib/sam.h>

#include "../common/reads_filter.h"
#include "bam_loader.h"
#include "ring_mem_pool.hpp"

enum stream_type { raw = 0, filtered = 1, downsampled = 2, transformed = 3 };
class ReadStream
{
    BamLoader* bam_loader_;
    stream_type stream_type_;

public:
    ReadStream(BamLoader* bam_loader) noexcept
        : bam_loader_(bam_loader)
        , stream_type_(stream_type::raw)
    {}
    ReadStream(ReadStream* other) noexcept
        : bam_loader_(other->bam_loader_)
    {}
    virtual ~ReadStream() {}
    virtual bool has_next() { return bam_loader_->has_next(); }
    virtual bam1_t* next(int& file_index)
    {
        while (has_next()) {
            bam1_t* item = bam_loader_->get_next_read(file_index);
            if (!item) continue;
            return item;
        }
        return nullptr;
    }
    virtual void set_type(stream_type type) { stream_type_ = type; }
    virtual stream_type get_type() { return stream_type_; }
    virtual void mem_recovery(bam1_t* read) { bam_loader_->read_recovery(read); }
    virtual void deallocate_bqsr_mem(__attribute((unused)) bam1_t* bam) {}  // Should not call parent method here!
    inline void set_target(const char* target) { bam_loader_->set_target(target); }
    inline void set_target(const p_bed_intervals target) { bam_loader_->set_target(target); }
};

class ReadsFilterIterator : public ReadStream
{
private:
    ReadFilter* read_filter_;
    ReadStream* reads_iterator_;

public:
    ReadsFilterIterator(ReadStream* reads_iterator, ReadFilter* read_filter)
        : ReadStream(reads_iterator)
        , read_filter_(read_filter)
        , reads_iterator_(reads_iterator)
    {
        set_type(stream_type::filtered);
    }
    ~ReadsFilterIterator() {}
    bool has_next() override { return reads_iterator_->has_next(); }
    bam1_t* next(int& file_index) override
    {
        while (has_next()) {
            bam1_t* item = reads_iterator_->next(file_index);
            if (item == NULL || item->data == NULL) continue;
            if (read_filter_->test(item))
                return item;
            else
                reads_iterator_->mem_recovery(item);
        }
        return nullptr;
    }
};

#include <condition_variable>
#include <memory>
#include <mutex>
#include <queue>
#include <utility>

template <typename T>
class TaskQueue
{
public:
    explicit TaskQueue(size_t max_size)
        : max_size_(max_size)
    {}

    void push(T value)
    {
        std::unique_lock<std::mutex> lock(mutex_);
        cond_not_full_.wait(lock, [this] { return queue_.size() < max_size_; });
        queue_.push(std::move(value));
        cond_not_empty_.notify_one();
    }

    bool pop(T& value)
    {
        std::unique_lock<std::mutex> lock(mutex_);
        if (queue_.empty()) return false;
        cond_not_empty_.wait(lock, [this] { return !queue_.empty(); });
        value = std::move(queue_.front());
        queue_.pop();
        cond_not_full_.notify_one();
        return true;
    }

    T pop()
    {
        std::unique_lock<std::mutex> lock(mutex_);
        cond_not_empty_.wait(lock, [this] { return !queue_.empty(); });
        T value = queue_.front();
        queue_.pop();
        cond_not_full_.notify_one();
        return value;
    }

    bool empty() const
    {
        std::unique_lock<std::mutex> lock(mutex_);
        return queue_.empty();
    }

    void wait_for_task()
    {
        std::unique_lock<std::mutex> lock(mutex_);
        cond_not_empty_.wait(lock, [this] { return !queue_.empty(); });
    }

private:
    size_t max_size_;
    mutable std::mutex mutex_;
    std::queue<T> queue_;
    std::condition_variable cond_not_empty_;
    std::condition_variable cond_not_full_;
};

#include <atomic>
#include <chrono>
#include <list>

struct ReducedQueue
{
private:
    std::unique_ptr<std::list<bam1_t*>> processed_reads_;
    std::mutex mtx;
    std::atomic_bool done;

public:
    ReducedQueue()
        : processed_reads_(std::make_unique<std::list<bam1_t*>>())
        , done(false)
    {}

    void push(bam1_t* item)
    {
        std::unique_lock<std::mutex> lock(mtx);
        processed_reads_->push_back(item);
    }

    bam1_t* pop()
    {
        std::unique_lock<std::mutex> lock(mtx);
        if (processed_reads_->empty()) {
            return nullptr;  // Return nullptr instead of undefined behavior
        }
        bam1_t* ret = processed_reads_->front();
        processed_reads_->pop_front();
        return ret;
    }

    void splice(std::list<bam1_t*>& list_item)
    {
        std::unique_lock<std::mutex> lock(mtx);
        processed_reads_->splice(processed_reads_->end(), list_item);
    }

    bool empty()
    {
        std::unique_lock<std::mutex> lock(mtx);
        bool empty = processed_reads_->empty();
        lock.unlock();
        return empty;
    }
    void emit_done_signal() { done.store(true, std::memory_order_release); }
    bool reads_done() { return done.load(std::memory_order_acquire); }
};

#include <memory_resource>
class ReadsTransformerIterator : public ReadStream
{
private:
    ReducedQueue* reduced_queue_;
    RingMemPool<bam1_t>* mem_;

public:
    ReadsTransformerIterator(ReadStream* reads_iterator, ReducedQueue* reduced_queue, RingMemPool<bam1_t>* mem)
        : ReadStream(reads_iterator)
        , reduced_queue_(reduced_queue)
        , mem_(mem)
    {
        set_type(stream_type::transformed);
    }
    ~ReadsTransformerIterator() {}

    bool has_next() override { return !reduced_queue_->reads_done() || !reduced_queue_->empty(); }

    bam1_t* next(__attribute((unused)) int& file_index) override { return reduced_queue_->pop(); }
    void deallocate_bqsr_mem(bam1_t* bam) override
    {
        if (!bam_get_mempolicy(bam)) {
            bam_destroy1(bam);
            return;
        }
        mem_->free(bam);
    }
};

#include "../common/downsampler.h"

class ReadsDownsampleIterator : public ReadStream
{
private:
    Downsampler* downsampler_;
    ReadStream* reads_iterator_;
    bam1_t* next_item;
    std::list<bam1_t*> cache;

public:
    ReadsDownsampleIterator(ReadStream* reads_iterator, Downsampler* downsampler)
        : ReadStream(reads_iterator)
        , downsampler_(downsampler)
        , reads_iterator_(reads_iterator)
        , next_item(nullptr)
        , cache({})
    {
        set_type(stream_type::downsampled);
    }
    ~ReadsDownsampleIterator() {}
    bool has_next() override { return reads_iterator_->has_next() && !downsampler_->has_finalized_items(); }
    bam1_t* next(int& file_index) override
    {
    file_index = 0;  // Temporary.
        load_next_item();
        bam1_t* ret = next_item;
        return ret;
    }

    bool load_next_item()
    {
        if (!cache.empty() || (cache.empty() && apply())) {
            next_item = cache.front();
            cache.pop_front();
        }
        else {
            next_item = nullptr;
        }
    // Prefer to pop reads from cache first.
        return next_item != nullptr;
    }

    bool apply()
    {
        while (has_next()) {
            int file_index;
            auto item = reads_iterator_->next(file_index);
            if (item == NULL || item->data == NULL) continue;
            // This interface returns reads that need to be downsampled.
            auto ret = downsampler_->submit(item);
            if (ret) reads_iterator_->mem_recovery(ret);
        }

        if (!reads_iterator_->has_next()) {
            downsampler_->input_end_signal();
        }
        downsampler_->consume_finalized_items(cache);
        return !cache.empty();
    }
};
#endif  // READS_ITERATOR_H