#include <chrono>
#include <condition_variable>
#include <mutex>
#include <queue>
#ifndef HC_BLOCKING_QUEUE
#define HC_BLOCKING_QUEUE
template <typename T>
class BlockingQueue
{
public:
    BlockingQueue(size_t fixed_size)
        : m_max_queue_size(fixed_size)
    {}

    void push(const T& value)
    {
        std::unique_lock<std::mutex> lock(m_mutex_);
        while (m_queue_.size() == m_max_queue_size) {
            m_not_full_condition_.wait(lock);
        }
        m_queue_.push(value);
        lock.unlock();
        m_empty_condition_.notify_one();
    }

    T pop(uint32_t sec)
    {
        T value = nullptr;
        std::chrono::seconds tt{sec};
        std::unique_lock<std::mutex> lock(m_mutex_);
        if (m_empty_condition_.wait_for(lock, tt, [&]() { return !m_queue_.empty(); })) {
            value = m_queue_.front();
            m_queue_.pop();
            m_not_full_condition_.notify_one();
        }
        return value;
    }

    T pop()
    {
        std::unique_lock<std::mutex> lock(m_mutex_);
        while (m_queue_.empty()) {
            m_empty_condition_.wait(lock);
        }
        T value = m_queue_.front();
        m_queue_.pop();
        m_not_full_condition_.notify_one();
        return value;
    }

    bool empty()
    {
        std::unique_lock<std::mutex> lock(m_mutex_);
        return m_queue_.empty();
    }

    size_t size()
    {
        std::unique_lock<std::mutex> lock(m_mutex_);
        return m_queue_.size();
    }

private:
    size_t m_max_queue_size;
    std::queue<T> m_queue_;
    std::mutex m_mutex_;
    std::condition_variable m_empty_condition_;
    std::condition_variable m_not_full_condition_;
};
#endif