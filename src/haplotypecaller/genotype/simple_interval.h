#ifndef ROVACA_HC_SIMPLE_INTERVAL_H_
#define ROVACA_HC_SIMPLE_INTERVAL_H_
#include "forward.h"
#include "genotype_macors.h"
#include "interface/interface_locatable.hpp"

namespace rovaca
{

class SimpleInterval : public InterfaceLocatable
{
    /* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
private:
    int32_t _tid{INVALID_INT};
    int64_t _start{INVALID_INT}, _stop{INVALID_INT};

public:
    DISALLOW_COPY_AND_ASSIGN(SimpleInterval);
    static pSimpleInterval create(const InterfaceLocatable& loc, pMemoryPool pool);
    static pSimpleInterval create(int32_t tid, int64_t start, int64_t stop, pMemoryPool pool);
    static pSimpleInterval create(int32_t tid, int64_t start, int64_t stop);

    int32_t get_tid() const override { return _tid; }
    int64_t get_start() const override { return _start; }
    int64_t get_stop() const override { return _stop; }
    void set_tid(int32_t tid) override { _tid = tid; }
    void set_start(int64_t start) override { _start = start; }
    void set_stop(int64_t stop) override { _stop = stop; }
    bool overlaps(const InterfaceLocatable& other) const override { return overlaps_with_margin(other, 0); }

    int64_t get_length() const { return _stop - _start + 1; }

    int64_t size() const { return _stop - _start + 1; }

    /*!
     * @brief 用于确定该区间是否与另一个Locatable对象重叠，并且两个区间之间的距离不超过指定的margin
     */
    bool overlaps_with_margin(const InterfaceLocatable& other, int64_t margin) const
    {
        return _tid == other.get_tid() && this->_start <= other.get_stop() + margin && this->_stop >= other.get_start() - margin;
    }

    /*!
     * @brief 将一个区间向两个方向扩展指定的距离
     */
    pSimpleInterval expand_within_contig(int64_t padding, int64_t sequence_length, pMemoryPool pool) const;

    /*!
     * @brief 计算两个区间的交集
     */
    pSimpleInterval intersect(const InterfaceLocatable& other, pMemoryPool pool) const;

    /* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
private:
    SimpleInterval() = default;
    ~SimpleInterval() override = default;
    SimpleInterval(int32_t tid, int64_t start, int64_t stop)
        : InterfaceLocatable()
        , _tid(tid)
        , _start(start)
        , _stop(stop)
    {}
    explicit SimpleInterval(const InterfaceLocatable& loc)
        : InterfaceLocatable()
        , _tid(loc.get_tid())
        , _start(loc.get_start())
        , _stop(loc.get_stop())
    {}
};

}  // namespace rovaca

#endif  // ROVACA_HC_SIMPLE_INTERVAL_H_
