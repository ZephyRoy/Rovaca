#ifndef ROVACA_HC_INDEX_RANGE_H_
#define ROVACA_HC_INDEX_RANGE_H_
#include <cstdint>
#include <functional>
#include "rovaca_logger.h"
#include "rovaca_logger.h"
#include "genotype_macors.h"

namespace rovaca
{

class IndexRange
{
    /* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
private:
    int32_t _from, _to;

    /* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
public:
    IndexRange(int32_t from, int32_t to)
        : _from(from)
        , _to(to)
    {}

    void shift(int32_t sht)
    {
        _to += sht;
        _from += sht;
        validate();
    }

    void shift_left(int32_t sht) { shift(-sht); }

    void shift_start(int32_t sht)
    {
        _from += sht;
        validate();
    }

    void shift_start_left(int32_t sht) { shift_start(-sht); }

    void shift_end(int32_t sht)
    {
        _to += sht;
        validate();
    }

    void shift_end_left(int32_t sht) { shift_end(-sht); }

    int32_t get_start() const { return _from; }
    int32_t get_end() const { return _to; }
    int32_t size() const { return _to - _from; }

    void for_each(const std::function<void(int32_t)>& func) const
    {
        for (int32_t i = _from; i < _to; ++i) {
            func(i);
        }
    }

    DoubleVector map_to_double(const std::function<double(int32_t)>& int2double_func, pMemoryPool pool) const
    {
        DoubleVector result(size(), pool);
        for (int32_t i = _from; i < _to; ++i) {
            result.at(i - _from) = int2double_func(i);
        }
        return result;
    }

    /**
     * sums the values of an int -> double function applied to this range
     * @param lambda the int -> double function
     */
    double sum(const std::function<double(int32_t)>& int2double_func) const
    {
        double result = 0;
        for (int32_t i = _from; i < _to; i++) {
            result += int2double_func(i);
        }
        return result;
    }

    Int32Vector map_to_integer(const std::function<int32_t(int32_t)>& func, pMemoryPool pool) const
    {
        Int32Vector result(size(), pool);
        for (int32_t i = _from; i < _to; i++) {
            result[i - _from] = func(i);
        }
        return result;
    }

    /* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
private:
    inline void validate() const
    {
        CHECK_CONDITION_EXIT(_from > _to, "the range size cannot be negative");
        CHECK_CONDITION_EXIT(_from < 0, "the range cannot contain negative indices");
    }
};

}  // namespace rovaca

#endif  // ROVACA_HC_INDEX_RANGE_H_
