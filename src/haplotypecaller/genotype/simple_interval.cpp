#include "simple_interval.h"

#include "rovaca_logger.h"

namespace rovaca
{

pSimpleInterval SimpleInterval::create(int32_t tid, int64_t start, int64_t stop, pMemoryPool pool)
{
    return new ALLOC_TYPE_IN_POOL(pool, SimpleInterval) SimpleInterval{tid, start, stop};
}

pSimpleInterval SimpleInterval::create(const InterfaceLocatable& loc, pMemoryPool pool)
{
    return new ALLOC_TYPE_IN_POOL(pool, SimpleInterval) SimpleInterval{loc};
}

pSimpleInterval SimpleInterval::create(int32_t tid, int64_t start, int64_t stop) { return new SimpleInterval{tid, start, stop}; }

pSimpleInterval SimpleInterval::expand_within_contig(int64_t padding, int64_t sequence_length, pMemoryPool pool) const
{
    CHECK_CONDITION_EXIT(padding < 0, "padding must be >= 0");
    int64_t bounded_start = std::max((int64_t)1, _start - padding);
    int64_t bounded_stop = std::min(sequence_length, _stop + padding);

    if (bounded_start > sequence_length || bounded_stop < 1) {
        return nullptr;
    }
    else {
        return new ALLOC_TYPE_IN_POOL(pool, SimpleInterval) SimpleInterval{_tid, bounded_start, bounded_stop};
    }
}

pSimpleInterval SimpleInterval::intersect(const InterfaceLocatable& other, pMemoryPool pool) const
{
    CHECK_CONDITION_EXIT(!this->overlaps(other), "the two intervals need to overlap");
    return SimpleInterval::create(this->get_tid(), std::max(this->get_start(), other.get_start()),
                                  std::min(this->get_stop(), other.get_stop()), pool);
}

}  // namespace rovaca
