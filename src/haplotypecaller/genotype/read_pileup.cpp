#include "read_pileup.h"

#include "pileup_element.h"

namespace rovaca
{

pReadPileup ReadPileup::create(int32_t contig, int64_t start, int64_t stop, pMemoryPool pool)
{
    return new ALLOC_TYPE_IN_POOL(pool, ReadPileup) ReadPileup{contig, start, stop, pool};
}

ReadPileup::ReadPileup(int32_t contig, int64_t start, int64_t stop, pMemoryPool pool)
    : InterfaceLocatable()
    , _tid(contig)
    , _start(start)
    , _stop(stop)
    , _pileup_elements(pool)
    , qual_his()
    , qual_max()
{}

}  // namespace rovaca