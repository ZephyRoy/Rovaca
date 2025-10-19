#ifndef ROVACA_HC_CLIPPING_OP_H_
#define ROVACA_HC_CLIPPING_OP_H_
#include "forward.h"
#include "genotype_enum.h"

namespace rovaca
{

class ClippingOp
{
    /* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
public:
    int64_t _start, _stop;  // inclusive
    pMemoryPool _pool;
    pBamDataPool _bam_pool;

    ClippingOp(int64_t start, int64_t stop, pMemoryPool pool, pBamDataPool bam_pool)
        : _start(start)
        , _stop(stop)
        , _pool(pool)
        , _bam_pool(bam_pool)
    {}

    /*!
     * @brief Clips the bases in read according to this operation's start and stop.
     * Uses the clipping representation used is the one provided by algorithm argument.
     */
    pReadRecord apply(const pReadRecord& original, ClippingRepresentation op);

    /* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
private:
    /*!
     * Hard clip bases from read, from start to stop in base coordinates
     * <p>
     * If start == 0, then we will clip from the front of the read, otherwise we clip from the right.  If start == 0 and stop == 10, this
     * would clip out the first 10 bases of the read.
     * <p>
     * Note that this function works with reads with negative alignment starts, in order to allow us to hardClip reads that have had their
     * soft clips reverted and so might have negative alignment starts
     * <p>
     * Works properly with reduced reads and insertion/deletion base qualities
     * <p>
     * Note: this method does not assume that the read is directly modifiable and makes a copy of it.
     *
     * @param read  a non-null read
     * @param start a start >= 0 and < read.length
     * @param stop  a stop >= 0 and < read.length.
     * @return a cloned version of read that has been properly trimmed down (Could be an empty, unmapped read)
     */
    pReadRecord apply_hard_clip_bases(const pReadRecord& read, int64_t start, int64_t stop) const;
};

}  // namespace rovaca

#endif  // ROVACA_HC_CLIPPING_OP_H_
