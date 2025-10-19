#ifndef ROVACA_HC_READ_CLIPPER_H_
#define ROVACA_HC_READ_CLIPPER_H_
#include <vector>

#include "clipping_op.h"
#include "forward.h"
#include "genotype_macors.h"

namespace rovaca
{

class ReadClipper
{
    /* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
private:
    pReadRecord _read;
    pMemoryPool _pool;
    pBamDataPool _bam_pool;
    std::pmr::vector<ClippingOp> _ops;

    /* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
public:
    DISALLOW_COPY_AND_ASSIGN(ReadClipper);
    ReadClipper(pReadRecord read, pMemoryPool pool, pBamDataPool bampool)
        : _read(read)
        , _pool(pool)
        , _bam_pool(bampool)
        , _ops(pool)
    {}

    /*! @brief Hard clip the read to the variable region (from refStart to refStop) */
    pReadRecord hard_clip_to_region(int64_t ref_start, int64_t ref_stop);

    pReadRecord clip_read(ClippingRepresentation algorithm);

    /* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
    // private:
    /*!
     * @brief 裁剪read，通过ref_start<0或者ref_stop<0指示裁剪的方向
     * @note 仅供hard_clip_by_reference_coordinates_left_tail和hard_clip_by_reference_coordinates_right_tail调用
     */
    pReadRecord clip_by_reference_coordinates(int64_t ref_start, int64_t ref_stop, ClippingRepresentation op);

    /*! @brief Hard clips the left tail of a read up to (and including) refStop using reference coordinates */
    pReadRecord hard_clip_by_reference_coordinates_left_tail(int64_t ref_stop)
    {
        return clip_by_reference_coordinates(-1, ref_stop, HARDCLIP_BASES);
    }
    /*! @brief Hard clips the right tail of a read starting at (and including) refStart using reference coordinates. */
    pReadRecord hard_clip_by_reference_coordinates_right_tail(int64_t ref_start)
    {
        return clip_by_reference_coordinates(ref_start, -1, HARDCLIP_BASES);
    }

    /*!
     * @brief Hard clips both tails of a read.
     *           Left tail goes from the beginning to the 'left' coordinate (inclusive)
     *           Right tail goes from the 'right' coordinate (inclusive) until the end of the read
     */
    pReadRecord hard_clip_both_ends_by_reference_coordinates(int64_t left, int64_t right);

    pReadRecord hard_clip_soft_clipped_bases();
};

}  // namespace rovaca

#endif  // ROVACA_HC_READ_CLIPPER_H_
