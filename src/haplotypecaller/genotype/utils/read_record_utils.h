#ifndef ROVACA_HC_READ_RECORD_UTILS_H_
#define ROVACA_HC_READ_RECORD_UTILS_H_
#include <cstdint>
#include <utility>

#include "forward.h"

namespace rovaca
{

namespace ReadRecordUtils
{

static constexpr int64_t s_read_index_not_found = -1;

pReadRecord copy_read(const pReadRecord& original, uint8_t* seq, uint8_t* qual, uint32_t seq_len, pCigar cigar, pBamDataPool bampool,
                      pMemoryPool gpool);
pReadRecord copy_read(const pReadRecord& original, pBamDataPool bampool, pMemoryPool gpool);

/*!
 * @brief Find the 0-based index within a read base array corresponding to a given 1-based position in the reference, along with the cigar
 * operator of the element containing that base.
 * If the reference coordinate occurs within a deletion, the first index after the deletion is returned.
 * @note Note that this treats soft-clipped bases as if they align with the reference, which is useful for hard-clipping reads with soft
 * clips.
 * @param alignment_start The soft start of the read on the reference
 * @param cigar_num
 * @param cigar
 * @param ref_coord The target reference coordinate
 * @return if the reference coordinate falls within an alignment block of the read's cigar, the corresponding read coordinate;
 *         -1 If the reference coordinate occurs before the read start or after the read end.
 *         -1 if the reference coordinate falls within a deletion, the first read coordinate after the deletion.  Note: if the last cigar
 *         element is a deletion (which isn't meaningful)
 */

Int64Uint32Pair get_read_index_for_reference_coordinate(int64_t alignment_start, uint32_t cigar_num, const uint32_t* cigar,
                                                        int64_t ref_coord);

Int64Uint32Pair get_read_index_for_reference_coordinate(const pReadRecord& read, int64_t ref_coord);

/*!
 * @brief pull out the sample names from a samfile_header
 * @note we use a tree_set so that they are sorted
 */
std::vector<std::string> get_samples_from_header(const std::vector<bam_hdr_t*>& headers);

OptionalUint8 get_read_base_quality_at_reference_coordinate(pReadRecord read, int64_t ref_coord);

}  // namespace ReadRecordUtils

}  // namespace rovaca

#endif  // ROVACA_HC_READ_RECORD_UTILS_H_
