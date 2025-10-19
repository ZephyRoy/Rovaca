#include "read_clipper.h"

#include "bam_data_pool.hpp"
#include "rovaca_logger.h"
#include "genotype_struct.h"
#include "read_record.h"
#include "utils/cigar_utils.h"
#include "utils/read_record_utils.h"

namespace rovaca
{

pReadRecord ReadClipper::hard_clip_to_region(int64_t ref_start, int64_t ref_stop)
{
    int64_t read_start = _read->get_start();
    int64_t read_stop = _read->get_stop();
    // check if the read is contained in region
    if (read_start <= ref_stop && read_stop >= ref_start) {
        if (read_start < ref_start && read_stop > ref_stop) {
            return hard_clip_both_ends_by_reference_coordinates(ref_start - 1, ref_stop + 1);
        }
        else if (read_start < ref_start) {
            return hard_clip_by_reference_coordinates_left_tail(ref_start - 1);
        }
        else if (read_stop > ref_stop) {
            return hard_clip_by_reference_coordinates_right_tail(ref_stop + 1);
        }
        return ReadRecordUtils::copy_read(_read, _bam_pool, _pool);
    }
    else {
        // TODO: 返回emptyRead会因为seq长度为0被过滤，索性不构建新read
        // return ReadUtils.emptyRead(read);
        return nullptr;
    }
}

pReadRecord ReadClipper::clip_read(ClippingRepresentation algorithm)
{
    if (_ops.empty()) {
        return _read;
    }

    pReadRecord clipped_read = _read;
    for (const ClippingOp& op : _ops) {
        int64_t read_len = clipped_read->seq_length();
        if (op._start < read_len) {
            ClippingOp fixed_op = op;
            if (op._start >= read_len) {
                fixed_op = ClippingOp{op._start, read_len - 1, _pool, _bam_pool};
            }
            clipped_read = fixed_op.apply(clipped_read, algorithm);
        }
    }
    _ops.clear();
    return (nullptr == clipped_read || clipped_read->is_empty()) ? nullptr : clipped_read;
}

pReadRecord ReadClipper::clip_by_reference_coordinates(int64_t ref_start, int64_t ref_stop, ClippingRepresentation op)
{
    if (_read->is_empty()) {
        return _read;
    }

    CHECK_CONDITION_EXIT(op == SOFTCLIP_BASES && _read->is_unmapped(),
                         "cannot soft-clip read by reference coordinates because it is unmapped");

    int64_t start, stop;

    if (ref_start < 0) {
        CHECK_CONDITION_EXIT(ref_stop < 0, "only one of ref_start or ref_stop must be < 0, not both ({}, {})", ref_start, ref_stop);
        start = 0;
        Int64Uint32Pair stop_pos_and_op = ReadRecordUtils::get_read_index_for_reference_coordinate(_read, ref_stop);

        // if the refStop falls in a deletion, the above method returns the position after the deletion
        // Since the stop we return here is inclusive, we decrement the stop to avoid overclipping by one base
        // As a result we do not clip the deletion, which is fine.
        stop = stop_pos_and_op.first - (consumes_read_bases(stop_pos_and_op.second) ? 0 : 1);
    }
    else {
        CHECK_CONDITION_EXIT(ref_stop >= 0, "either ref_start or ref_stop must be < 0 ({}, {})", ref_start, ref_stop);
        // unlike the above case where we clip the start fo the read, here we clip the end and returning the base to the right of a deletion
        // avoids overclipping
        start = ReadRecordUtils::get_read_index_for_reference_coordinate(_read, ref_start).first;
        stop = _read->seq_length() - 1;
    }

    if (start == ReadRecordUtils::s_read_index_not_found || stop == ReadRecordUtils::s_read_index_not_found) {
        return _read;
    }

    CHECK_CONDITION_EXIT(start < 0 || stop > _read->seq_length() - 1, "trying to clip before the start or after the end of a read");
    CHECK_CONDITION_EXIT(start > stop, "START ({}) > ({}) STOP -- this should never happen, please check read", start, stop);
    CHECK_CONDITION_EXIT(start > 0 && stop < _read->seq_length() - 1, "trying to clip the middle of the read: start {}, stop {}", start,
                         stop);

    _ops.emplace_back(start, stop, _pool, _bam_pool);

    pReadRecord clipped_read = clip_read(op);

    _ops.clear();

    return clipped_read;
}

pReadRecord ReadClipper::hard_clip_both_ends_by_reference_coordinates(int64_t left, int64_t right)
{
    if (_read->is_empty() || left == right) {
        return nullptr;
    }

    pReadRecord left_tail_read = clip_by_reference_coordinates(right, -1, HARDCLIP_BASES);
    if (left > left_tail_read->get_stop()) {
        return nullptr;
    }

    ReadClipper clipper(left_tail_read, _pool, _bam_pool);
    return clipper.hard_clip_by_reference_coordinates_left_tail(left);
}

pReadRecord ReadClipper::hard_clip_soft_clipped_bases()
{
    if (_read->is_empty()) {
        return _read;
    }

    int32_t read_index = 0;
    int32_t cut_left = -1;
    int32_t cut_right = -1;
    bool right_trail = false;

    for (size_t c = 0; c < _read->cigar_length(); c++) {
        uint32_t cigar_element = _read->cigar_i(c);
        uint32_t cigar_op = bam_cigar_op(cigar_element);
        uint32_t cigar_op_len = bam_cigar_oplen(cigar_element);
        if (cigar_op_is_soft_clip(cigar_op)) {
            if (right_trail) {
                cut_right = read_index;
            }
            else {
                cut_left = read_index + int32_t(cigar_op_len) - 1;
            }
        }
        else if (!cigar_op_is_hard_clip(cigar_op)) {
            right_trail = true;
        }

        if (consumes_read_bases(cigar_op)) {
            read_index += int32_t(cigar_op_len);
        }
    }

    // It is extremely important that we cut the end first otherwise the read coordinates change.
    if (cut_right >= 0) {
        ClippingOp c0(cut_right, _read->seq_length() - 1, _pool, _bam_pool);
        this->_ops.push_back(c0);
    }

    if (cut_left >= 0) {
        ClippingOp c1(0, cut_left, _pool, _bam_pool);
        this->_ops.push_back(c1);
    }

    return clip_read(HARDCLIP_BASES);
}

}  // namespace rovaca