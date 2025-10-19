#include "read_pos_rank_sum_test.h"

#include <algorithm>

#include "rovaca_logger.h"
#include "info_data.hpp"
#include "read_record.h"
#include "utils/cigar_utils.h"
#include "utils/read_record_utils.h"
#include "variant.h"

namespace rovaca
{

void ReadPosRankSumTest::set_values_to_info(pInfoData info, void* data)
{
    double value = *((double*)data);
    info->set_read_pos_rank_sum(value);
}

OptionalDouble ReadPosRankSumTest::get_read_position(pReadRecord read, pVariant vc)
{
    CHECK_CONDITION_EXIT(nullptr == read, "nullptr == read");

    uint32_t* cigar = read->cigar();
    uint32_t cigar_num = read->cigar_length();
    if (read->get_start() == vc->get_stop() + 1) {
        bool first_non_clipping_is_ins = false;
        uint32_t op;
        for (uint32_t i = 0; i < cigar_num; ++i) {
            op = bam_cigar_op(cigar[i]);
            if (!cigar_op_is_clipping(op)) {
                if (cigar_op_is_ins(op)) {
                    first_non_clipping_is_ins = true;
                }
                break;
            }
        }
        if (first_non_clipping_is_ins) {
            return {true, 0.0};
        }
    }

    auto offset = ReadRecordUtils::get_read_index_for_reference_coordinate(read, vc->get_start());
    if (offset.first == ReadRecordUtils::s_read_index_not_found) {
        return {false, NAN};
    }

    uint32_t first_element = cigar[0];
    uint32_t last_element = cigar[cigar_num - 1];
    uint32_t leading_hard_clips = cigar_op_is_hard_clip(bam_cigar_op(first_element)) ? bam_cigar_oplen(first_element) : 0;
    uint32_t trailing_hard_clips = cigar_op_is_hard_clip(bam_cigar_op(last_element)) ? bam_cigar_oplen(last_element) : 0;
    uint32_t left_distance = leading_hard_clips + (uint32_t)offset.first;
    uint32_t right_distance = (uint32_t)read->seq_length() - 1 - (uint32_t)offset.first + trailing_hard_clips;

    return {true, (double)std::min(left_distance, right_distance)};
}

}  // namespace rovaca