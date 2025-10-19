#include "base_quality_rank_sum_test.h"

#include "rovaca_logger.h"
#include "info_data.hpp"
#include "utils/read_record_utils.h"
#include "variant.h"

namespace rovaca
{

void BaseQualityRankSumTest::set_values_to_info(pInfoData info, void *data)
{
    double value = *((double *)data);
    info->set_base_qrank_sum(value);
}

OptionalDouble BaseQualityRankSumTest::get_read_base_quality(pReadRecord read, pVariant vc)
{
    CHECK_CONDITION_EXIT(nullptr == read, "nullptr == read");
    OptionalUint8 ref_bases_qual = ReadRecordUtils::get_read_base_quality_at_reference_coordinate(read, vc->get_start());
    if (ref_bases_qual.first) {
        return {true, (double)ref_bases_qual.second};
    }
    return {false, NAN};
}

}  // namespace rovaca