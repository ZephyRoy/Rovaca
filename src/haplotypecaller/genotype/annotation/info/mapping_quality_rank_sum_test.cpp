#include "mapping_quality_rank_sum_test.h"

#include "rovaca_logger.h"
#include "genotype_macors.h"
#include "info_data.hpp"
#include "read_record.h"

namespace rovaca
{

void MappingQualityRankSumTest::set_values_to_info(pInfoData info, void *data)
{
    double value = *((double *)data);
    info->set_mq_rank_sum(value);
}

OptionalDouble MappingQualityRankSumTest::get_element_for_read(pReadRecord read, [[maybe_unused]] pVariant vc)
{
    return {true, (double)read->mapping_quality()};
}

}  // namespace rovaca