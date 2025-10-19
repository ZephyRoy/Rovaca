#include "rms_mapping_quality.h"

#include "allele_likelihoods.hpp"
#include "rovaca_logger.h"
#include "info_data.hpp"
#include "quality_utils.h"
#include "read_record.h"

namespace rovaca
{

void RMSMappingQuality::annotate([[maybe_unused]] pRefFragment ref, [[maybe_unused]] pVariant vc, pRALikelihoods likelihoods,
                                 pInfoData target, pMemoryPool pool)
{
    if (likelihoods->evidence_count() < 1) {
        return;
    }

    long square_sum = 0;
    long num_reads_used = 0;
    uint8_t mq;
    size_t si, sample_count;
    for (si = 0, sample_count = likelihoods->number_of_samples(); si < sample_count; ++si) {
        const auto& sample_evidences = likelihoods->sample_evidence(si);
        for (auto r : sample_evidences) {
            mq = r->mapping_quality();
            if (mq != QualityUtils::s_mapping_quality_unavailable) {
                square_sum += mq * mq;
                num_reads_used++;
            }
        }
    }

    if (_use_raw) {
        // gvcf
        std::pmr::vector<long> RAW_MQandDP({square_sum, num_reads_used}, pool);
        target->set_raw_mq_and_dp(std::move(RAW_MQandDP));
    }
    else {
        // vcf
        double MQ = std::sqrt(double(square_sum) / double(num_reads_used));
        target->set_mq(MQ);
    }
}

}  // namespace rovaca