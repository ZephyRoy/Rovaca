#ifndef ROVACA_HC_CONSTANTS_STR_H_
#define ROVACA_HC_CONSTANTS_STR_H_
#include <cstring>
#include <memory>
#include <memory_resource>
#include <string>

#include "forward.h"
#include "genotype_macors.h"
#include "genotype_struct.h"

namespace rovaca
{

class ConstantsStr
{
    /* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
private:
    static constexpr size_t _buffer_size = 4096;
    uint8_t _buffer[_buffer_size];
    std::pmr::monotonic_buffer_resource* _pool;

    /* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
public:
    // germline_info
    const std::string k_mle_allele_count_key = "MLEAC";
    const std::string k_mle_allele_frequency_key = "MLEAF";
    const std::string k_base_qual_rank_sum_key = "BaseQRankSum";
    const std::string k_allele_count_key = "AC";
    const std::string k_allele_frequency_key = "AF";
    const std::string k_allele_number_key = "AN";
    const std::string k_depth_key = "DP";
    const std::string k_excess_het_key = "ExcessHet";
    const std::string k_fisher_strand_key = "FS";
    const std::string k_inbreeding_coefficient_key = "InbreedingCoeff";
    const std::string k_map_qual_rank_sum_key = "MQRankSum";
    const std::string k_qual_by_depth_key = "QD";
    const std::string k_rms_mapping_quality_key = "MQ";
    const std::string k_raw_mapping_quality_with_depth_key = "RAW_MQandDP";
    const std::string k_read_pos_rank_sum_key = "ReadPosRankSum";
    const std::string k_strand_odds_ratio_key = "SOR";
    const std::string k_end_key = "END";

    // germline_format
    const std::string k_genotype_key = "GT";
    const std::string k_genotype_allele_depths = "AD";
    const std::string k_genotype_quality_key = "GQ";
    const std::string k_haplotype_caller_phasing_gt_key = "PGT";
    const std::string k_haplotype_caller_phasing_id_key = "PID";
    const std::string k_genotype_pl_key = "PL";
    const std::string k_phase_set_key = "PS";
    const std::string k_strand_bias_key = "SB";
    const std::string k_min_dp_key = "MIN_DP";

    // special alleles
    const uint8_t k_spanning_deletion_allele = '*';
    const uint8_t k_no_call_allele = '.';
    const uint8_t k_null_allele = '-';

    // 原来Allele内部的静态数据
    const uint8_t k_single_breakend_indicator = '.';
    const uint8_t k_breakend_extending_right = '[';
    const uint8_t k_breakend_extending_left = ']';
    const uint8_t k_symbolic_allele_start = '<';
    const uint8_t k_symbolic_allele_end = '>';

    const char k_phased = '|';
    const char k_unphased = '/';
    const char k_empty_allele = '.';

    const char k_empty_id_field = '.';
    const char k_field_separator = '\t';
    const char k_missing_value_v4 = '.';
    const char k_info_field_separator = ';';
    const char k_genotype_field_separator = ':';
    const char k_info_field_array_separator = ',';
    const char k_empty_alternate_allele_field = '.';

    const char k_unfiltered = '.';
    const std::string k_passes_filters_v4 = "pass";

    // 预先申请好的Bases
    pBases k_bases_a;
    pBases k_bases_t;
    pBases k_bases_g;
    pBases k_bases_c;
    pBases k_bases_n;
    pBases k_bases_span_del;
    pBases k_bases_no_call;
    pBases k_bases_non_ref_allele;
    pBases k_bases_unspecified_alternate_allele;
    pBases k_spanning_deletion_symbolic_allele_deprecated;
    pBases k_empty_allele_bases;

    ~ConstantsStr()
    {
        if (_pool) {
            delete _pool;
            _pool = nullptr;
        }
    }

    static pConstantsStr get_instance()
    {
        static std::unique_ptr<ConstantsStr> c(new ConstantsStr{});
        return c.get();
    }

    /* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
private:
    ConstantsStr()
        : _buffer()
        , _pool(new std::pmr::monotonic_buffer_resource(_buffer, _buffer_size, std::pmr::null_memory_resource()))
        , k_bases_a(new ALLOC_FLEXIBLE_IN_POOL(_pool, Bases, 1, uint8_t) Bases{1})
        , k_bases_t(new ALLOC_FLEXIBLE_IN_POOL(_pool, Bases, 1, uint8_t) Bases{1})
        , k_bases_g(new ALLOC_FLEXIBLE_IN_POOL(_pool, Bases, 1, uint8_t) Bases{1})
        , k_bases_c(new ALLOC_FLEXIBLE_IN_POOL(_pool, Bases, 1, uint8_t) Bases{1})
        , k_bases_n(new ALLOC_FLEXIBLE_IN_POOL(_pool, Bases, 1, uint8_t) Bases{1})
        , k_bases_span_del(new ALLOC_FLEXIBLE_IN_POOL(_pool, Bases, 1, uint8_t) Bases{1})
        , k_bases_no_call(new ALLOC_FLEXIBLE_IN_POOL(_pool, Bases, 1, uint8_t) Bases{1})
        , k_bases_non_ref_allele(new ALLOC_FLEXIBLE_IN_POOL(_pool, Bases, 9, uint8_t) Bases{9})
        , k_bases_unspecified_alternate_allele(new ALLOC_FLEXIBLE_IN_POOL(_pool, Bases, 3, uint8_t) Bases{3})
        , k_spanning_deletion_symbolic_allele_deprecated(new ALLOC_FLEXIBLE_IN_POOL(_pool, Bases, 7, uint8_t) Bases{7})
        , k_empty_allele_bases(new ALLOC_FLEXIBLE_IN_POOL(_pool, Bases, 0, uint8_t) Bases{0})
    {
        k_bases_a->data[0] = 'A';
        k_bases_t->data[0] = 'T';
        k_bases_g->data[0] = 'G';
        k_bases_c->data[0] = 'C';
        k_bases_n->data[0] = 'N';
        k_bases_span_del->data[0] = '*';
        k_bases_no_call->data[0] = '.';
        memcpy(k_bases_non_ref_allele->data, "<NON_REF>", sizeof(uint8_t) * strlen("<NON_REF>"));
        memcpy(k_bases_unspecified_alternate_allele->data, "<*>", sizeof(uint8_t) * strlen("<*>"));
        memcpy(k_spanning_deletion_symbolic_allele_deprecated->data, "<*:DEL>", sizeof(uint8_t) * strlen("<*:DEL>"));
    }
};

}  // namespace rovaca

#endif  // ROVACA_HC_CONSTANTS_STR_H_
