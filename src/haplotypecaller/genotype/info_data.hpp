#ifndef ROVACA_HC_INFO_DATA_H_
#define ROVACA_HC_INFO_DATA_H_
#include <bitset>

#include "constants_str.hpp"
#include "forward.h"
#include "genotype_macors.h"
#include "htslib/vcf.h"
#include "rovaca_logger.h"
#include "utils/debug_utils.h"

namespace rovaca
{

#define INFO_FLAG_NUM            (32)
#define INFO_DATA_EXIST_FLAG_NUM (14)

enum InfoDataExistFlag {
    IDEF_MLEAC = 0,
    IDEF_MLEAF,
    IDEF_BaseQualityRankSumTest,
    IDEF_ChromosomeCounts,
    IDEF_Coverage,
    IDEF_ExcessHet,
    IDEF_FisherStrand,
    IDEF_InbreedingCoeff,
    IDEF_MappingQualityRankSumTest,
    IDEF_QualByDepth,
    IDEF_RMSMappingQuality_MQ,
    IDEF_RMSMappingQuality_RAW_MQandDP,
    IDEF_ReadPosRankSumTest,
    IDEF_StrandOddsRatio
};

/*!
 * @brief vcf文件INFO列存储, 有Info时由MemoryPool申请内存,将指针传递给Variant
 * @note 仅实际变异点有 Info 信息，gvcf 模式 NON_REF 点无info，仅有 Variant::_block_end
 * @param _mleac: MLEAC 最大似然估计的等位基因计数
 * @param _mleaf: MLEAF 最大似然估计的等位基因频率
 * @param _base_qrank_sum: BaseQualityRankSumTest 用于检测参考和变异碱基的质量分数之间的差异
 * @param _an: ChromosomeCounts 每个等位基因在每个染色体上的计数
 * @param _ac: ChromosomeCounts 每个等位基因在每个染色体上的计数
 * @param _af: ChromosomeCounts 每个等位基因在每个染色体上的计数
 * @param _dp: Coverage 该位点的测序深度
 * @param _excess_het: ExcessHet 杂合子过量的统计检验
 * @param _fs: FisherStrand 用于检测参考和变异碱基在正负链上的分布是否均匀
 * @param _mq_rank_sum: InbreedingCoeff 群体内近亲交配系数
 * @param _qd: MappingQualityRankSumTest 用于检测参考和变异碱基的比对质量之间的差异
 * @param _mq: QualByDepth 每个样本的深度和质量的比率
 * @param _raw_mq_and_dp: RMSMappingQuality 该位点的比对质量的均方根
 * @param _read_pos_rank_sum: ReadPosRankSumTest 用于检测参考和变异碱基的读取位置之间的差异
 * @param _sor: StrandOddsRatio 用于检测参考和变异碱基在正负链上的分布是否均匀
 */
class InfoData
{
    /* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
private:
    Int32Vector _mleac;
    DoubleVector _mleaf;                        // %.2f         2
    double _base_qrank_sum{NEGATIVE_INFINITY};  // %.3f         1
    int32_t _an{INVALID_INT};
    Int32Vector _ac;
    DoubleVector _af;  // %.2f      2 标记为2的属性按照vcf_format
    int32_t _dp{INVALID_INT};
    double _excess_het{NEGATIVE_INFINITY};   // %.4f      1 标记为1的属性按照自身format
    double _fs{NEGATIVE_INFINITY};           // %.3f      1
    double _mq_rank_sum{NEGATIVE_INFINITY};  // %.3f      1
    double _qd{NEGATIVE_INFINITY};           // %.2f      1
    double _mq{NEGATIVE_INFINITY};           // %.2f      1
    LongVector _raw_mq_and_dp;
    double _read_pos_rank_sum{NEGATIVE_INFINITY};  // %.3f      1
    double _sor{NEGATIVE_INFINITY};                // %.3f      1

    std::bitset<INFO_FLAG_NUM> _exist_attributes{};

    pConstantsStr s_constant;

    /* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
public:
    static pInfoData create(pMemoryPool pool) { return new ALLOC_TYPE_IN_POOL(pool, InfoData) InfoData{pool}; }

    std::string& to_string(bool is_vcf_model, std::string& s) const { return is_vcf_model ? write_vcf_data(s) : write_gvcf_data(s); }

    void info2bcf(bcf_hdr_t* vcf_header, bcf1_t* b)
    {
        size_t af_size = _af.size();
        size_t mleaf_size = _mleaf.size();
        std::vector<float> convert(mleaf_size > af_size ? mleaf_size : af_size, 0.0f);
        float tmp;

        int32_t ret;
        if (_exist_attributes.test(IDEF_ChromosomeCounts)) {
            ret = bcf_update_info_int32(vcf_header, b, s_constant->k_allele_count_key.c_str(), _ac.data(), _ac.size());
            CHECK_CONDITION_EXIT(ret != 0, "write_chromosome_counts: ac");
            for (size_t i = 0; i < af_size; ++i) {
                convert[i] = static_cast<float>(_af[i]);
            }
            ret = bcf_update_info_float(vcf_header, b, s_constant->k_allele_frequency_key.c_str(), convert.data(), af_size);
            CHECK_CONDITION_EXIT(ret != 0, "write_chromosome_counts: af");
            ret = bcf_update_info_int32(vcf_header, b, s_constant->k_allele_number_key.c_str(), &_an, 1);
            CHECK_CONDITION_EXIT(ret != 0, "write_chromosome_counts: an");
        }
        if (_exist_attributes.test(IDEF_BaseQualityRankSumTest)) {
            tmp = static_cast<float>(_base_qrank_sum);
            ret = bcf_update_info_float(vcf_header, b, s_constant->k_base_qual_rank_sum_key.c_str(), &tmp, 1);
            CHECK_CONDITION_EXIT(ret != 0, "write_base_quality_rank_sum_test");
        }
        if (_exist_attributes.test(IDEF_Coverage)) {
            ret = bcf_update_info_int32(vcf_header, b, s_constant->k_depth_key.c_str(), &_dp, 1);
            CHECK_CONDITION_EXIT(ret != 0, "write_coverage");
        }
        if (_exist_attributes.test(IDEF_ExcessHet)) {
            tmp = static_cast<float>(_excess_het);
            ret = bcf_update_info_float(vcf_header, b, s_constant->k_excess_het_key.c_str(), &tmp, 1);
            CHECK_CONDITION_EXIT(ret != 0, "write_excesshet");
        }
        if (_exist_attributes.test(IDEF_FisherStrand)) {
            tmp = static_cast<float>(_fs);
            ret = bcf_update_info_float(vcf_header, b, s_constant->k_fisher_strand_key.c_str(), &tmp, 1);
            CHECK_CONDITION_EXIT(ret != 0, "write_fisher_strand");
        }
        if (_exist_attributes.test(IDEF_MLEAC)) {
            ret = bcf_update_info_int32(vcf_header, b, s_constant->k_mle_allele_count_key.c_str(), _mleac.data(), _mleac.size());
            CHECK_CONDITION_EXIT(ret != 0, "write_mleac");
        }
        if (_exist_attributes.test(IDEF_MLEAF)) {
            for (size_t i = 0; i < mleaf_size; ++i) {
                convert[i] = static_cast<float>(_mleaf[i]);
            }
            ret = bcf_update_info_float(vcf_header, b, s_constant->k_mle_allele_frequency_key.c_str(), convert.data(), mleaf_size);
            CHECK_CONDITION_EXIT(ret != 0, "write_mleaf");
        }
        if (_exist_attributes.test(IDEF_RMSMappingQuality_MQ)) {
            tmp = static_cast<float>(_mq);
            ret = bcf_update_info_float(vcf_header, b, s_constant->k_rms_mapping_quality_key.c_str(), &tmp, 1);
            CHECK_CONDITION_EXIT(ret != 0, "write_rms_mapping_quality_mq");
        }
        if (_exist_attributes.test(IDEF_InbreedingCoeff)) {
            // 这个tag是合并多个gvcf文件时使用的，hc中永远不会出现
            RovacaLogger::error("write_inbreeding_coeff: invalid code");
        }
        if (_exist_attributes.test(IDEF_MappingQualityRankSumTest)) {
            tmp = static_cast<float>(_mq_rank_sum);
            ret = bcf_update_info_float(vcf_header, b, s_constant->k_map_qual_rank_sum_key.c_str(), &tmp, 1);
            CHECK_CONDITION_EXIT(ret != 0, "write_mapping_quality_rank_sum_test");
        }
        if (_exist_attributes.test(IDEF_QualByDepth)) {
            tmp = static_cast<float>(_qd);
            ret = bcf_update_info_float(vcf_header, b, s_constant->k_qual_by_depth_key.c_str(), &tmp, 1);
            CHECK_CONDITION_EXIT(ret != 0, "write_qual_by_depth");
        }
        if (_exist_attributes.test(IDEF_RMSMappingQuality_RAW_MQandDP)) {
            int32_t tt[2]{};
            tt[0] = int32_t(_raw_mq_and_dp.at(0));
            tt[1] = int32_t(_raw_mq_and_dp.at(1));
            ret = bcf_update_info_int32(vcf_header, b, s_constant->k_raw_mapping_quality_with_depth_key.c_str(), tt, 2);
            CHECK_CONDITION_EXIT(ret != 0, "write_rms_mapping_quality_raw_mq_and_dp");
        }
        if (_exist_attributes.test(IDEF_ReadPosRankSumTest)) {
            tmp = static_cast<float>(_read_pos_rank_sum);
            ret = bcf_update_info_float(vcf_header, b, s_constant->k_read_pos_rank_sum_key.c_str(), &tmp, 1);
            CHECK_CONDITION_EXIT(ret != 0, "write_read_pos_rank_sum_test");
        }
        if (_exist_attributes.test(IDEF_StrandOddsRatio)) {
            tmp = static_cast<float>(_sor);
            ret = bcf_update_info_float(vcf_header, b, s_constant->k_strand_odds_ratio_key.c_str(), &tmp, 1);
            CHECK_CONDITION_EXIT(ret != 0, "write_strand_odds_ratio");
        }
    }

    void set_mleac(Int32Vector&& mleac)
    {
        _mleac = std::move(mleac);
        _exist_attributes.set(IDEF_MLEAC);
    }

    void set_mleac(const Int32Vector& mleac)
    {
        std::copy(mleac.begin(), mleac.end(), std::back_inserter(_mleac));
        _exist_attributes.set(IDEF_MLEAC);
    }

    void set_mleaf(DoubleVector&& mleaf)
    {
        _mleaf = std::move(mleaf);
        _exist_attributes.set(IDEF_MLEAF);
    }

    void set_base_qrank_sum(double base_qrank_sum)
    {
        _base_qrank_sum = base_qrank_sum;
        _exist_attributes.set(IDEF_BaseQualityRankSumTest);
    }

    void set_chromosome_counts(int32_t an, Int32Vector&& ac, DoubleVector&& af)
    {
        _an = an;
        _ac = std::move(ac);
        _af = std::move(af);
        _exist_attributes.set(IDEF_ChromosomeCounts);
    }

    void set_dp(int32_t dp)
    {
        _dp = dp;
        _exist_attributes.set(IDEF_Coverage);
    }

    void set_excess_het(double excess_het)
    {
        _excess_het = excess_het;
        _exist_attributes.set(IDEF_ExcessHet);
    }

    void set_fs(double fs)
    {
        _fs = fs;
        _exist_attributes.set(IDEF_FisherStrand);
    }

    void set_mq_rank_sum(double mq_rank_sum)
    {
        _mq_rank_sum = mq_rank_sum;
        _exist_attributes.set(IDEF_MappingQualityRankSumTest);
    }

    void set_qd(double qd)
    {
        _qd = qd;
        _exist_attributes.set(IDEF_QualByDepth);
    }

    void set_mq(double mq)
    {
        _mq = mq;
        _exist_attributes.set(IDEF_RMSMappingQuality_MQ);
    }

    void set_raw_mq_and_dp(LongVector&& raw_mq_and_dp)
    {
        _raw_mq_and_dp = std::move(raw_mq_and_dp);
        _exist_attributes.set(IDEF_RMSMappingQuality_RAW_MQandDP);
    }

    void set_read_pos_rank_sum(double read_pos_rank_sum)
    {
        _read_pos_rank_sum = read_pos_rank_sum;
        _exist_attributes.set(IDEF_ReadPosRankSumTest);
    }

    void set_sor(double sor)
    {
        _sor = sor;
        _exist_attributes.set(IDEF_StrandOddsRatio);
    }

    /* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
private:
    explicit InfoData(pMemoryPool pool)
        : _mleac(pool)
        , _mleaf(pool)
        , _ac(pool)
        , _af(pool)
        , _raw_mq_and_dp(pool)
        , s_constant(ConstantsStr::get_instance())
    {}

    std::string& write_vcf_data(std::string& s) const
    {
        write_chromosome_counts(s);
        write_base_quality_rank_sum_test(s);
        write_coverage(s);
        write_excesshet(s);
        write_fisher_strand(s);
        write_mleac(s);
        write_mleaf(s);
        write_rms_mapping_quality_mq(s);  // vcf 写出 MQ
        write_mapping_quality_rank_sum_test(s);
        write_qual_by_depth(s);
        write_read_pos_rank_sum_test(s);
        write_strand_odds_ratio(s);
        return s;
    }

    std::string& write_gvcf_data(std::string& s) const
    {
        write_base_quality_rank_sum_test(s);
        write_coverage(s);
        write_excesshet(s);
        write_mleac(s);
        write_mleaf(s);
        write_inbreeding_coeff(s);
        write_mapping_quality_rank_sum_test(s);
        write_rms_mapping_quality_raw_mq_and_dp(s);  // gvcf 写出 RAW_MQandDP
        write_read_pos_rank_sum_test(s);
        return s;
    }

    void write_mleac(std::string& s) const
    {
        if (_exist_attributes.test(IDEF_MLEAC)) {
            s.append(s_constant->k_mle_allele_count_key).append("=");
            DebugUtils::write_int_vector(_mleac, s).append(1, s_constant->k_info_field_separator);
        }
    }

    void write_mleaf(std::string& s) const
    {
        if (_exist_attributes.test(IDEF_MLEAF)) {
            s.append(s_constant->k_mle_allele_frequency_key).append("=");
            DebugUtils::write_double_vector(_mleaf, s).append(1, s_constant->k_info_field_separator);
        }
    }

    void write_base_quality_rank_sum_test(std::string& s) const
    {
        if (_exist_attributes.test(IDEF_BaseQualityRankSumTest)) {
            s.append(s_constant->k_base_qual_rank_sum_key).append("=");
            DebugUtils::format_double(_base_qrank_sum, "%.3f", s);
            s.append(1, s_constant->k_info_field_separator);
        }
    }

    void write_chromosome_counts(std::string& s) const
    {
        if (_exist_attributes.test(IDEF_ChromosomeCounts)) {
            s.append(s_constant->k_allele_count_key).append("=");
            DebugUtils::write_int_vector(_ac, s).append(1, s_constant->k_info_field_separator);
            s.append(s_constant->k_allele_frequency_key).append("=");
            DebugUtils::write_double_vector(_af, s).append(1, s_constant->k_info_field_separator);
            s.append(s_constant->k_allele_number_key).append("=").append(std::to_string(_an)).append(1, s_constant->k_info_field_separator);
        }
    }

    void write_coverage(std::string& s) const
    {
        if (_exist_attributes.test(IDEF_Coverage)) {
            s.append(s_constant->k_depth_key).append("=").append(std::to_string(_dp)).append(1, s_constant->k_info_field_separator);
        }
    }

    void write_excesshet(std::string& s) const
    {
        if (_exist_attributes.test(IDEF_ExcessHet)) {
            s.append(s_constant->k_excess_het_key).append("=");
            DebugUtils::format_double(_excess_het, "%.4f", s);
            s.append(1, s_constant->k_info_field_separator);
        }
    }

    void write_fisher_strand(std::string& s) const
    {
        if (_exist_attributes.test(IDEF_FisherStrand)) {
            s.append(s_constant->k_fisher_strand_key).append("=");
            DebugUtils::format_double(_fs, "%.3f", s);
            s.append(1, s_constant->k_info_field_separator);
        }
    }

    void write_inbreeding_coeff([[maybe_unused]] std::string& s) const
    {
        if (_exist_attributes.test(IDEF_InbreedingCoeff)) {
            RovacaLogger::error("invalid code");
        }
    }

    void write_mapping_quality_rank_sum_test(std::string& s) const
    {
        if (_exist_attributes.test(IDEF_MappingQualityRankSumTest)) {
            s.append(s_constant->k_map_qual_rank_sum_key).append("=");
            DebugUtils::format_double(_mq_rank_sum, "%.3f", s);
            s.append(1, s_constant->k_info_field_separator);
        }
    }

    void write_qual_by_depth(std::string& s) const
    {
        if (_exist_attributes.test(IDEF_QualByDepth)) {
            s.append(s_constant->k_qual_by_depth_key).append("=");
            DebugUtils::format_double(_qd, "%.2f", s);
            s.append(1, s_constant->k_info_field_separator);
        }
    }

    void write_rms_mapping_quality_mq(std::string& s) const
    {
        if (_exist_attributes.test(IDEF_RMSMappingQuality_MQ)) {
            s.append(s_constant->k_rms_mapping_quality_key).append("=");
            DebugUtils::format_double(_mq, "%.2f", s);
            s.append(1, s_constant->k_info_field_separator);
        }
    }

    void write_rms_mapping_quality_raw_mq_and_dp(std::string& s) const
    {
        if (_exist_attributes.test(IDEF_RMSMappingQuality_RAW_MQandDP)) {
            s.append(s_constant->k_raw_mapping_quality_with_depth_key).append("=");
            DebugUtils::write_long_vector(_raw_mq_and_dp, s).append(1, s_constant->k_info_field_separator);
        }
    }

    void write_read_pos_rank_sum_test(std::string& s) const
    {
        if (_exist_attributes.test(IDEF_ReadPosRankSumTest)) {
            s.append(s_constant->k_read_pos_rank_sum_key).append("=");
            DebugUtils::format_double(_read_pos_rank_sum, "%.3f", s);
            s.append(1, s_constant->k_info_field_separator);
        }
    }

    void write_strand_odds_ratio(std::string& s) const
    {
        if (_exist_attributes.test(IDEF_StrandOddsRatio)) {
            s.append(s_constant->k_strand_odds_ratio_key).append("=");
            DebugUtils::format_double(_sor, "%.3f", s);
            s.append(1, s_constant->k_info_field_separator);
        }
    }
};

}  // namespace rovaca

#endif  // ROVACA_HC_INFO_DATA_H_
