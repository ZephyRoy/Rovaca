#include "genotype.h"

#include "allele.h"
#include "constants_str.hpp"
#include "rovaca_logger.h"
#include "genotype_likelihoods.h"
#include "utils/debug_utils.h"

namespace rovaca
{

pGenotypeLikelihoods Genotype::get_likelihoods() const { return GenotypeLikelihoods::create(_pl, _pl.get_allocator().resource()); }

GenotypeType Genotype::get_type()
{
    if (_type == GT_UNINITIALIZED) {
        _type = determine_type();
    }
    return static_cast<GenotypeType>(_type);
}

int32_t Genotype::count_allele(pAllele a) const
{
    int32_t count = 0;
    for (pAllele allele : _alleles) {
        if (allele->equals(*a)) {
            count += 1;
        }
    }
    return count;
}

bool Genotype::is_het_non_ref() { return is_het() && !_alleles.at(0)->is_reference() && !_alleles.at(1)->is_reference(); }

std::string& Genotype::append_genotype_data(std::string& s, const std::pmr::map<pAllele, char>& mm)
{
    pConstantsStr con = ConstantsStr::get_instance();

    bool saw_good_gt = is_available();
    bool saw_good_qual = has_gq();
    bool saw_dp = has_dp();
    bool saw_ad = has_ad();
    bool saw_pl = has_likelihoods();

    if (saw_good_gt) s.append(con->k_genotype_key).append(1, con->k_genotype_field_separator);
    if (saw_ad) s.append(con->k_genotype_allele_depths).append(1, con->k_genotype_field_separator);
    if (saw_dp) s.append(con->k_depth_key).append(1, con->k_genotype_field_separator);
    if (saw_good_qual) s.append(con->k_genotype_quality_key).append(1, con->k_genotype_field_separator);
    if (_phased)
        s.append(con->k_haplotype_caller_phasing_gt_key)
            .append(1, con->k_genotype_field_separator)
            .append(con->k_haplotype_caller_phasing_id_key)
            .append(1, con->k_genotype_field_separator);
    if (saw_pl) s.append(con->k_genotype_pl_key).append(1, con->k_genotype_field_separator);
    if (_phased) s.append(con->k_phase_set_key).append(1, con->k_genotype_field_separator);
    if (!_sb.empty()) s.append(con->k_strand_bias_key).append(1, con->k_genotype_field_separator);

    // 最后一个 : 替换为 \t
    s.back() = con->k_field_separator;

    // "GT:AD:DP:GQ:PGT:PID:PL:PS:SB"
    if (saw_good_gt) {
        s.append(1, mm.at(_alleles.at(0)));
        for (int32_t i = 1, ploidy = get_ploidy(); i < ploidy; ++i) {
            s.append(1, is_phased() ? con->k_phased : con->k_unphased);
            s.append(1, mm.at(_alleles.at(i)));
        }
        s.append(1, con->k_genotype_field_separator);
    }
    if (saw_ad) {
        DebugUtils::write_int_vector(_ad, s).append(1, con->k_genotype_field_separator);
    }
    if (saw_dp) {
        s.append(std::to_string(_dp)).append(1, con->k_genotype_field_separator);
    }
    if (saw_good_qual) {
        s.append(std::to_string(std::min(MAX_GENOTYPE_QUAL, _gq))).append(1, con->k_genotype_field_separator);
    }
    if (_phased) {
        s.append(_phased_data->pgt)
            .append(1, con->k_genotype_field_separator)
            .append((char*)_phased_data->pid->data, _phased_data->pid->num)
            .append(1, con->k_genotype_field_separator);
    }
    if (saw_pl) {
        DebugUtils::write_int_vector(_pl, s).append(1, con->k_genotype_field_separator);
    }
    if (_phased) {
        s.append(std::to_string(_phased_data->ps)).append(1, con->k_genotype_field_separator);
    }
    if (!_sb.empty()) {
        DebugUtils::write_int_vector(_sb, s).append(1, con->k_genotype_field_separator);
    }

    return s;
}

void Genotype::genotype2bcf(bcf_hdr_t* vcf_header, bcf1_t* rec, const std::pmr::map<pAllele, int32_t>& mm)
{
    pConstantsStr k_constan = ConstantsStr::get_instance();

    bool saw_dp = false;
    bool saw_ad = false;
    bool saw_pl = false;
    bool saw_good_gt = false;
    bool saw_good_qual = false;

    if (has_dp()) saw_dp = true;
    if (has_ad()) saw_ad = true;
    if (has_gq()) saw_good_qual = true;
    if (has_likelihoods()) saw_pl = true;
    if (is_available()) saw_good_gt = true;

    int32_t ret;
    // "GT:AD:DP:GQ:PGT:PID:PL:PS:SB"
    if (saw_good_gt) {
        int32_t gt[HOM_GT_CUNT];
        if (ROVACA_UNLIKELY(_alleles.at(0) == _alleles.at(1) && _alleles.at(0) == StaticAllele::get_instance()->_no_call.get())) {
            gt[0] = gt[1] = _phased ? 1 : 0;
        }
        else {
            int32_t idx0 = mm.at(_alleles.at(0));
            int32_t idx1 = mm.at(_alleles.at(1));
            gt[0] = _phased ? bcf_gt_phased(idx0) : bcf_gt_unphased(idx0);
            gt[1] = _phased ? bcf_gt_phased(idx1) : bcf_gt_unphased(idx1);
        }
        ret = bcf_update_genotypes(vcf_header, rec, gt, bcf_hdr_nsamples(vcf_header) * 2);
        CHECK_CONDITION_EXIT(ret != 0, "bcf_update_genotypes");
    }
    if (saw_ad) {
        ret = bcf_update_format_int32(vcf_header, rec, k_constan->k_genotype_allele_depths.c_str(), _ad.data(), _ad.size());
        CHECK_CONDITION_EXIT(ret != 0, "bcf_update_info_int32: {}", k_constan->k_genotype_allele_depths.c_str());
    }
    if (saw_dp) {
        ret = bcf_update_format_int32(vcf_header, rec, k_constan->k_depth_key.c_str(), &_dp, 1);
        CHECK_CONDITION_EXIT(ret != 0, "bcf_update_info_int32: {}", k_constan->k_depth_key.c_str());
    }
    if (saw_good_qual) {
        int32_t filter_gq = std::min(_gq, MAX_GENOTYPE_QUAL);
        ret = bcf_update_format_int32(vcf_header, rec, k_constan->k_genotype_quality_key.c_str(), &filter_gq, 1);
        CHECK_CONDITION_EXIT(ret != 0, "bcf_update_info_int32: {}", k_constan->k_genotype_quality_key.c_str());
    }
    if (_phased) {
        ret = bcf_update_format_string(vcf_header, rec, k_constan->k_haplotype_caller_phasing_gt_key.c_str(), &_phased_data->pgt, 1);
        CHECK_CONDITION_EXIT(ret != 0, "bcf_update_info_string: {}", k_constan->k_haplotype_caller_phasing_gt_key.c_str());

        const char* tmp = (const char*)_phased_data->pid->data;
        ret = bcf_update_format_string(vcf_header, rec, k_constan->k_haplotype_caller_phasing_id_key.c_str(), &tmp, 1);
        CHECK_CONDITION_EXIT(ret != 0, "bcf_update_info_string: {}", k_constan->k_haplotype_caller_phasing_id_key.c_str());
    }
    if (saw_pl) {
        ret = bcf_update_format_int32(vcf_header, rec, k_constan->k_genotype_pl_key.c_str(), _pl.data(), _pl.size());
        CHECK_CONDITION_EXIT(ret != 0, "bcf_update_info_int32: {}", k_constan->k_genotype_pl_key.c_str());
    }
    if (_phased) {
        ret = bcf_update_format_int32(vcf_header, rec, k_constan->k_phase_set_key.c_str(), &_phased_data->ps, 1);
        CHECK_CONDITION_EXIT(ret != 0, "bcf_update_info_int32: {}", k_constan->k_phase_set_key.c_str());
    }
    if (!_sb.empty()) {
        ret = bcf_update_format_int32(vcf_header, rec, k_constan->k_strand_bias_key.c_str(), _sb.data(), _sb.size());
        CHECK_CONDITION_EXIT(ret != 0, "bcf_update_info_int32: {}", k_constan->k_strand_bias_key.c_str());
    }
}

Genotype::Genotype(pMemoryPool pool)
    : _alleles(pool)
    , _ad(pool)
    , _pl(pool)
    , _sb(pool)
    , _phased(0)
    , _type(GT_UNINITIALIZED)
    , _sample_id(INVALID_INT)
{}

GenotypeType Genotype::determine_type()
{
    if (_alleles.empty()) {
        return GT_UNAVAILABLE;
    }
    bool saw_no_call = false, saw_multiple_alleles = false;
    pAllele first_call_allele = nullptr;

    for (pAllele a : _alleles) {
        if (!a->is_called()) {
            saw_no_call = true;
        }
        else if (nullptr == first_call_allele) {
            first_call_allele = a;
        }
        else if (!a->equals(*first_call_allele)) {
            saw_multiple_alleles = true;
        }
    }

    if (saw_no_call) {
        if (first_call_allele == nullptr) {
            return GT_NO_CALL;
        }
        return GT_MIXED;
    }

    CHECK_CONDITION_EXIT(first_call_allele == nullptr,
                         "bug: there are no alleles present in this genotype but the alleles list is not nullptr");

    return saw_multiple_alleles ? GT_HET : first_call_allele->is_reference() ? GT_HOM_REF : GT_HOM_VAR;
}

}  // namespace rovaca