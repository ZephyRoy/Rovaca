#include "variant_annotator_engine.h"

#include "annotation/format/depth_per_allele_by_sample.h"
#include "annotation/format/depth_per_sample_hc.h"
#include "annotation/format/strand_bias_by_sample.h"
#include "annotation/info/base_quality_rank_sum_test.h"
#include "annotation/info/chromosome_counts.h"
#include "annotation/info/coverage.h"
#include "annotation/info/excess_het.h"
#include "annotation/info/fisher_strand.h"
#include "annotation/info/inbreeding_coeff.h"
#include "annotation/info/mapping_quality_rank_sum_test.h"
#include "annotation/info/qual_by_depth.h"
#include "annotation/info/read_pos_rank_sum_test.h"
#include "annotation/info/rms_mapping_quality.h"
#include "annotation/info/strand_odds_ratio.h"
#include "genotype.h"
#include "genotypes_context.hpp"
#include "rovaca_logger.h"
#include "utils/rovaca_variant_context_utils.h"
#include "variant.h"

namespace rovaca
{

VariantAnnotatorEngine::VariantAnnotatorEngine(bool use_raw)
{
    if (use_raw) {
        gvcf_init_info_annotations();
        gvcf_init_genotype_annotations();
    }
    else {
        vcf_init_info_annotations();
        vcf_init_genotype_annotations();
    }
}

VariantAnnotatorEngine::~VariantAnnotatorEngine()
{
    for (auto* p : _info_annotations) {
        delete p;
    }

    for (auto* p : _genotype_annotations) {
        delete p;
    }
}

void VariantAnnotatorEngine::gvcf_init_info_annotations()
{
    _info_annotations.push_back(new BaseQualityRankSumTest{});
    _info_annotations.push_back(new Coverage{});
    _info_annotations.push_back(new ExcessHet{});
    _info_annotations.push_back(new InbreedingCoeff{});
    _info_annotations.push_back(new MappingQualityRankSumTest{});
    _info_annotations.push_back(new RMSMappingQuality{true});
    _info_annotations.push_back(new ReadPosRankSumTest{});
}

void VariantAnnotatorEngine::gvcf_init_genotype_annotations()
{
    _genotype_annotations.push_back(new DepthPerAlleleBySample{});
    _genotype_annotations.push_back(new DepthPerSampleHC{});
    _genotype_annotations.push_back(new StrandBiasBySample{});
}

void VariantAnnotatorEngine::vcf_init_info_annotations()
{
    _info_annotations.push_back(new BaseQualityRankSumTest{});
    _info_annotations.push_back(new ChromosomeCounts{});
    _info_annotations.push_back(new Coverage{});
    _info_annotations.push_back(new ExcessHet{});
    _info_annotations.push_back(new FisherStrand{});
    _info_annotations.push_back(new InbreedingCoeff{});
    _info_annotations.push_back(new MappingQualityRankSumTest{});
    _info_annotations.push_back(new QualByDepth{});
    _info_annotations.push_back(new RMSMappingQuality{false});
    _info_annotations.push_back(new ReadPosRankSumTest{});
    _info_annotations.push_back(new StrandOddsRatio{});
}

void VariantAnnotatorEngine::vcf_init_genotype_annotations()
{
    _genotype_annotations.push_back(new DepthPerAlleleBySample{});
    _genotype_annotations.push_back(new DepthPerSampleHC{});
}

pVariant VariantAnnotatorEngine::annotate_context(pRefFragment ref, pVariant vc, pRALikelihoods likelihoods, int32_t db_offset,
                                                  std::vector<bcf1_t*>* db_data, pMemoryPool pool)
{
    CHECK_CONDITION_EXIT(nullptr == vc, "nullptr == vc");
    CHECK_CONDITION_EXIT(nullptr == ref, "nullptr == ref");
    CHECK_CONDITION_EXIT(nullptr == pool, "nullptr == pool");
    CHECK_CONDITION_EXIT(nullptr == likelihoods, "nullptr == likelihoods");

    annotate_genotypes(ref, vc, likelihoods, pool);
    annotate_info(ref, vc, likelihoods, pool);

    if (db_data) {
        std::pmr::vector<bcf1_t*> match_bcf{pool};
        int32_t data_len = static_cast<int32_t>(db_data->size());
        int32_t vc_start = static_cast<int32_t>(vc->get_start());
        int32_t pos_1base;
        for (; db_offset < data_len; ++db_offset) {
            pos_1base = db_data->at(db_offset)->pos + 1;
            if (pos_1base > vc_start) {
                break;
            }
            if (pos_1base == vc_start) {
                match_bcf.push_back(db_data->at(db_offset));
            }
        }

        if (!match_bcf.empty()) {
            uint32_t max_possibility_rsid_len = 0;
            std::for_each(match_bcf.begin(), match_bcf.end(), [&](bcf1_t* b) {
                bcf_unpack(b, BCF_UN_STR);
                max_possibility_rsid_len += strlen(b->d.id);
            });

            pBases id = new ALLOC_FLEXIBLE_IN_POOL(pool, Bases, max_possibility_rsid_len, uint8_t) Bases{0};
            id->num = 0;
            id->data[0] = 0;

            std::pmr::vector<Event> vc_annotate_list = ROVACAVariantContextUtils::split_variant_context_to_biallelics_event(vc, true, pool);

            for (bcf1_t* b : match_bcf) {
                if (b->d.n_flt) {
                    continue;
                }

                pVariant rs_id_source_vc = Variant::create_dbsnp_source(b, pool);
                CHECK_CONDITION_EXIT(rs_id_source_vc->get_tid() != vc->get_tid(), "dbsnp tid not equal at: rs_start=%ld, vc_start=%ld.",
                                     rs_id_source_vc->get_tid(), vc->get_tid());

                bool add_this_id = false;
                std::pmr::vector<Event> vc_comp_list =
                    ROVACAVariantContextUtils::split_variant_context_to_biallelics_event(rs_id_source_vc, true, pool);

                for (const Event& vc_comp : vc_comp_list) {
                    for (const Event& vc_annotate : vc_annotate_list) {
                        if (vc_comp == vc_annotate) {
                            add_this_id = true;
                            break;
                        }
                    }

                    if (add_this_id) {
                        pBases current_id = rs_id_source_vc->db_id();
                        if (0 == id->num) {
                            memcpy(id->data, current_id->data, current_id->num * sizeof(uint8_t));
                            id->num = current_id->num;
                        }
                        else {
                            id->data[id->num] = ';';
                            memcpy(id->data + id->num + 1, current_id->data, current_id->num * sizeof(uint8_t));
                            id->num += current_id->num + 1;
                        }
                        id->data[id->num] = 0;
                        break;
                    }
                }
            }

            if (0 != id->num) {
                vc->set_id(id);
            }
        }
    }

    return vc;
}

void VariantAnnotatorEngine::annotate_genotypes(pRefFragment ref, pVariant vc, pRALikelihoods likelihoods, pMemoryPool pool)
{
    pGenotypesContext genotypes = vc->genotype();
    pGenotype g;
    for (size_t i = 0, len = genotypes->size(); i < len; ++i) {
        g = genotypes->at(i);
        for (pInterfaceGenotypeAnnotation a : _genotype_annotations) {
            a->annotate(ref, vc, g, likelihoods, pool);
        }
    }
}

void VariantAnnotatorEngine::annotate_info(pRefFragment ref, pVariant vc, pRALikelihoods likelihoods, pMemoryPool pool)
{
    pInfoData info = vc->info();
    for (pInterfaceInfoFieldAnnotation a : _info_annotations) {
        a->annotate(ref, vc, likelihoods, info, pool);
    }
}

}  // namespace rovaca