#ifndef ROVACA_HC_VARIANT_ANNOTATOR_ENGINE_H_
#define ROVACA_HC_VARIANT_ANNOTATOR_ENGINE_H_
#include "forward.h"
#include "genotype_macors.h"
#include "htslib/vcf.h"

namespace rovaca
{

class VariantAnnotatorEngine
{
    /* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
private:
    std::vector<pInterfaceInfoFieldAnnotation> _info_annotations;
    std::vector<pInterfaceGenotypeAnnotation> _genotype_annotations;

    /* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
public:
    DISALLOW_COPY_AND_ASSIGN(VariantAnnotatorEngine);
    explicit VariantAnnotatorEngine(bool use_raw);
    ~VariantAnnotatorEngine();

    pVariant annotate_context(pRefFragment ref, pVariant vc, pRALikelihoods likelihoods, int32_t db_offset, std::vector<bcf1_t*>* db_data,
                              pMemoryPool pool);

    /* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
private:
    void gvcf_init_info_annotations();
    void gvcf_init_genotype_annotations();

    void vcf_init_info_annotations();
    void vcf_init_genotype_annotations();

    void annotate_genotypes(pRefFragment ref, pVariant vc, pRALikelihoods likelihoods, pMemoryPool pool);
    void annotate_info(pRefFragment ref, pVariant vc, pRALikelihoods likelihoods, pMemoryPool pool);
};

}  // namespace rovaca

#endif  // ROVACA_HC_VARIANT_ANNOTATOR_ENGINE_H_
