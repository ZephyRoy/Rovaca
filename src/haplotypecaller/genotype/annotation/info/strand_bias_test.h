#ifndef ROVACA_HC_STRAND_BIAS_TEST_H_
#define ROVACA_HC_STRAND_BIAS_TEST_H_
#include "interface/interface_info_field_annotation.hpp"

namespace rovaca
{

class StrandBiasTest : public InterfaceInfoFieldAnnotation
{
    /* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
public:
    void annotate(pRefFragment ref, pVariant vc, pRALikelihoods likelihoods, pInfoData target, pMemoryPool pool) override;

    /**
     * allocate and fill a 2x2 strand contingency table.  in the end, it'll look something like this:
     *             fw      rc
     *   allele1   #       #
     *   allele2   #       #
     * @return a 2x2 contingency table
     */
    static Int32Vector2D get_contingency_table(pVariant vc, pRALikelihoods likelihoods, int32_t min_count, pMemoryPool pool);
    static Int32Vector2D get_contingency_table(pRALikelihoods likelihoods, pAllele ref, const AlleleVector& alt_alleles, int32_t min_count,
                                               pMemoryPool pool);

    ~StrandBiasTest() override = default;

    /* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
protected:
    virtual void calculate_annotation_from_gt_field(pGenotypesContext gc, pInfoData target, pMemoryPool pool) = 0;

    virtual void calculate_annotation_from_likelihoods(pVariant vc, pRALikelihoods likelihoods, pInfoData target, pMemoryPool pool) = 0;

    /**
     * does this strand data array pass the minimum threshold for inclusion?
     * @param data  the array
     * @param min_count the minimum threshold of counts in the array
     * @return true if it passes the minimum threshold, false otherwise
     */
    static bool passes_minimum_threshold(const Int32Vector& data, int32_t min_count);

    /* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
private:
    static Int32Vector2D get_contingency_table(pRALikelihoods likelihoods, pAllele ref, const AlleleVector& alt_alleles, int32_t min_count,
                                               size_t sample_count, pMemoryPool pool);

    static void update_table(Int32Vector& table, pAllele a, pReadRecord read, pAllele ref, const AlleleSet& all_alts);

    /**
     * helper method to copy the per-sample table to the main table
     * @param per_sample_table   per-sample table (single dimension)
     * @param main_table        main table (two dimensions)
     */
    static void copy_to_main_table(const Int32Vector& per_sample_table, Int32Vector2D& main_table);
};

}  // namespace rovaca

#endif  // ROVACA_HC_STRAND_BIAS_TEST_H_
