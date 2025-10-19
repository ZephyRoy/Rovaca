#include "strand_bias_test.h"

#include "allele_likelihoods.hpp"
#include "rovaca_logger.h"
#include "genotype.h"
#include "genotype_macors.h"
#include "genotypes_context.hpp"
#include "read_record.h"
#include "variant.h"

namespace rovaca
{

static constexpr int32_t s_array_dim = 2;

void StrandBiasTest::annotate([[maybe_unused]] pRefFragment ref, pVariant vc, pRALikelihoods likelihoods, pInfoData target,
                              pMemoryPool pool)
{
    if (!vc->is_variant()) {
        return;
    }

    if (vc->has_genotypes()) {
        pGenotypesContext gc = vc->genotype();
        for (size_t i = 0, len = gc->size(); i < len; ++i) {
            if (gc->at(i)->has_sb()) {
                calculate_annotation_from_gt_field(gc, target, pool);
                return;
            }
        }
    }

    if (nullptr != likelihoods) {
        calculate_annotation_from_likelihoods(vc, likelihoods, target, pool);
    }
}

Int32Vector2D StrandBiasTest::get_contingency_table(pVariant vc, pRALikelihoods likelihoods, int32_t min_count, pMemoryPool pool)
{
    pAllele ref = vc->ref_allele();
    const AlleleVector& alt_alleles = vc->alt_alleles();

    return get_contingency_table(likelihoods, ref, alt_alleles, min_count, pool);
}

Int32Vector2D StrandBiasTest::get_contingency_table(pRALikelihoods likelihoods, pAllele ref, const AlleleVector& alt_alleles,
                                                    int32_t min_count, pMemoryPool pool)
{
    size_t sample_count = likelihoods->number_of_samples();
    return get_contingency_table(likelihoods, ref, alt_alleles, min_count, sample_count, pool);
}

bool StrandBiasTest::passes_minimum_threshold(const Int32Vector& data, int32_t min_count)
{
    // the ref and alt totals must be greater than min_count
    return data[0] + data[1] + data[2] + data[3] > min_count;
}

Int32Vector2D StrandBiasTest::get_contingency_table(pRALikelihoods likelihoods, pAllele ref, const AlleleVector& alt_alleles,
                                                    int32_t min_count, size_t sample_count, pMemoryPool pool)
{
    AlleleSet all_alts(alt_alleles.begin(), alt_alleles.end(), pool);
    Int32Vector2D table(s_array_dim, Int32Vector(s_array_dim), pool);
    Int32Vector sample_table(s_array_dim * 2, pool);
    for (int32_t si = 0; si < (int32_t)sample_count; ++si) {
        sample_table.clear();
        sample_table.resize(s_array_dim * 2);
        auto best_alleles = likelihoods->best_alleles_breaking_ties(si);

        std::for_each(best_alleles.begin(), best_alleles.end(), [&](RALikelihoods::BestAllele* b) {
            if (b->is_informative()) {
                update_table(sample_table, b->best_allele, b->evidence, ref, all_alts);
            }
        });

        if (passes_minimum_threshold(sample_table, min_count)) {
            copy_to_main_table(sample_table, table);
        }
    }

    return table;
}

void StrandBiasTest::update_table(Int32Vector& table, pAllele a, pReadRecord read, pAllele ref, const AlleleSet& all_alts)
{
    bool matches_ref = a->equals(*ref, true);
    bool matches_any_alt = all_alts.count(a);
    if (matches_ref || matches_any_alt) {
        int32_t offset = matches_ref ? 0 : s_array_dim;
        bool is_forward = !read->is_reverse_strand();  // a normal read with an actual strand
        table[offset + (is_forward ? 0 : 1)]++;
    }
}

void StrandBiasTest::copy_to_main_table(const Int32Vector& per_sample_table, Int32Vector2D& main_table)
{
    main_table[0][0] += per_sample_table[0];
    main_table[0][1] += per_sample_table[1];
    main_table[1][0] += per_sample_table[2];
    main_table[1][1] += per_sample_table[3];
}

}  // namespace rovaca