#include "strand_bias_by_sample.h"

#include "annotation/info/fisher_strand.h"
#include "rovaca_logger.h"
#include "genotype.h"

namespace rovaca
{

static size_t s_array_dim = 2;

void StrandBiasBySample::annotate([[maybe_unused]] pRefFragment ref, pVariant vc, pGenotype g, pRALikelihoods likelihoods, pMemoryPool pool)
{
    if (g->has_sb() && nullptr == likelihoods) {
        return;
    }
    if (nullptr == likelihoods || !g->is_called()) {
        return;
    }

    Int32Vector2D table = FisherStrand::get_contingency_table(vc, likelihoods, 0, pool);
    Int32Vector sb = get_contingency_array(table, pool);
    g->set_sb(std::move(sb));
}

Int32Vector StrandBiasBySample::get_contingency_array(const Int32Vector2D& table, pMemoryPool pool)
{
    CHECK_CONDITION_EXIT(table.size() != s_array_dim || table.at(0).size() != s_array_dim, "must equal 2");
    Int32Vector result(pool);
    result.reserve(s_array_dim * 2);
    for (size_t i = 0; i < s_array_dim; ++i) {
        for (size_t j = 0; j < s_array_dim; ++j) {
            result.push_back(table[i][j]);
        }
    }
    return result;
}

}  // namespace rovaca