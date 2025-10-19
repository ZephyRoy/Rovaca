#ifndef ROVACA_HC_GENOTYPES_CONTEXT_H_
#define ROVACA_HC_GENOTYPES_CONTEXT_H_
#include <vector>

#include "allele.h"
#include "rovaca_logger.h"
#include "forward.h"
#include "genotype.h"
#include "genotype_macors.h"

namespace rovaca
{

class GenotypesContext
{
    static constexpr size_t s_default_size = 2;

    /* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
private:
    std::pmr::vector<pGenotype> _gs;
    int32_t _max_ploidy{INVALID_INT};

    /* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
public:
    DISALLOW_COPY_AND_ASSIGN(GenotypesContext);
    ~GenotypesContext() = default;

    static pGenotypesContext create(pMemoryPool pool) { return new ALLOC_TYPE_IN_POOL(pool, GenotypesContext) GenotypesContext(pool); }

    void add(pGenotype g) { _gs.push_back(g); }
    size_t size() const { return _gs.size(); }
    bool empty() const { return _gs.empty(); }
    pGenotype at(size_t i) { return _gs.at(i); }

    int32_t get_max_ploidy(int32_t default_ploidy)
    {
        CHECK_CONDITION_EXIT(default_ploidy < 0, "default_ploidy must be greater than or equal to 0");
        if (INVALID_INT == _max_ploidy) {
            _max_ploidy = 0;
            for (pGenotype g : _gs) {
                _max_ploidy = std::max(_max_ploidy, g->get_ploidy());
            }
            if (0 == _max_ploidy) {
                _max_ploidy = default_ploidy;
            }
        }
        return _max_ploidy;
    }

    int32_t allele_num() const
    {
        int32_t count = 0;
        for (pGenotype g : _gs) {
            for (pAllele a : g->alleles()) {
                count += a->is_called() ? 1 : 0;
            }
        }
        return count;
    }

    /* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
private:
    explicit GenotypesContext(pMemoryPool pool)
        : _gs(pool)
    {
        _gs.reserve(s_default_size);
    }
};

}  // namespace rovaca

#endif  // ROVACA_HC_GENOTYPES_CONTEXT_H_
