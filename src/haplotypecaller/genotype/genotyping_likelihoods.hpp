#ifndef ROVACA_HC_GENOTYPING_LIKELIHOODS_H_
#define ROVACA_HC_GENOTYPING_LIKELIHOODS_H_
#include "forward.h"
#include "genotype_macors.h"
#include "interface/interface_allele_list.hpp"
#include "interface/interface_ploidy_model.hpp"
#include "interface/interface_sample_list.hpp"

namespace rovaca
{

template <typename A>
class GenotypingLikelihoods : public InterfaceSampleList, public InterfaceAlleleList<A>
{
    /* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
public:
    typedef InterfaceAlleleList<A>* pInterfaceAlleleList;

    typedef typename InterfaceAlleleList<A>::Ta Ta;
    typedef typename InterfaceAlleleList<A>::TaVector TaVector;

    /* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
private:
    pInterfaceAlleleList _alleles;
    pInterfacePloidyModel _ploidy_model;
    std::pmr::vector<pGenotypeLikelihoods> _likelihoods;

    /* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
public:
    static GenotypingLikelihoods* create(pInterfaceAlleleList allele, pInterfacePloidyModel ploidy_model,
                                         std::pmr::vector<pGenotypeLikelihoods>&& likelihoods)
    {
        pMemoryPool pool = likelihoods.get_allocator().resource();
        return new ALLOC_TYPE_IN_POOL(pool, GenotypingLikelihoods<A>)
            GenotypingLikelihoods<A>{allele, ploidy_model, std::forward<std::pmr::vector<pGenotypeLikelihoods>&&>(likelihoods)};
    }

    ~GenotypingLikelihoods() override = default;
    size_t number_of_alleles() const override { return _alleles->number_of_alleles(); }
    Ta get_allele(size_t index) const override { return _alleles->get_allele(index); }
    const TaVector& get_alleles() const override { return _alleles->get_alleles(); }
    int32_t index_of_allele(Ta allele) const override { return _alleles->index_of_allele(allele); }
    size_t number_of_samples() const override { return _ploidy_model->number_of_samples(); }
    const std::string& get_sample(size_t index) const override { return _ploidy_model->get_sample(index); }
    int32_t index_of_sample(const std::string& name) const override { return _ploidy_model->index_of_sample(name); }

    int32_t sample_ploidy(size_t index) const { return _ploidy_model->sample_ploidy(index); }

    pGenotypeLikelihoods sample_likelihoods(int32_t sample_index) const { return _likelihoods.at(sample_index); }

    /* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
private:
    GenotypingLikelihoods(pInterfaceAlleleList allele, pInterfacePloidyModel ploidy_model,
                          std::pmr::vector<pGenotypeLikelihoods>&& likelihoods)
        : _alleles(allele)
        , _ploidy_model(ploidy_model)
        , _likelihoods(std::move(likelihoods))
    {}
};

}  // namespace rovaca

#endif  // ROVACA_HC_GENOTYPING_LIKELIHOODS_H_
