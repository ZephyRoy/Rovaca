#ifndef ROVACA_HC_ALLELE_LIKELIHOOD_MATRIX_MAPPER_H_
#define ROVACA_HC_ALLELE_LIKELIHOOD_MATRIX_MAPPER_H_
#include "forward.h"
#include "genotype_macors.h"
#include "interface/interface_likelihood_matrix.hpp"

namespace rovaca
{

template <typename EVIDENCE, typename A>
class PermutedAlleleLikelihood : public interfaceLikelihoodMatrix<EVIDENCE, A>
{
    typedef typename interfaceLikelihoodMatrix<EVIDENCE, A>::Ta Ta;
    typedef typename interfaceLikelihoodMatrix<EVIDENCE, A>::Te Te;
    typedef typename interfaceLikelihoodMatrix<EVIDENCE, A>::TaVector TaVector;
    typedef typename interfaceLikelihoodMatrix<EVIDENCE, A>::TeVector TeVector;

    typedef interfaceLikelihoodMatrix<EVIDENCE, A>* LikelihoodMatrixPtr;

    /* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
private:
    LikelihoodMatrixPtr _original;
    const InterfaceAlleleListPermutation<A>* _permutation;

    /* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
public:
    PermutedAlleleLikelihood(LikelihoodMatrixPtr original, const InterfaceAlleleListPermutation<A>* permutation)
        : _original(original)
        , _permutation(permutation)
    {}

    ~PermutedAlleleLikelihood() override = default;
    size_t number_of_alleles() const override { return _permutation->number_of_alleles(); }
    Ta get_allele(size_t index) const override { return _original->get_allele(_permutation->from_index(index)); }
    const TaVector& get_alleles() const override { return _permutation->get_alleles(); }
    int32_t index_of_allele(Ta allele) const override { return _permutation->index_of_allele(allele); }
    const TeVector& evidence() const override { return _original->evidence(); }
    size_t evidence_count() const override { return _original->evidence_count(); }
    int32_t index_of_evidence(Te e) const override { return _original->index_of_evidence(e); }
    Te get_evidence(size_t idx) const override { return _original->get_evidence(idx); }
    double get(size_t a_idx, size_t e_idx) const override { return _original->get(_permutation->from_index(a_idx), e_idx); }
    void set(size_t a_idx, size_t e_idx, double value) override { _original->set(_permutation->from_index(a_idx), e_idx, value); }
    void copy_allele_likelihoods(size_t allele_index, size_t offset, DoubleVector& dest) override
    {
        _original->copy_allele_likelihoods(_permutation->from_index(allele_index), offset, dest);
    }
};

template <typename A>
class AlleleLikelihoodMatrixMapper
{
    /* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
private:
    const InterfaceAlleleListPermutation<A>* _permutation;

    /* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
public:
    explicit AlleleLikelihoodMatrixMapper(const InterfaceAlleleListPermutation<A>* permutation)
        : _permutation(permutation)
    {}

    template <typename EVIDENCE>
    interfaceLikelihoodMatrix<EVIDENCE, A>* map_alleles(interfaceLikelihoodMatrix<EVIDENCE, A>* original, pMemoryPool pool)
    {
        if (_permutation->is_non_permuted()) {
            return original;
        }
        auto* p = pool->allocate(sizeof(PermutedAlleleLikelihood<EVIDENCE, A>));
        return new (p) PermutedAlleleLikelihood<EVIDENCE, A>(original, _permutation);
    }
};

}  // namespace rovaca

#endif  // ROVACA_HC_ALLELE_LIKELIHOOD_MATRIX_MAPPER_H_
