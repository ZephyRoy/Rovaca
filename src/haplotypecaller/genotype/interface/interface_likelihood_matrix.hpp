#ifndef ROVACA_HC_INTERFACE_LIKELIHOOD_MATRIX_H_
#define ROVACA_HC_INTERFACE_LIKELIHOOD_MATRIX_H_
#include <vector>

#include "forward.h"
#include "interface_allele_list.hpp"

namespace rovaca
{

template <typename EVIDENCE, typename A>
class interfaceLikelihoodMatrix : public InterfaceAlleleList<A>
{
    /* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
public:
    typedef EVIDENCE Te;
    typedef typename InterfaceAlleleList<A>::Ta Ta;
    typedef std::pmr::vector<Te> TeVector;
    typedef typename InterfaceAlleleList<A>::TaVector TaVector;

    /* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
public:
    virtual const TeVector& evidence() const = 0;
    virtual size_t evidence_count() const = 0;
    virtual int32_t index_of_evidence(Te e) const = 0;
    virtual Te get_evidence(size_t idx) const = 0;

    /*! @brief returns the likelihood of a unit of evidence given a haplotype */
    virtual double get(size_t a_idx, size_t e_idx) const = 0;

    /*! @brief set the likelihood of a unit of evidence given an allele through their indices */
    virtual void set(size_t a_idx, size_t e_idx, double value) = 0;

    /*! @brief copies the likelihood of all the evidence for a given allele into an array from a particular offset */
    virtual void copy_allele_likelihoods(size_t allele_index, size_t offset, DoubleVector& dest) = 0;
};

}  // namespace rovaca

#endif  // ROVACA_HC_INTERFACE_LIKELIHOOD_MATRIX_H_
