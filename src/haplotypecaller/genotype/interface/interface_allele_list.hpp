#ifndef ROVACA_HC_INTERFACE_ALLELE_LIST_H_
#define ROVACA_HC_INTERFACE_ALLELE_LIST_H_
#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <memory_resource>
#include <vector>

#include "rovaca_logger.h"
#include "genotype_macors.h"
#include "interface_permutation.hpp"

namespace rovaca
{

typedef std::pmr::memory_resource MemoryPool, *pMemoryPool;

template <typename A>
class InterfaceAlleleListPermutation;

template <typename A>
class InterfaceAlleleList
{
    class NonPermutation;
    class ActualPermutation;

    /* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
public:
    typedef A Ta;
    typedef std::pmr::vector<Ta> TaVector;

    /* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
public:
    virtual ~InterfaceAlleleList() = default;

    /*!
     * @brief returns the number of alleles in this AlleleList
     * @return
     */
    virtual size_t number_of_alleles() const = 0;

    /*!
     * @brief returns the allele at the given index in this AlleleList
     * @param index
     * @return
     */
    virtual Ta get_allele(size_t index) const = 0;

    /*!
     * @brief returns all allele
     * @return
     */
    virtual const TaVector& get_alleles() const = 0;

    /*!
     * @brief returns the index of the given Allele in this AlleleList
     * @param allele
     * @return -1 if the given allele is not present in this AlleleList
     */
    virtual int32_t index_of_allele(Ta allele) const = 0;

    /*!
     * @brief
     * @param allele
     * @return true if this AlleleList contains the specified allele
     */
    bool contains_allele(Ta allele) const { return index_of_allele(allele) >= 0; }

    /*!
     * @brief Resolves the index of the reference allele in an AlleleList
     * @return -1 if there is no reference allele
     */
    int32_t index_of_reference() const
    {
        for (int32_t i = 0, allele_count = (int32_t)number_of_alleles(); i < allele_count; ++i) {
            if (get_allele(i)->is_reference()) {
                return i;
            }
        }
        return -1;
    }

    /*!
     * @brief checks whether two allele lists are in fact the same
     * @param first
     * @param second
     * @return true if both list are equal
     */
    static bool equals(const InterfaceAlleleList<A>& first, const InterfaceAlleleList<A>& second)
    {
        if (first.number_of_alleles() != second.number_of_alleles()) {
            return false;
        }

        A a1, a2;
        for (size_t i = 0, len = first.number_of_alleles(); i < len; ++i) {
            a1 = first.get_allele(i);
            a2 = second.get_allele(i);
            if (!a1->equals(*a2)) {
                return false;
            }
        }

        return true;
    }

    InterfaceAlleleListPermutation<Ta>* permutation(const InterfaceAlleleList<Ta>* target, pMemoryPool pool)
    {
        if (equals(*this, *target)) {
            return new (pool->allocate(sizeof(NonPermutation))) NonPermutation{this};
        }
        else {
            return new (pool->allocate(sizeof(ActualPermutation))) ActualPermutation{this, target, pool};
        }
    }

    /* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
private:
    class NonPermutation : public InterfaceAlleleListPermutation<Ta>
    {
        /* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
    private:
        const InterfaceAlleleList<Ta>* _list;

        /* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
    public:
        explicit NonPermutation(const InterfaceAlleleList<Ta>* original)
            : _list(original)
        {}

        ~NonPermutation() override = default;
        size_t number_of_alleles() const override { return _list->number_of_alleles(); }
        Ta get_allele(size_t index) const override { return _list->get_allele(index); }
        const TaVector& get_alleles() const override { return _list->get_alleles(); }
        int32_t index_of_allele(Ta allele) const override { return _list->index_of_allele(allele); }
        bool is_partial() const override { return false; }
        bool is_non_permuted() const override { return true; }
        int32_t to_index(int32_t from_index) const override { return from_index; }
        int32_t from_index(int32_t to_index) const override { return to_index; }
        bool is_kept(int32_t from_index) const override
        {
            un_used(from_index);
            return true;
        }
        int32_t from_size() const override { return (int32_t)_list->number_of_alleles(); }
        int32_t to_size() const override { return (int32_t)_list->number_of_alleles(); }
        const TaVector& from_list() const override { return _list->get_alleles(); }
        const TaVector& to_list() const override { return _list->get_alleles(); }
    };

    class ActualPermutation : public InterfaceAlleleListPermutation<A>
    {
        /* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
    private:
        const InterfaceAlleleList<Ta>* _from;
        const InterfaceAlleleList<Ta>* _to;
        std::pmr::vector<int32_t> _from_index;
        std::pmr::vector<bool> _kept_from_indices;
        bool _non_permuted;
        bool _is_partial;

        /* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
    public:
        ActualPermutation(const InterfaceAlleleList<Ta>* original, const InterfaceAlleleList<Ta>* target, pMemoryPool pool)
            : _from(original)
            , _to(target)
            , _from_index(pool)
            , _kept_from_indices(pool)
            , _non_permuted(false)
            , _is_partial(false)
        {
            size_t to_size = target->number_of_alleles();
            size_t from_size = original->number_of_alleles();
            CHECK_CONDITION_EXIT(from_size < to_size, "target allele list is not a permutation of the original allele list");
            bool non_permuted = from_size == to_size;
            _is_partial = !non_permuted;

            _from_index.resize(to_size);
            std::fill(_from_index.begin(), _from_index.end(), 0);
            _kept_from_indices.resize(from_size);
            std::fill(_kept_from_indices.begin(), _kept_from_indices.end(), false);

            int32_t original_index;
            for (size_t i = 0; i < to_size; ++i) {
                original_index = original->index_of_allele(target->get_allele(i));
                CHECK_CONDITION_EXIT(original_index < 0, "target allele list is not a permutation of the original allele list");
                _kept_from_indices[original_index] = true;
                _from_index[i] = original_index;
                non_permuted &= ((size_t)original_index == i);
            }
            _non_permuted = non_permuted;
        }

        ~ActualPermutation() override = default;
        size_t number_of_alleles() const override { return _to->number_of_alleles(); }
        Ta get_allele(size_t index) const override { return _to->get_allele(index); }
        const TaVector& get_alleles() const override { return _to->get_alleles(); }
        int32_t index_of_allele(Ta allele) const override { return _to->index_of_allele(allele); }
        bool is_partial() const override { return _is_partial; }
        bool is_non_permuted() const override { return _non_permuted; }
        int32_t to_index(int32_t from_index) const override { return _to->index_of_allele(_from->get_allele(from_index)); }
        int32_t from_index(int32_t to_index) const override { return _from_index.at(to_index); }
        bool is_kept(int32_t from_index) const override { return _kept_from_indices.at(from_index); }
        int32_t from_size() const override { return _from->number_of_alleles(); }
        int32_t to_size() const override { return _to->number_of_alleles(); }
        const TaVector& from_list() const override { return _from->get_alleles(); }
        const TaVector& to_list() const override { return _to->get_alleles(); }
    };
};

template <typename A>
class InterfaceAlleleListPermutation : public InterfacePermutation<A>, public InterfaceAlleleList<A>
{
    /* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
public:
    typedef typename InterfaceAlleleList<A>::Ta Ta;
    typedef typename InterfaceAlleleList<A>::TaVector TaVector;
};

}  // namespace rovaca

#endif  // ROVACA_HC_INTERFACE_ALLELE_LIST_H_
