#ifndef ROVACA_HC_INDEXED_ALLELE_LIST_H_
#define ROVACA_HC_INDEXED_ALLELE_LIST_H_
#include <vector>

#include "rovaca_logger.h"
#include "forward.h"
#include "genotype_macors.h"
#include "interface/interface_allele_list.hpp"

namespace rovaca
{

template <typename A>
class IndexedAlleleList : public InterfaceAlleleList<A>
{
    /* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
public:
    typedef typename InterfaceAlleleList<A>::Ta Ta;
    typedef typename InterfaceAlleleList<A>::TaVector TaVector;

    /* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
private:
    std::pmr::vector<Ta> _alleles;
    std::pmr::map<Ta, size_t> _type2idx_map;

    /* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
public:
    /*! @brief 对象的唯一生成接口，确保所有对象均源于pMemoryPool */
    static InterfaceAlleleList<Ta>* create(const std::pmr::vector<Ta>& alleles, pMemoryPool pool)
    {
        return new ALLOC_TYPE_IN_POOL(pool, IndexedAlleleList<Ta>) IndexedAlleleList<Ta>{alleles, pool};
    }
    static InterfaceAlleleList<Ta>* create(pMemoryPool pool)
    {
        return new ALLOC_TYPE_IN_POOL(pool, IndexedAlleleList<Ta>) IndexedAlleleList<Ta>{pool};
    }

    ~IndexedAlleleList() override = default;

    size_t number_of_alleles() const override { return _alleles.size(); }

    Ta get_allele(size_t index) const override
    {
        CHECK_CONDITION_EXIT(index >= _alleles.size(), "out of range");
        return _alleles.at(index);
    }

    const TaVector& get_alleles() const override { return _alleles; }

    int32_t index_of_allele(Ta allele) const override
    {
        if (ROVACA_UNLIKELY(!_type2idx_map.count(allele))) {
            return -1;
        }
        return (int32_t)_type2idx_map.at(allele);
    }

    /* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
private:
    explicit IndexedAlleleList(pMemoryPool pool)
        : _alleles(pool)
        , _type2idx_map(pool)
    {}
    IndexedAlleleList(const std::pmr::vector<Ta>& alleles, pMemoryPool pool)
        : _alleles(pool)
        , _type2idx_map(pool)
    {
        _alleles.reserve(alleles.size());
        for (size_t i = 0, index = 0, an = alleles.size(); i < an; ++i) {
            if (!_type2idx_map.count(alleles.at(i))) {
                _alleles.emplace_back(alleles.at(i));
                _type2idx_map.insert({alleles.at(i), index});
                ++index;
            }
        }
    }
};

}  // namespace rovaca

#endif  // ROVACA_HC_INDEXED_ALLELE_LIST_H_
