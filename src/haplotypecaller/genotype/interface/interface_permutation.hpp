#ifndef ROVACA_HC_INTERFACE_PERMUTATION_H_
#define ROVACA_HC_INTERFACE_PERMUTATION_H_
#include <cstdint>
#include <vector>

namespace rovaca
{

template <typename E>
class InterfacePermutation
{
    /* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
public:
    virtual ~InterfacePermutation() = default;

    /*!
     * @brief checks whether this permutation is a partial one of the original list
     * @note a partial permutation is one in that not all original elements take part in
     * @return true if this is a partial permutation
     */
    virtual bool is_partial() const = 0;

    /*!
     * @brief checks whether this is a trivial permutation where the resulting element list is the same as original
     * @return true if the resulting element list is the same as the original
     */
    virtual bool is_non_permuted() const = 0;

    /*!
     * @brief given an index on the original list, returns the position of tha element in the resulting list
     * @param from_index the query original element index
     * @return -1 if that element is not part of the result (partial) permutation, otherwise some number between 0 and [to_size()-1]
     */
    virtual int32_t to_index(int32_t from_index) const = 0;

    /*!
     * @brief given an index on the resulting list, it gives you the index of that element on the original list
     * @param to_index the query resulting list index
     * @return a value between 0 and [from_size()-1]
     */
    virtual int32_t from_index(int32_t to_index) const = 0;

    /*!
     * @brief given an index of the original list, return whether this index is found at any position of the permuted list
     * @param from_index this is trivial if the permutation is not partial
     * @return
     */
    virtual bool is_kept(int32_t from_index) const = 0;

    /**
     * @brief length of the original element list.
     * @return 0 or greater.
     */
    virtual int32_t from_size() const = 0;

    /**
     * @brief length of the resulting element list.
     * @return 0 or greater.
     */
    virtual int32_t to_size() const = 0;

    /*!
     * @brief returns an unmodifiable view to the original element list
     * @return
     */
    virtual const std::pmr::vector<E>& from_list() const = 0;

    /*!
     * @brief returns an unmodifiable view to the original element list
     * @return
     */
    virtual const std::pmr::vector<E>& to_list() const = 0;
};

}  // namespace rovaca

#endif  // ROVACA_HC_INTERFACE_PERMUTATION_H_
