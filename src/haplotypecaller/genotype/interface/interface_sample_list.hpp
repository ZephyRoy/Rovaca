#ifndef ROVACA_HC_INTERFACE_SAMPLE_LIST_H_
#define ROVACA_HC_INTERFACE_SAMPLE_LIST_H_
#include <string>

namespace rovaca
{

class InterfaceSampleList
{
    /* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
public:
    virtual ~InterfaceSampleList() = default;

    /*! @brief returns number of elements in the list */
    virtual size_t number_of_samples() const = 0;

    /*! @brief returns the element given its index within the set */
    virtual const std::string& get_sample(size_t index) const = 0;

    /*!
     * @brief returns the index of an object
     * @return -1 if such a sample is not an element of this set
     */
    virtual int32_t index_of_sample(const std::string& name) const = 0;
};

}  // namespace rovaca

#endif  // ROVACA_HC_INTERFACE_SAMPLE_LIST_H_
