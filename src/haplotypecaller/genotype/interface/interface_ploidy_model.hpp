#ifndef ROVACA_HC_INTERFACE_PLOIDY_MODEL_H_
#define ROVACA_HC_INTERFACE_PLOIDY_MODEL_H_
#include "interface_sample_list.hpp"

namespace rovaca
{

class InterfacePloidyModel : public InterfaceSampleList
{
    /* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
public:
    ~InterfacePloidyModel() override = default;

    /*!
     * @brief return the assumed ploidy for a sample given its index
     * @param sample_index target sample index
     * @return 0 or greater
     */
    virtual int32_t sample_ploidy(size_t index) const = 0;

    /*!
     * @brief checks whether the ploidy is homogeneous across all samples
     * @return true if all samples has the same ploidy
     */
    virtual bool is_homogeneous() const = 0;

    /*!
     * @brief sum of all ploidy across all samples
     * @return 0 or greater
     */
    virtual int32_t total_ploidy() const = 0;
};

}  // namespace rovaca

#endif  // ROVACA_HC_INTERFACE_PLOIDY_MODEL_H_
