#ifndef ROVACA_HC_HOMOGENEOUS_PLOIDY_MODEL_H_
#define ROVACA_HC_HOMOGENEOUS_PLOIDY_MODEL_H_
#include "forward.h"
#include "genotype_macors.h"
#include "interface/interface_ploidy_model.hpp"

namespace rovaca
{

class HomogeneousPloidyModel : public InterfacePloidyModel
{
    /* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
private:
    int32_t _ploidy;
    pInterfaceSampleList _sample_list;

    /* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
public:
    DISALLOW_COPY_AND_ASSIGN(HomogeneousPloidyModel);
    static pInterfacePloidyModel create(int32_t ploidy, pInterfaceSampleList sampleList)
    {
        return new HomogeneousPloidyModel(ploidy, sampleList);
    }

    ~HomogeneousPloidyModel() override = default;
    int32_t sample_ploidy([[maybe_unused]] size_t index) const override { return _ploidy; }
    bool is_homogeneous() const override { return true; }
    int32_t total_ploidy() const override { return _ploidy * (int32_t)_sample_list->number_of_samples(); }
    size_t number_of_samples() const override { return _sample_list->number_of_samples(); }
    const std::string& get_sample(size_t index) const override { return _sample_list->get_sample(index); }
    int32_t index_of_sample(const std::string& name) const override { return _sample_list->index_of_sample(name); }

    /* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
private:
    HomogeneousPloidyModel(int32_t ploidy, pInterfaceSampleList sampleList)
        : _ploidy(ploidy)
        , _sample_list(sampleList)
    {}
};

}  // namespace rovaca

#endif  // ROVACA_HC_HOMOGENEOUS_PLOIDY_MODEL_H_
