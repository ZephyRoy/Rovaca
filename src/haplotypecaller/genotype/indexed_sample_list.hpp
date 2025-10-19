#ifndef ROVACA_HC_INDEXED_SAMPLE_LIST_H_
#define ROVACA_HC_INDEXED_SAMPLE_LIST_H_
#include <map>
#include <string>
#include <vector>

#include "rovaca_logger.h"
#include "genotype_macors.h"
#include "interface/interface_sample_list.hpp"

namespace rovaca
{

/*!
 * @brief 此类在程序运行期间仅需生成一次
 */
class IndexedSampleList : public InterfaceSampleList
{
    /* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
private:
    std::vector<std::string> _samples;
    std::map<std::string, size_t> _name2idx_map;

    /* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
public:
    static InterfaceSampleList* create(const std::vector<std::string>& samples) { return new IndexedSampleList{samples}; }

    ~IndexedSampleList() override = default;
    DISALLOW_COPY_AND_ASSIGN(IndexedSampleList);

    size_t number_of_samples() const override { return _samples.size(); }

    const std::string& get_sample(size_t index) const override
    {
        CHECK_CONDITION_EXIT(index >= _samples.size(), "out of range");
        return _samples.at(index);
    }

    int32_t index_of_sample(const std::string& name) const override
    {
        if (ROVACA_UNLIKELY(!_name2idx_map.count(name))) {
            return -1;
        }
        return (int32_t)_name2idx_map.at(name);
    }

    /* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
private:
    explicit IndexedSampleList(const std::vector<std::string>& samples)
        : _samples()
        , _name2idx_map()
    {
        for (size_t i = 0, index = 0, sn = samples.size(); i < sn; ++i) {
            _samples.reserve(samples.size());
            if (!_name2idx_map.count(samples.at(i))) {
                _samples.push_back(samples.at(i));
                _name2idx_map.insert({samples.at(i), index});
                index++;
            }
        }
    }
};

}  // namespace rovaca

#endif  // ROVACA_HC_INDEXED_SAMPLE_LIST_H_
