#ifndef ROVACA_HC_GENOTYPE_ALLELE_COUNTS_MANGER_H_
#define ROVACA_HC_GENOTYPE_ALLELE_COUNTS_MANGER_H_
#include <algorithm>
#include <memory>

#include "forward.h"
#include "genotype_allele_counts.h"
#include "genotype_macors.h"

namespace rovaca
{

static constexpr int32_t s_maximum_ploidy = 2;
static constexpr int32_t s_maximum_allele = 50;
static constexpr int32_t s_genotype_count_overflow = -1;
static constexpr int32_t s_maximum_strong_ref_genotype_per_ploidy = 1275;

class GenotypeAlleleCountsManger
{
    static constexpr size_t s_buffer_size = 1024 * 1024 * 100lu;

    /* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
private:
    uint8_t _data[s_buffer_size];
    std::pmr::monotonic_buffer_resource _pool;

    std::pmr::vector<std::pmr::vector<int32_t>> _allele_first_genotype_offset_by_ploidy;
    std::pmr::vector<std::pmr::vector<pGenotypeAlleleCounts>> _genotype_table_by_ploidy;

    /* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
public:
    DISALLOW_COPY_AND_ASSIGN(GenotypeAlleleCountsManger);

    static pGenotypeAlleleCountsManger get_instance()
    {
        static std::unique_ptr<GenotypeAlleleCountsManger> u(new GenotypeAlleleCountsManger{});
        return u.get();
    }

    const std::pmr::vector<int32_t>& allele_first_genotype_offset_by_ploidy(int32_t ploidy)
    {
        return _allele_first_genotype_offset_by_ploidy.at(ploidy);
    }

    const std::pmr::vector<pGenotypeAlleleCounts>& genotype_table_by_ploidy(int32_t ploidy) { return _genotype_table_by_ploidy.at(ploidy); }

    /* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
private:
    GenotypeAlleleCountsManger()
        : _data()
        , _pool(_data, s_buffer_size, std::pmr::null_memory_resource())
        , _allele_first_genotype_offset_by_ploidy(&_pool)
        , _genotype_table_by_ploidy(&_pool)
    {
        build_allele_first_genotype_offset_table();
        build_genotype_allele_counts_array();
    }

    /*!
     * @brief [ploidy][allele_count]: (C上p下a + C上1下a) -> 组合出来的基因型总数
     */
    void build_allele_first_genotype_offset_table()
    {
        int32_t row_count = s_maximum_ploidy + 1;
        int32_t col_count = s_maximum_allele + 1;
        _allele_first_genotype_offset_by_ploidy.resize(row_count);
        std::for_each(_allele_first_genotype_offset_by_ploidy.begin(), _allele_first_genotype_offset_by_ploidy.end(),
                      [&](std::pmr::vector<int32_t>& v) { v.resize(col_count); });
        std::fill(_allele_first_genotype_offset_by_ploidy.at(0).begin() + 1, _allele_first_genotype_offset_by_ploidy.at(0).end(), 1);
        auto& result = _allele_first_genotype_offset_by_ploidy;
        for (int32_t ploidy = 1; ploidy < row_count; ploidy++) {
            for (int32_t allele = 1; allele < col_count; allele++) {
                result[ploidy][allele] = result[ploidy][allele - 1] + result[ploidy - 1][allele];
                if (result[ploidy][allele] < result[ploidy][allele - 1]) {
                    result[ploidy][allele] = s_genotype_count_overflow;
                }
            }
        }
    }

    void build_genotype_allele_counts_array()
    {
        int32_t row_count = s_maximum_ploidy + 1;
        int32_t allele_count = s_maximum_allele;
        _genotype_table_by_ploidy.reserve(row_count);

        int32_t length, strong_ref_length;
        for (int32_t ploidy = 0; ploidy <= s_maximum_ploidy; ploidy++) {
            length = _allele_first_genotype_offset_by_ploidy.at(ploidy).at(allele_count);
            strong_ref_length = length == s_genotype_count_overflow ? s_maximum_strong_ref_genotype_per_ploidy
                                                                    : std::min(length, s_maximum_strong_ref_genotype_per_ploidy);
            std::pmr::vector<pGenotypeAlleleCounts> result_at_ploidy(&_pool);
            result_at_ploidy.reserve(strong_ref_length);
            result_at_ploidy.push_back(GenotypeAlleleCounts::first(ploidy, &_pool));
            for (int32_t genotype_index = 1; genotype_index < strong_ref_length; genotype_index++) {
                result_at_ploidy.push_back(result_at_ploidy.back()->next(&_pool));
            }
            _genotype_table_by_ploidy.push_back(std::move(result_at_ploidy));
        }
    }
};

}  // namespace rovaca

#endif  // ROVACA_HC_GENOTYPE_ALLELE_COUNTS_MANGER_H_
