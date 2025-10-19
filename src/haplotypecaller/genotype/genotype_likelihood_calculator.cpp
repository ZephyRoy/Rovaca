#include "genotype_likelihood_calculator.h"

#include <memory>
#include <memory_resource>

#include "genotype_allele_counts.h"
#include "genotype_allele_counts_manger.hpp"
#include "rovaca_logger.h"

namespace rovaca
{

pGenotypeAlleleCountsManger GenotypeLikelihoodCalculator::s_manger = GenotypeAlleleCountsManger::get_instance();

pGenotypeLikelihoodCalculator GenotypeLikelihoodCalculator::create(int32_t ploidy, int32_t allele_count, pMemoryPool pool)
{
    CHECK_CONDITION_EXIT(ploidy != 2 || allele_count >= s_maximum_allele,
                         "the number of genotypes is too large for ploidy {} and allele {}", ploidy, allele_count);

    int32_t genotype_count = s_manger->allele_first_genotype_offset_by_ploidy(ploidy).at(allele_count);
    return new ALLOC_TYPE_IN_POOL(pool, GenotypeLikelihoodCalculator)
        GenotypeLikelihoodCalculator{ploidy, allele_count, genotype_count, pool};
}

int32_t GenotypeLikelihoodCalculator::genotype_count(int32_t ploidy, int32_t allele_count)
{
    CHECK_CONDITION_EXIT(allele_count >= (s_maximum_ploidy + 1), "the number of genotypes is too large");
    return s_manger->allele_first_genotype_offset_by_ploidy(ploidy).at(allele_count);
}

pGenotypeAlleleCounts GenotypeLikelihoodCalculator::genotype_allele_counts_at(int32_t index)
{
    CHECK_CONDITION_EXIT(index < 0 || index >= _genotype_count,
                         "invalid index. index = {}, genotype_count = {}, ploidy = {}, allele_count = {}", index, _genotype_count, _ploidy,
                         _allele_count);
    if (index < s_maximum_strong_ref_genotype_per_ploidy) {
        return s_manger->genotype_table_by_ploidy(_ploidy).at(index);
    }
    else if (nullptr == _last_overhead_counts || _last_overhead_counts->index() > index) {
        _last_overhead_counts = s_manger->genotype_table_by_ploidy(_ploidy).at(s_maximum_strong_ref_genotype_per_ploidy - 1)->copy(_pool);
        _last_overhead_counts->increase(index - s_maximum_strong_ref_genotype_per_ploidy + 1);
    }
    else {
        _last_overhead_counts->increase(index - _last_overhead_counts->index());
    }
    return _last_overhead_counts->copy(_pool);
}

GenotypeLikelihoodCalculator::GenotypeLikelihoodCalculator(int32_t ploidy, int32_t allele_count, int32_t genotype_count, pMemoryPool pool)
    : _ploidy(ploidy)
    , _allele_count(allele_count)
    , _genotype_count(genotype_count)
    , _maximum_distinct_alleles_in_genotype(std::min(ploidy, allele_count))
    , _allele_heap(std::less<>(), Int32Vector(pool))
    , _pool(pool)
{}

void GenotypeLikelihoodCalculator::genotype_likelihood_by_read(const DoubleVector& pre_result, size_t read_count, DoubleVector2D& result,
                                                               DoubleVector& buffer)
{
    // here we don't use the convenience of {@link #genotype_allele_counts_at(int)} within the loop to spare instantiations of
    // genotype_allele_counts class when we are dealing with many genotypes
    pGenotypeAlleleCounts gac = s_manger->genotype_table_by_ploidy(_ploidy).at(0);
    int32_t component_count;
    for (int32_t genotype_index = 0; genotype_index < _genotype_count; genotype_index++) {
        DoubleVector& read_likelihoods = result.at(genotype_index);
        component_count = gac->distinct_allele_count();
        switch (component_count) {
            case 1: single_component_genotype_likelihood_by_read(gac, pre_result, (int32_t)read_count, read_likelihoods); break;
            case 2: two_component_genotype_likelihood_by_read(gac, pre_result, (int32_t)read_count, read_likelihoods); break;
            default: many_component_genotype_likelihood_by_read(gac, pre_result, (int32_t)read_count, read_likelihoods, buffer); break;
        }
        if (ROVACA_LIKELY(genotype_index < _genotype_count - 1)) {
            gac = next_genotype_allele_counts(gac);
        }
    }
}

void GenotypeLikelihoodCalculator::genotype_likelihoods(const DoubleVector2D& pre_result, size_t read_count,
                                                        DoubleVector& genotype_likelihoods_result) const
{
    double denominator = (double)read_count * MathUtils::log10(_ploidy);
    // instead of dividing each read likelihood by ploidy ( so subtract log10(ploidy) )
    // we multiply them all and the divide by ploidy^read_count (so substract read_count * log10(ploidy))
    for (int32_t g = 0; g < _genotype_count; g++) {
        genotype_likelihoods_result[g] = std::accumulate(pre_result.at(g).begin(), pre_result.at(g).end(), 0.0) - denominator;
    }
}

void GenotypeLikelihoodCalculator::single_component_genotype_likelihood_by_read(pGenotypeAlleleCounts gac, const DoubleVector& pre_result,
                                                                                int32_t read_count, DoubleVector& result) const
{
    int32_t allele = gac->allele_index_at(0);
    int32_t offset = (allele * (_ploidy + 1) + _ploidy) * read_count;
    for (int32_t r = 0; r < read_count; r++) {
        result.at(r) = pre_result.at(offset++);
    }
}

void GenotypeLikelihoodCalculator::two_component_genotype_likelihood_by_read(pGenotypeAlleleCounts gac, const DoubleVector& pre_result,
                                                                             int32_t read_count, DoubleVector& result) const
{
    int32_t allele0 = gac->allele_index_at(0);
    int32_t freq0 = gac->allele_count_at(0);
    int32_t allele1 = gac->allele_index_at(1);
    int32_t freq1 = _ploidy - freq0;
    int32_t allele0ln_lk_offset = read_count * ((_ploidy + 1) * allele0 + freq0);
    int32_t allele1ln_lk_offset = read_count * ((_ploidy + 1) * allele1 + freq1);
    double ln_lk0, ln_lk1;
    for (int32_t r = 0; r < read_count; r++) {
        ln_lk0 = pre_result.at(allele0ln_lk_offset++);
        ln_lk1 = pre_result.at(allele1ln_lk_offset++);
        result.at(r) = MathUtils::approximate_log10sum_log10(ln_lk0, ln_lk1);
    }
}

void GenotypeLikelihoodCalculator::many_component_genotype_likelihood_by_read(pGenotypeAlleleCounts gac, const DoubleVector& pre_result,
                                                                              int32_t read_count, DoubleVector& result,
                                                                              DoubleVector& buffer) const
{
    Int32Vector genotype_alleles_and_counts(_maximum_distinct_alleles_in_genotype, buffer.get_allocator());
    int32_t component_count = gac->distinct_allele_count();
    int32_t allele_data_size = (_ploidy + 1) * read_count;

    int32_t allele_index, allele_count, allele_data_offset;
    for (int32_t c = 0, cc = 0; c < component_count; c++) {
        allele_index = genotype_alleles_and_counts.at(cc++);
        allele_count = genotype_alleles_and_counts.at(cc++);
        // allele_data_offset will point to the index of the first read likelihood for that allele and allele count.
        allele_data_offset = allele_data_size * allele_index + allele_count * read_count;
        for (int32_t r = 0, read_data_offset = c; r < read_count; r++, read_data_offset += _maximum_distinct_alleles_in_genotype) {
            buffer.at(read_data_offset) = pre_result.at(allele_data_offset++);
        }
    }

    // calculate the likelihood per read.
    for (int32_t r = 0, read_data_offset = 0; r < read_count; r++, read_data_offset += _maximum_distinct_alleles_in_genotype) {
        result.at(r) = MathUtils::approximate_log10sum_log10(buffer, read_data_offset, read_data_offset + component_count);
    }
}

pGenotypeAlleleCounts GenotypeLikelihoodCalculator::next_genotype_allele_counts(pGenotypeAlleleCounts gac)
{
    pGenotypeAlleleCounts result;
    int32_t index = gac->index();
    int32_t cmp = index - s_maximum_strong_ref_genotype_per_ploidy + 1;
    if (cmp < 0) {
        result = s_manger->genotype_table_by_ploidy(_ploidy).at(index + 1);
    }
    else if (cmp == 0) {
        result = s_manger->genotype_table_by_ploidy(_ploidy).at(index)->copy(_pool);
        result->increase();
    }
    else {
        gac->increase();
        result = gac;
    }
    return result;
}

Int32Vector GenotypeLikelihoodCalculator::genotype_index_map(const Int32Vector& old_to_new_allele_index_map)
{
    int32_t result_allele_count = (int32_t)old_to_new_allele_index_map.size();
    CHECK_CONDITION_EXIT(result_allele_count > _allele_count, "this calculator does not have enough capacity");
    int32_t result_length = result_allele_count == _allele_count
                                ? _genotype_count
                                : s_manger->allele_first_genotype_offset_by_ploidy(_ploidy).at(result_allele_count);

    Int32Vector result(result_length, old_to_new_allele_index_map.get_allocator());
    Int32Vector sorted_allele_counts(std::max(_ploidy, result_allele_count) << 1, old_to_new_allele_index_map.get_allocator());
    while (!_allele_heap.empty()) {  // 居然没有 clear ?
        _allele_heap.pop();
    }
    pGenotypeAlleleCounts gac = s_manger->genotype_table_by_ploidy(_ploidy).at(0);
    for (int32_t i = 0; i < result_length; ++i) {
        genotype_index_map_per_genotype_index(i, gac, old_to_new_allele_index_map, result, sorted_allele_counts);
        if (i < result_length - 1) {
            gac = next_genotype_allele_counts(gac);
        }
    }
    return result;
}

int32_t GenotypeLikelihoodCalculator::allele_counts_to_index(const Int32Vector& allele_count_array)
{
    CHECK_CONDITION_EXIT((allele_count_array.size() & 1) != 0, "the allele counts array cannot have odd length");
    while (!_allele_heap.empty()) {  // 居然没有 clear ?
        _allele_heap.pop();
    }
    int index, count;
    for (int i = 0, len = (int32_t)allele_count_array.size(); i < len; i += 2) {
        index = allele_count_array[i];
        count = allele_count_array[i + 1];
        CHECK_CONDITION_EXIT(count < 0, "no allele count can be less than 0");
        for (int32_t j = 0; j < count; j++) {
            _allele_heap.push(index);
        }
    }
    return allele_heap_to_index();
}

int32_t GenotypeLikelihoodCalculator::alleles_to_index(const Int32Vector& allele_count_array)
{
    if (0 == _ploidy) {
        return 0;
    }

    while (!_allele_heap.empty()) {  // 居然没有 clear ?
        _allele_heap.pop();
    }
    std::for_each(allele_count_array.begin(), allele_count_array.end(), [&](int32_t i) { _allele_heap.push(i); });
    return allele_heap_to_index();
}

void GenotypeLikelihoodCalculator::genotype_index_map_per_genotype_index(int32_t idx, pGenotypeAlleleCounts gac,
                                                                         const Int32Vector& old2new_map, Int32Vector& dest,
                                                                         Int32Vector& buffer)
{
    int32_t distinct_allele_count = gac->distinct_allele_count();
    gac->copy_allele_counts(0, buffer);
    int32_t old_index, repeats, new_index;
    for (int32_t j = 0, jj = 0; j < distinct_allele_count; ++j) {
        old_index = buffer.at(jj++);
        repeats = buffer.at(jj++);
        new_index = old2new_map.at(old_index);
        CHECK_CONDITION_EXIT(new_index < 0 || new_index >= _allele_count, "found invalid new allele index");
        for (int32_t k = 0; k < repeats; k++) {
            _allele_heap.push(new_index);
        }
    }
    int32_t genotype_index = allele_heap_to_index();
    dest.at(idx) = genotype_index;
}

int32_t GenotypeLikelihoodCalculator::allele_heap_to_index()
{
    CHECK_CONDITION_EXIT(_allele_heap.size() != (size_t)_ploidy, "the sum of allele counts must be equal to the ploidy of the calculator");
    CHECK_CONDITION_EXIT(_allele_heap.top() >= _allele_count, "invalid allele {} more than the maximum {}.", _allele_heap.top(),
                         _allele_count);
    int32_t top, result = 0;
    for (int32_t p = _ploidy; p > 0; --p) {
        top = _allele_heap.top();
        _allele_heap.pop();
        CHECK_CONDITION_EXIT(top < 0, "invalid allele index: {} must be equal or greater than 0", top);
        result += s_manger->allele_first_genotype_offset_by_ploidy(p).at(top);
    }
    return result;
}

}  // namespace rovaca