#ifndef ROVACA_HC_GENOTYPE_LIKELIHOOD_CALCULATOR_H_
#define ROVACA_HC_GENOTYPE_LIKELIHOOD_CALCULATOR_H_
#include <queue>

#include "forward.h"
#include "genotype_likelihoods.h"
#include "genotype_macors.h"
#include "interface/interface_likelihood_matrix.hpp"
#include "math_utils.h"

namespace rovaca
{

class GenotypeLikelihoodCalculator
{
    /* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
private:
    static pGenotypeAlleleCountsManger s_manger;

    int32_t _ploidy;
    int32_t _allele_count;
    int32_t _genotype_count;
    int32_t _maximum_distinct_alleles_in_genotype;
    pGenotypeAlleleCounts _last_overhead_counts{nullptr};
    std::priority_queue<int32_t, std::pmr::vector<int32_t>, std::less<>> _allele_heap;

    pMemoryPool _pool{nullptr};

    /* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
public:
    static pGenotypeLikelihoodCalculator create(int32_t ploidy, int32_t allele_count, pMemoryPool pool);

    static int32_t genotype_count(int32_t ploidy, int32_t allele_count);

    /*! @brief 返回cache中对应index的pGenotypeAlleleCounts，如果index大于cache，则返回一个副本 */
    pGenotypeAlleleCounts genotype_allele_counts_at(int32_t index);

    int32_t ploidy() const { return _ploidy; };
    int32_t allele_count() const { return _allele_count; };
    int32_t genotype_count() const { return _genotype_count; };

    /*! @brief calculate the likelihoods given the list of alleles and the likelihood map */
    template <typename E, typename A>
    pGenotypeLikelihoods genotype_likelihoods(interfaceLikelihoodMatrix<E, A>* likelihoods)
    {
        DoubleVector read_likelihoods = get_read_raw_read_likelihoods_by_genotype_index(likelihoods);
        return GenotypeLikelihoods::create(std::move(read_likelihoods), _pool);
    }

    /**
     * Composes a genotype index map given a allele index recoding.
     * @param old_to_new_allele_index_map allele recoding. The ith entry indicates the index of the allele in original encoding that
     *                                    corresponds to the ith allele index in the final encoding.
     * @return never nullptr
     */
    Int32Vector genotype_index_map(const Int32Vector& old_to_new_allele_index_map);

    /**
     * Returns the likelihood index given the allele counts.
     * @param alleleCountArray the query allele counts. This must follow the format returned by
     *  {@link GenotypeAlleleCounts#copyAlleleCounts} with 0 offset.
     *
     *  <ul>
     *      <li>is {@code null},</li>
     *      <li>or its length is not even,</li>
     *      <li>or it contains any negatives,
     *      <li>or the count sum does not match the calculator ploidy,</li>
     *      <li>or any of the alleles therein is negative or greater than the maximum allele index.</li>
     *  </ul>
     *
     * @return 0 or greater but less than {@link #genotypeCount}.
     */
    int32_t allele_counts_to_index(const Int32Vector& allele_count_array);

    /*!
     * @brief give a list of alleles, returns the likelihood array index.
     * @param allele_count_array
     * @return
     */
    int32_t alleles_to_index(const Int32Vector& allele_count_array);

    /* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
private:
    GenotypeLikelihoodCalculator(int32_t ploidy, int32_t allele_count, int32_t genotype_count, pMemoryPool pool);

    /*!
     * @brief a helper method that actually does the matrix operations but returns the raw values.
     * @return the raw array (in log10 likelihoods space) of the gl for each genotype
     */
    template <typename E, typename A>
    DoubleVector get_read_raw_read_likelihoods_by_genotype_index(interfaceLikelihoodMatrix<E, A>* likelihoods)
    {
        CHECK_CONDITION_EXIT((size_t)_allele_count != likelihoods->number_of_alleles(), "mismatch between allele list and allele_count");
        size_t evidence_count = likelihoods->evidence_count();

        /*!
         * @brief 此处需要用pool开辟三个buffer
         * read_allele_likelihood_by_allele_count[evidence_count * _allele_count * (_ploidy + 1)]
         * read_likelihoods_by_genotype_index[_genotype_count][evidence_count]
         * read_genotype_likelihood_components[_ploidy * evidence_count]
         * genotype_likelihoods_result[_genotype_count]
         */
        DoubleVector read_allele_likelihood_by_allele_count(evidence_count * _allele_count * (_ploidy + 1), _pool);
        DoubleVector2D read_likelihoods_by_genotype_index(_genotype_count, DoubleVector(evidence_count), _pool);
        DoubleVector read_genotype_likelihood_components(_ploidy * evidence_count, _pool);

        DoubleVector genotype_likelihoods_result(_genotype_count, _pool);

        /// [x][y][z] = z * LnLk(Read_x | Allele_y)
        read_likelihood_components_by_allele_count(likelihoods, read_allele_likelihood_by_allele_count);

        genotype_likelihood_by_read(read_allele_likelihood_by_allele_count, evidence_count, read_likelihoods_by_genotype_index,
                                    read_genotype_likelihood_components);
        genotype_likelihoods(read_likelihoods_by_genotype_index, evidence_count, genotype_likelihoods_result);

        return genotype_likelihoods_result;
    }

    /**
     * @brief 将read_likelihood_components_by_allele_count分为 _ploidy + 1 大段，每大段又分为 _ploidy + 1 小段
     * @example 当_ploidy==2时，read_likelihood_components_by_allele_count数组逻辑上被分为3大段，每大段被分为3小段
     * 每个大段中有三个小段，循环至数组结尾
     *      第一小段为默认值
     *      第二段在likelihoods中copy的对应allele_index的likelihoods
     *      第三段为第二段基础上每个值加MathUtils::log10(frequency)
     */
    template <typename E, typename A>
    void read_likelihood_components_by_allele_count(interfaceLikelihoodMatrix<E, A>* likelihoods,
                                                    DoubleVector& read_likelihood_components_by_allele_count)
    {
        int32_t read_count = (int32_t)likelihoods->evidence_count();
        int32_t allele_data_size = read_count * (_ploidy + 1);

        // frequency1offset = read_count to skip the useless frequency == 0. so now we are at the start frequency == 1
        // frequency1offset += allele_data_size to skip to the next allele index data location (+ read_count) at each iteration.
        for (int32_t a = 0, frequency1offset = read_count; a < _allele_count; a++, frequency1offset += allele_data_size) {
            likelihoods->copy_allele_likelihoods(a, frequency1offset, read_likelihood_components_by_allele_count);

            // p = 2 because the frequency == 1 we already have it.
            for (int32_t frequency = 2, destination_offset = frequency1offset + read_count; frequency <= _ploidy; frequency++) {
                double log10frequency = MathUtils::log10(frequency);
                for (int32_t r = 0, source_offset = frequency1offset; r < read_count; r++) {
                    read_likelihood_components_by_allele_count[destination_offset++] =
                        read_likelihood_components_by_allele_count[source_offset++] + log10frequency;
                }
            }
        }
    }

    void genotype_likelihood_by_read(const DoubleVector& pre_result, size_t read_count, DoubleVector2D& result, DoubleVector& buffer);

    void genotype_likelihoods(const DoubleVector2D& pre_result, size_t read_count, DoubleVector& genotype_likelihoods_result) const;

    /*!
     * @brief calculates the likelihood component by read for a given genotype allele count assuming that there are exactly one allele
     * present in the genotype.
     */
    void single_component_genotype_likelihood_by_read(pGenotypeAlleleCounts gac, const DoubleVector& pre_result, int32_t read_count,
                                                      DoubleVector& result) const;
    /*!
     * @brief calculates the likelihood component by read for a given genotype allele count assuming that there are exactly two alleles
     * present in the genotype (with arbitrary non-zero counts each).
     */
    void two_component_genotype_likelihood_by_read(pGenotypeAlleleCounts gac, const DoubleVector& pre_result, int32_t read_count,
                                                   DoubleVector& result) const;
    /*!
     * @brief general genotype likelihood component by read calculator. it does not make any assumption in the exact number of alleles
     * present in the genotype.
     */
    void many_component_genotype_likelihood_by_read(pGenotypeAlleleCounts gac, const DoubleVector& pre_result, int32_t read_count,
                                                    DoubleVector& result, DoubleVector& buffer) const;

    pGenotypeAlleleCounts next_genotype_allele_counts(pGenotypeAlleleCounts gac);

    /*!
     * @brief Performs the genotype mapping per new genotype index.
     * @param idx
     * @param gac
     * @param old2new_map
     * @param dest
     * @param buffer
     */
    void genotype_index_map_per_genotype_index(int32_t idx, pGenotypeAlleleCounts gac, const Int32Vector& old2new_map, Int32Vector& dest,
                                               Int32Vector& buffer);

    /**
     * Transforms the content of the heap into an index.
     * <p>
     *     The heap contents are flushed as a result, so is left ready for another use.
     * </p>
     * @return a valid likelihood index.
     */
    int32_t allele_heap_to_index();
};

}  // namespace rovaca

#endif  // ROVACA_HC_GENOTYPE_LIKELIHOOD_CALCULATOR_H_
