#ifndef ROVACA_HC_GENOTYPE_ALLELE_COUNTS_H_
#define ROVACA_HC_GENOTYPE_ALLELE_COUNTS_H_
#include <cstdint>
#include <functional>
#include <vector>

#include "forward.h"
#include "genotype_macors.h"

namespace rovaca
{

/**
 * Collection of allele counts for a genotype. It encompasses what alleles are present in the genotype and in what number.</p>
 *
 * <p>Alleles are represented herein by their indices running from <b>0</b> to <b>N-1</b> where <i>N</i> is the number of alleles.</p>
 *
 * <p>Genotypes are represented as a single array of alternating alleles and counts, where only alleles with non-zero counts are included:
 * [allele 1, count1, allele 2, count2. . .]</p>
 *
 * <p>Each allele present in a genotype (count != 0) has a <i>rank</i>, that is the 0-based ordinal of
 * that allele amongst the ones present in the genotype as sorted by their index.</p>
 *
 * <p>For example:</p>
 *
 * <p><b>[0,1,2,1]</b> has two alleles with indices <b>0</b> and <b>2</b>, both with count 1, corresponding to diploid genotype 0/2.
 * The rank of <b>0</b> is <i>0</i> whereas the rank of <b>2</b> is <i>1</i>.</p>
 *
 * <p><b>[2,1,4,2,7,1]</b> has three alleles with indices <b>2</b>, <b>4</b> and <b>7</b>. <b>2</b> and <b>7</b> have count 1 whereas
 * <b>4</b> has count 2. It corresponds to tetraploid genotype 2/4/4/7 The rank of <b>2</b> is <i>0</i>, the rank of <b>4</b> is <i>1</i>.
 * and the rank of <b>7</b> is <i>2</i>.</p>
 *
 * <p>In contrast, in both examples above both <b>3</b> and <b>10</b> (and many others) are absent thus they have no rank (represented by
 * <i>-1</i> whenever applies).</p>
 *
 * <p><b>[0,0,1,2]</b> is not valid because allele 0 has a count of 0 and should be absent from the array.</p>
 *
 * <p><b>[1,1,0,1]</b> is not valid because allele 1 comes before allele 0.</p>
 *
 * <p>{@link GenotypeAlleleCounts} instances have themselves their own index (returned by {@link #index() index()}, that indicate their
 * 0-based ordinal within the possible genotype combinations with the same ploidy.</p>
 *
 * <p>For example, for ploidy 3:</p>
 *
 * <table>
 *     <th>Index</th><th>Genotype</th>
 *     <tr><td>0</td><td><b>0/0/0</b></td></tr>
 *     <tr><td>1</td><td><b>0/0/1</b></td></tr>
 *     <tr><td>2</td><td><b>0/1/1</b></td></tr>
 *     <tr><td>3</td><td><b>1/1/1</b></td></tr>
 *     <tr><td>4</td><td><b>0/0/2</b></td></tr>
 *     <tr><td>6</td><td><b>0/1/2</b></td></tr>
 *     <tr><td>7</td><td><b>1/1/2</b></td></tr>
 *     <tr><td>8</td><td><b>0/2/2</b></td></tr>
 *     <tr><td>9</td><td><b>1/2/2</b></td></tr>
 *     <tr><td>10</td><td><b>2/2/2</b></td></tr>
 *     <tr><td>11</td><td><b>0/0/3</b></td></tr>
 *     <tr><td>12</td><td><b>0/1/3</b></td></tr>
 *     <tr><td>13</td><td><b>1/1/3</b></td></tr>
 *     <tr><td>14</td><td><b>0/2/3</b></td></tr>
 *     <tr><td>15</td><td><b>1/2/3</b></td></tr>
 *     <tr><td>16</td><td><b>2/2/3</b></td></tr>
 *     <tr><td>17</td><td><b>0/3/3</b></td></tr>
 *     <tr><td>...</td><td>...</td></tr>
 * </table>
 *
 * The total number of possible genotypes is only bounded by the maximum allele index.
 */
class GenotypeAlleleCounts
{
    /* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
private:
    int32_t _index;                  // index of this genotype within genotypes of the same ploidy and number of alleles
    int32_t _ploidy;
    int32_t _distinct_allele_count;  // number of different alleles in the genotype
    double _log10combination_count;
    Int32Vector _sorted_allele_counts;

    /* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
public:
    DISALLOW_COPY_AND_ASSIGN(GenotypeAlleleCounts);

    /*! @brief GenotypeAlleleCounts 的原始构建方法，有了 first，其它的通过调用 next 构造 */
    static pGenotypeAlleleCounts first(int32_t ploidy, pMemoryPool pool);
    /*! @brief calculates the next genotype in likelihood indexing order */
    pGenotypeAlleleCounts next(pMemoryPool pool);

    int32_t ploidy() const { return _ploidy; }
    int32_t index() const { return _index; }
    double log10combination_count();
    int32_t distinct_allele_count() const { return _distinct_allele_count; }
    int32_t allele_index_at(int32_t rank) const { return _sorted_allele_counts.at(rank << 1); }
    int32_t allele_count_at(int32_t rank) const { return _sorted_allele_counts.at((rank << 1) + 1); }
    int32_t allele_rank_for(int32_t index) const { return allele_index_to_rank(index, 0, _distinct_allele_count); }
    int32_t allele_count_for(int32_t index) const;
    bool contains_allele(int32_t index) const { return allele_rank_for(index) >= 0; }
    void copy_allele_counts(int32_t offset, Int32Vector& dest) const;
    Int32Vector allele_counts_by_index(int32_t maximum_allele_index, pMemoryPool pool);

    /*!
     * @brief updates the genotype counts to match the next genotype according to the canonical ordering of pls
     * @note 调用此方法的 GenotypeAlleleCounts 对象必须是 copy 出来的，不能使用 cache 中的 GenotypeAlleleCounts 调用
     */
    void increase();

    /*!
     * @brief increases the allele counts a number of times
     * @note 调用此方法的 GenotypeAlleleCounts 对象必须是 copy 出来的，不能使用 cache 中的 GenotypeAlleleCounts 调用
     */
    void increase(int32_t times);

    /**
     * returns the largest allele index present in the genotype.
     * @return -1 if there is no alleles (ploidy == 0), 0 or greater otherwise.
     */
    int32_t maximum_allele_index() const
    {
        return _distinct_allele_count == 0 ? -1 : _sorted_allele_counts.at((_distinct_allele_count - 1) << 1);
    }

    /**
     * returns the smallest allele index present in the genotype.
     * @return -1 if there is no allele (ploidy == 0), 0 or greater otherwise.
     */
    int32_t minimum_allele_index() const { return _distinct_allele_count == 0 ? -1 : _sorted_allele_counts.at(0); }

    /*! @brief creates an independent copy of this genotype_allele_counts */
    pGenotypeAlleleCounts copy(pMemoryPool pool) const;

    double sum_over_allele_indices_and_counts(const std::function<double(int32_t, int32_t)>& func);

    void for_each_allele_index_and_count(const std::function<void(int32_t, int32_t)>& func);

    /**
     * perform an action for every allele index not represented in this genotype.  for example if the total allele count is 4 and {@code
     * sorted_allele_counts} is [0,1,2,1] then alleles 0 and 2 are present, each with a count of 1, while alleles 1 and 3 are absent, so we
     * perform {@code action} on 1 and 3.
     */
    void for_each_absent_allele_index(const std::function<void(int32_t)>& func, int32_t allele_count);

    /**
     * composes a list with the alleles, possibly containing repeats i.e. if internally this stores
     * allele 0 count = 1, allele 2 count = 2, the output is [allele0, allele2, allele2]
     * @param alleles_to_use alleles to use.
     * @return never null, but it might be restricted (unmodifiable or non-expandable).
     */
    AlleleVector as_allele_list(const AlleleVector& alleles_to_use);

    /* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
private:
    GenotypeAlleleCounts(int32_t index, int32_t ploidy, const Int32Vector& sorted_allele_counts, pMemoryPool pool)
        : GenotypeAlleleCounts(index, ploidy, sorted_allele_counts, (int32_t)sorted_allele_counts.size() / 2, pool)
    {}

    GenotypeAlleleCounts(int32_t index, int32_t ploidy, Int32Vector&& sorted_allele_counts)
        : GenotypeAlleleCounts(index, ploidy, (int32_t)sorted_allele_counts.size() / 2, std::forward<Int32Vector>(sorted_allele_counts))
    {}

    GenotypeAlleleCounts(int32_t index, int32_t ploidy, int32_t distinct_allele_count, Int32Vector&& sorted_allele_counts)
        : _index(index)
        , _ploidy(ploidy)
        , _distinct_allele_count(distinct_allele_count)
        , _log10combination_count(NEGATIVE_INFINITY)
        , _sorted_allele_counts(std::forward<Int32Vector>(sorted_allele_counts))
    {}

    GenotypeAlleleCounts(int32_t index, int32_t ploidy, const Int32Vector& sorted_allele_counts, int32_t distinct_allele_count,
                         pMemoryPool pool)
        : _index(index)
        , _ploidy(ploidy)
        , _distinct_allele_count(distinct_allele_count)
        , _log10combination_count(NEGATIVE_INFINITY)
        , _sorted_allele_counts(sorted_allele_counts, pool)
    {}

    /**
     * Implements binary search across allele indexes.
     * @param index the target index.
     * @param from first inclusive possible rank.
     * @param to last exclusive possible rank.
     * @return -1 or less if the allele index is not in the genotype false otherwise. You can obtain the potential insertion point
     * (within the interval [from,to]) as {@code -result - 1}
     */
    int32_t allele_index_to_rank(int32_t index, int32_t from, int32_t to) const;

    void copy_allele_counts_by_index(Int32Vector& dest, int32_t offset, int32_t minimum_allele_index, int32_t maximum_allele_index);
};

}  // namespace rovaca

#endif  // ROVACA_HC_GENOTYPE_ALLELE_COUNTS_H_
