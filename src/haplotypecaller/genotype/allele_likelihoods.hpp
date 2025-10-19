#ifndef ROVACA_HC_ALLELE_LIKELIHOODS_H_
#define ROVACA_HC_ALLELE_LIKELIHOODS_H_
#include <algorithm>
#include <functional>
#include <memory_resource>
#include <vector>

#include "allele.h"
#include "forward.h"
#include "genotype_macors.h"
#include "haplotype.h"
#include "indexed_allele_list.hpp"
#include "interface/interface_likelihood_matrix.hpp"
#include "interface/interface_sample_list.hpp"
#include "rovaca_logger.h"
#include "simple_interval.h"

namespace rovaca
{

template <typename EVIDENCE, typename A>
class AlleleLikelihoods : public InterfaceSampleList, public InterfaceAlleleList<A>
{
    /* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
public:
    static constexpr int32_t s_missing_index = -1;
    static constexpr double s_log_10_informative_threshold = 0.2000l;
    static constexpr double s_natural_log_informative_threshold = 0.4605170186l;

    typedef interfaceLikelihoodMatrix<EVIDENCE, A> LikelihoodMatrix;
    typedef typename LikelihoodMatrix::Te Te;  // TemplateEvidence
    typedef typename LikelihoodMatrix::Ta Ta;  // TemplateAllele
    typedef typename LikelihoodMatrix::TaVector TaVector;
    typedef typename LikelihoodMatrix::TeVector TeVector;
    typedef std::pmr::vector<TeVector> TeVector2D;
    typedef std::pmr::unordered_map<Te, int32_t> TeToInt32Map;
    typedef std::pmr::vector<TeToInt32Map> TeToInt32MapVector;
    typedef std::pmr::map<int32_t, TeVector> Int32ToTeVectorMap;
    typedef InterfaceAlleleList<A>* AlleleListPtr;
    typedef LikelihoodMatrix* LikelihoodMatrixPtr;
    typedef AlleleLikelihoods<Te, Ta>* AlleleLikelihoodsPtr;
    typedef std::pmr::vector<LikelihoodMatrixPtr> LikelihoodMatrixPtrVector;

    /* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
public:
    struct BestAllele
    {
        int32_t sample_index;
        Ta best_allele;
        Ta second_best_allele;
        Te evidence;
        double likelihood;
        double second_best_likelihood;
        double confidence;

        BestAllele(int32_t idx, Ta best, Ta second, Te evi, double best_likeli, double second_likeli)
            : sample_index(idx)
            , best_allele(best)
            , second_best_allele(second)
            , evidence(evi)
            , likelihood(best_likeli)
            , second_best_likelihood(second_likeli)
            , confidence(std::abs(best_likeli - second_likeli) < std::numeric_limits<double>::epsilon() ? 0.0 : best_likeli - second_likeli)
        {}

        bool is_informative() { return confidence > s_log_10_informative_threshold; }
    };

    typedef struct BestAllele* pBestAllele;

    /* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
private:
    class SampleMatrix : public LikelihoodMatrix
    {
        /* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
    private:
        size_t _sample_index;
        AlleleLikelihoodsPtr _self;

        /* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
    public:
        SampleMatrix(size_t si, AlleleLikelihoodsPtr self)
            : _sample_index(si)
            , _self(self)
        {}

        ~SampleMatrix() override = default;
        size_t number_of_alleles() const override { return _self->number_of_alleles(); }
        Ta get_allele(size_t index) const override { return _self->get_allele(index); }
        const TaVector& get_alleles() const override { return _self->get_alleles(); }
        int32_t index_of_allele(Ta allele) const override { return _self->index_of_allele(allele); }
        const TeVector& evidence() const override { return _self->sample_evidence(_sample_index); }
        size_t evidence_count() const override { return _self->sample_evidence_count(_sample_index); }
        int32_t index_of_evidence(Te e) const override { return _self->evidence_index(_sample_index, e); }
        Te get_evidence(size_t idx) const override { return _self->sample_evidence(_sample_index).at(idx); }
        double get(size_t a_idx, size_t e_idx) const override { return _self->sample_likelihoods(_sample_index).at(a_idx).at(e_idx); }
        void set(size_t a_idx, size_t e_idx, double value) override
        {
            _self->sample_likelihoods(_sample_index).at(a_idx).at(e_idx) = value;
        }
        void copy_allele_likelihoods(size_t allele_index, size_t offset, DoubleVector& dest) override
        {
            long evidence_count = (long)_self->sample_evidence_count(_sample_index);
            const DoubleVector& allele_likelihoods = _self->sample_likelihoods(_sample_index).at(allele_index);
            std::copy_n(allele_likelihoods.begin(), evidence_count, dest.begin() + (long)offset);
        }
    };

    /* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
private:
    pMemoryPool _pool;
    TeVector2D _evidence_by_sample_index;
    TeVector2D _filtered_evidence_by_sample_index;
    DoubleVector3D _values_by_sample_index;  // [s][a][r]
    Uint32Vector _likelihoods_matrix_evidence_capacity_by_sample_index;
    Uint32Vector _number_of_evidences;
    pInterfaceSampleList _samples;
    AlleleListPtr _alleles;
    bool _is_natural_log{false};
    pSimpleInterval _subsetted_genomic_loc{nullptr};
    TeToInt32MapVector _evidence_index_by_sample_index;
    int32_t _reference_allele_index;
    LikelihoodMatrixPtrVector _sample_matrices;

    /* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
public:
    template <typename e, typename a>
    static AlleleLikelihoods<e, a>* create(pMemoryPool pool, pInterfaceSampleList samples, AlleleListPtr alleles,
                                           const Int32ToTeVectorMap& evidence_by_sample, DoubleVector3D&& likelihoods)
    {
        return new (pool->allocate(sizeof(AlleleLikelihoods<e, a>)))
            AlleleLikelihoods<e, a>{pool, samples, alleles, evidence_by_sample, std::forward<DoubleVector3D>(likelihoods)};
    }

    template <typename e, typename a>
    static AlleleLikelihoods<e, a>* create(pMemoryPool pool, pInterfaceSampleList samples, AlleleListPtr alleles,
                                           const Int32ToTeVectorMap& evidence_by_sample)
    {
        return new (pool->allocate(sizeof(AlleleLikelihoods<e, a>))) AlleleLikelihoods<e, a>{pool, samples, alleles, evidence_by_sample};
    }

    template <typename e, typename a>
    static AlleleLikelihoods<e, a>* create(pMemoryPool pool, pInterfaceSampleList samples, AlleleListPtr alleles,
                                           TeVector2D&& evidence_by_sample_index, TeVector2D&& filtered_evidence_by_sample_index,
                                           DoubleVector3D&& values)
    {
        return new (pool->allocate(sizeof(AlleleLikelihoods<e, a>)))
            AlleleLikelihoods<e, a>{pool,
                                    samples,
                                    alleles,
                                    std::forward<TeVector2D>(evidence_by_sample_index),
                                    std::forward<TeVector2D>(filtered_evidence_by_sample_index),
                                    std::forward<DoubleVector3D>(values)};
    }

    ~AlleleLikelihoods() override = default;
    size_t number_of_alleles() const override { return _alleles->number_of_alleles(); }
    Ta get_allele(size_t index) const override { return _alleles->get_allele(index); }
    const TaVector& get_alleles() const override { return _alleles->get_alleles(); }
    int32_t index_of_allele(Ta allele) const override { return _alleles->index_of_allele(allele); }
    size_t number_of_samples() const override { return _samples->number_of_samples(); }
    const std::string& get_sample(size_t index) const override { return _samples->get_sample(index); }
    int32_t index_of_sample(const std::string& name) const override { return _samples->index_of_sample(name); }

    int32_t evidence_index(size_t si, Te evi)
    {
        const auto& evidence_at_sample = evidence_index_by_sample_index(si);
        if (evidence_at_sample.count(evi)) {
            return evidence_at_sample.at(evi);
        }
        return INVALID_INT;
    }
    uint32_t sample_evidence_count(size_t index) const { return _number_of_evidences.at(index); }
    const TeVector& sample_evidence(size_t index) const { return _evidence_by_sample_index.at(index); }
    DoubleVector2D& sample_likelihoods(size_t index) { return _values_by_sample_index.at(index); }

    bool is_natural_log() const { return _is_natural_log; }
    void set_is_natural_log(bool is_natural_log) { _is_natural_log = is_natural_log; }

    size_t evidence_count() const
    {
        size_t count = 0;
        std::for_each(_evidence_by_sample_index.begin(), _evidence_by_sample_index.end(), [&](const TeVector& v) { count += v.size(); });
        return count;
    }

    /*! @brief 从一个等位基因集到另一个（较小的一个）进行边缘化，取原始等位基因子集中每个证据的最大值 */
    template <typename B>
    AlleleLikelihoods<Te, B>* marginalize(const std::pmr::vector<B>& new_alleles,
                                          const std::pmr::map<B, std::pmr::vector<Ta>>& new_to_old_allele_map)
    {
        size_t new_allele_count = new_alleles.size();
        size_t sample_count = _samples->number_of_samples();

        Int32Vector old2new_idx_map = old_to_new_allele_index_map(new_alleles, new_allele_count, new_to_old_allele_map);

        DoubleVector3D new_likelihood_values = marginal_likelihoods_direct(new_allele_count, old2new_idx_map);

        TeVector2D new_evidence_by_sample_index(sample_count, _pool);
        for (size_t si = 0; si < sample_count; ++si) {
            new_evidence_by_sample_index.at(si) = _evidence_by_sample_index.at(si);
        }

        InterfaceAlleleList<B>* new_allele_list = IndexedAlleleList<B>::create(new_alleles, _pool);

        TeVector2D new_filtered_evidence_by_sample_index(_filtered_evidence_by_sample_index, _pool);

        auto* new_allelelikelihoosd = new (_pool->allocate(sizeof(AlleleLikelihoods<Te, B>)))
            AlleleLikelihoods<Te, B>{_pool,
                                     _samples,
                                     new_allele_list,
                                     std::move(new_evidence_by_sample_index),
                                     std::move(new_filtered_evidence_by_sample_index),
                                     std::move(new_likelihood_values)};
        new_allelelikelihoosd->set_is_natural_log(this->is_natural_log());
        return new_allelelikelihoosd;
    }

    /*!
     * @brief
     * @note this method modifies the current read-likelihoods collection
     * @param interval 根据 merged_vc 左右扩展一部分得到区间
     */
    void retain_evidence(pSimpleInterval interval)
    {
        size_t si, sample_count = _samples->number_of_samples();
        for (si = 0; si < sample_count; ++si) {
            const TeVector& sample_evidence = _evidence_by_sample_index.at(si);
            std::pmr::vector<size_t> remove_indices(_pool);
            for (size_t ei = 0, e_len = sample_evidence.size(); ei < e_len; ++ei) {
                if (!interval->overlaps(*sample_evidence.at(ei))) {
                    remove_indices.emplace_back(ei);
                }
            }
            if (!remove_indices.empty()) {
                remove_evidence_by_index(si, remove_indices);
            }

            TeVector new_sample_filtered(_pool);
            for (const auto& filtered_read : _filtered_evidence_by_sample_index.at(si)) {
                if (interval->overlaps(*filtered_read)) {
                    new_sample_filtered.emplace_back(filtered_read);
                }
            }
            if (!new_sample_filtered.empty()) {
                _filtered_evidence_by_sample_index.at(si) = std::move(new_sample_filtered);
            }
        }
    }

    void set_variant_calling_subset_used(pSimpleInterval loc) { _subsetted_genomic_loc = loc; }
    pSimpleInterval get_variant_calling_subset_applied() const { return _subsetted_genomic_loc; }

    /*! @brief 添加 non_ref_allele, 并更新对应数据 */
    void add_non_reference_allele()
    {
        pAllele non_ref_allele = StaticAllele::get_instance()->_non_ref_allele.get();
        if (_alleles->contains_allele(non_ref_allele)) {
            return;
        }

        AlleleVector new_alleles{_alleles->get_alleles(), _pool};
        new_alleles.push_back(StaticAllele::get_instance()->_non_ref_allele.get());
        _alleles = IndexedAlleleList<pAllele>::create(new_alleles, _pool);

        for (size_t si = 0, sample_count = _samples->number_of_samples(); si < sample_count; ++si) {
            size_t sample_evidence_count = _evidence_by_sample_index.at(si).size();
            uint32_t sample_evidence_capacity = _likelihoods_matrix_evidence_capacity_by_sample_index.at(si);
            DoubleVector new_allele_values(sample_evidence_capacity, NEGATIVE_INFINITY, _pool);
            std::fill(new_allele_values.begin() + (long)sample_evidence_count, new_allele_values.end(), NAN);

            _values_by_sample_index.at(si).push_back(std::move(new_allele_values));
        }
        update_non_ref_allele_likelihoods();
    }

    void update_non_ref_allele_likelihoods() { update_non_ref_allele_likelihoods(_alleles); }
    void update_non_ref_allele_likelihoods(AlleleListPtr alleles_to_consider)
    {
        int32_t non_ref_allele_index = index_of_allele(StaticAllele::get_instance()->_non_ref_allele.get());
        if (non_ref_allele_index == s_missing_index) {
            return;
        }

        pBestAllele best_allele;
        double allele_likelihood, non_ref_likelihood;
        size_t sample_count = _samples->number_of_samples();
        size_t allele_count = _alleles->number_of_alleles();
        size_t non_symbolic_allele_count = allele_count - 1, evidence_count;
        long number_of_qualified_allele_likelihoods;
        DoubleVector qualified_allele_likelihoods(non_symbolic_allele_count, _pool);
        for (size_t si = 0; si < sample_count; ++si) {
            std::fill(qualified_allele_likelihoods.begin(), qualified_allele_likelihoods.end(), 0.0);
            DoubleVector2D& sample_values = _values_by_sample_index.at(si);
            evidence_count = _evidence_by_sample_index.at(si).size();
            for (size_t r = 0; r < evidence_count; ++r) {
                best_allele = search_best_allele(si, r, true);
                number_of_qualified_allele_likelihoods = 0;
                for (size_t ai = 0; ai < allele_count; ++ai) {
                    allele_likelihood = sample_values.at(ai).at(r);
                    if (ai != (size_t)non_ref_allele_index && !std::isnan(allele_likelihood) &&
                        allele_likelihood < best_allele->likelihood &&
                        (alleles_to_consider == _alleles ||
                         alleles_to_consider->index_of_allele(_alleles->get_allele(ai)) != s_missing_index)) {
                        qualified_allele_likelihoods[number_of_qualified_allele_likelihoods++] = allele_likelihood;
                    }
                }
                non_ref_likelihood = evaluate(qualified_allele_likelihoods, number_of_qualified_allele_likelihoods);
                // when the median is NaN that means that all applicable likekihoods are the same as the best so the evidence is not
                // informative at all given the existing alleles. Unless there is only one (or zero) concrete alleles with give the same
                // (the best) likelihood to the NON-REF. When there is only one (or zero) concrete alleles we set the NON-REF likelihood to
                // NaN
                sample_values.at(non_ref_allele_index).at(r) = !std::isnan(non_ref_likelihood)  ? non_ref_likelihood
                                                               : non_symbolic_allele_count <= 1 ? NAN
                                                                                                : best_allele->likelihood;
            }
        }
    }

    LikelihoodMatrixPtr sample_matrix(size_t sample_index)
    {
        if (nullptr == _sample_matrices.at(sample_index)) {
            _sample_matrices.at(sample_index) = new ALLOC_TYPE_IN_POOL(_pool, SampleMatrix) SampleMatrix{sample_index, this};
        }
        return _sample_matrices.at(sample_index);
    }

    void add_evidence(const std::pmr::map<int32_t, std::pmr::vector<Te>>& evidence_by_sample, double initial_likelihood)
    {
        for (const auto& tup : evidence_by_sample) {
            int32_t sample_index = tup.first;
            const std::pmr::vector<Te>& extra_evidence = tup.second;
            CHECK_CONDITION_EXIT(_samples->get_sample(sample_index).empty(), "null sample");
            if (extra_evidence.empty()) {
                continue;
            }

            size_t old_count = _evidence_by_sample_index.at(sample_index).size();
            append_evidence(extra_evidence, sample_index);
            size_t new_count = _evidence_by_sample_index.at(sample_index).size();

            extends_likelihood_arrays(initial_likelihood, sample_index, old_count, new_count);
        }
    }

    std::pmr::vector<pBestAllele> best_alleles_breaking_ties()
    {
        std::pmr::vector<pBestAllele> result(_pool);
        auto lambda_func = [](void* a) -> double { return static_cast<pAllele>(a)->is_reference() ? 1.0 : 0.0; };

        for (int32_t s = 0, n = (int32_t)number_of_samples(); s < n; ++s) {
            std::pmr::vector<pBestAllele> best = best_alleles_breaking_ties(s, lambda_func);
            std::copy(best.begin(), best.end(), std::back_inserter(result));
        }
        return result;
    }

    std::pmr::vector<pBestAllele> best_alleles_breaking_ties(int32_t sample_index)
    {
        auto lambda_func = [](void* a) -> double { return static_cast<pAllele>(a)->is_reference() ? 1.0 : 0.0; };
        return best_alleles_breaking_ties(sample_index, lambda_func);
    }

    std::pmr::vector<pBestAllele> best_alleles_breaking_ties(int32_t sample_index, const std::function<double(void* a)>& func)
    {
        CHECK_CONDITION_EXIT(sample_index >= (int32_t)number_of_samples(), "sample_index >= number_of_samples()");

        const TaVector& alleles = _alleles->get_alleles();
        DoubleVector priorities{_pool};
        priorities.reserve(alleles.size());
        for (Ta aa : alleles) {
            priorities.emplace_back(func(aa));
        }

        size_t evidence_count = _evidence_by_sample_index.at(sample_index).size();
        std::pmr::vector<pBestAllele> result{_pool};
        result.reserve(evidence_count);
        for (size_t r = 0; r < evidence_count; ++r) {
            result.emplace_back(search_best_allele((size_t)sample_index, r, true, priorities));
        }

        return result;
    }

    static double evaluate(DoubleVector& values, long length)
    {
        if (length == 0) {
            return NAN;
        }
        if (length == 1) {
            return values[0];
        }

        std::sort(values.begin(), values.begin() + length);
        if (length & 1) {  // 奇数
            return values[int(length / 2)];
        }
        else {
            double rigth = values[int(length / 2)];
            double left = values[int(length / 2) - 1];
            return (rigth + left) / 2;
        }
    }

    /*! @brief 构造函数共有是不得已而为之，但请不要随意使用，尽量使用 create 方法 */
    AlleleLikelihoods(pMemoryPool pool, pInterfaceSampleList samples, AlleleListPtr alleles, TeVector2D&& evidence_by_sample_index,
                      TeVector2D&& filtered_evidence_by_sample_index, DoubleVector3D&& values)
        : _pool(pool)
        , _evidence_by_sample_index(std::move(evidence_by_sample_index), pool)
        , _filtered_evidence_by_sample_index(std::move(filtered_evidence_by_sample_index), pool)
        , _values_by_sample_index(std::move(values), pool)
        , _likelihoods_matrix_evidence_capacity_by_sample_index(pool)
        , _number_of_evidences(pool)
        , _samples(samples)
        , _alleles(alleles)
        , _evidence_index_by_sample_index(pool)
        , _reference_allele_index(alleles->index_of_reference())
        , _sample_matrices(pool)
    {
        size_t capacity, sample_count = samples->number_of_samples();

        _sample_matrices.resize(sample_count);
        _number_of_evidences.resize(sample_count);
        _evidence_index_by_sample_index.resize(sample_count);
        _likelihoods_matrix_evidence_capacity_by_sample_index.resize(sample_count);
        for (size_t si = 0; si < sample_count; ++si) {
            _number_of_evidences.at(si) = _evidence_by_sample_index.at(si).size();

            // we take the shortest allele's values array as the maximum evidence capacity for each sample
            capacity = INT32_MAX;
            const auto& sample_values = _values_by_sample_index.at(si);
            for (const auto& allele_values : sample_values) {
                capacity = std::min(capacity, allele_values.size());
            }
            _likelihoods_matrix_evidence_capacity_by_sample_index.at(si) = capacity;
        }
    }

    /*! @brief 构造函数共有是不得已而为之，但请不要随意使用，尽量使用 create 方法 */
    AlleleLikelihoods(pMemoryPool pool, pInterfaceSampleList samples, AlleleListPtr alleles, const Int32ToTeVectorMap& evidence_by_sample,
                      DoubleVector3D&& likelihoods)
        : _pool(pool)
        , _evidence_by_sample_index(pool)
        , _filtered_evidence_by_sample_index(pool)
        , _values_by_sample_index(std::move(likelihoods))
        , _likelihoods_matrix_evidence_capacity_by_sample_index(pool)
        , _number_of_evidences(pool)
        , _samples(samples)
        , _alleles(alleles)
        , _evidence_index_by_sample_index(pool)
        , _reference_allele_index(alleles->index_of_reference())
        , _sample_matrices(pool)
    {
        CHECK_CONDITION_EXIT(nullptr == samples, "samples is nullptr");
        CHECK_CONDITION_EXIT(nullptr == alleles, "alleles is nullptr");
        size_t sample_count = samples->number_of_samples();
        _evidence_by_sample_index.resize(sample_count);
        // _values_by_sample_index.resize(sample_count);
        _likelihoods_matrix_evidence_capacity_by_sample_index.resize(sample_count);
        _number_of_evidences.resize(sample_count);
        _evidence_index_by_sample_index.resize(sample_count);
        _filtered_evidence_by_sample_index.resize(sample_count);
        _sample_matrices.resize(sample_count);

        for (size_t si = 0; si < sample_count; ++si) {
            const auto& sample_evidences = evidence_by_sample.at(si);

            _number_of_evidences[si] = sample_evidences.empty() ? 0 : (uint32_t)sample_evidences.size();

            TeVector& evidence_at_sample = _evidence_by_sample_index.at(si);
            size_t evidence_at_sample_count = sample_evidences.size();
            evidence_at_sample.reserve(evidence_at_sample_count);
            std::for_each(sample_evidences.begin(), sample_evidences.end(),
                          [&evidence_at_sample](const Te& e) { evidence_at_sample.push_back(e); });

            _likelihoods_matrix_evidence_capacity_by_sample_index[si] = (uint32_t)evidence_at_sample_count;

            // _values_by_sample_index[si] = DoubleVector2D(allele_count, DoubleVector(evidence_at_sample_count, NEGATIVE_INFINITY), pool);
        }
    }

    /*! @brief 构造函数共有是不得已而为之，但请不要随意使用，尽量使用 create 方法 */
    AlleleLikelihoods(pMemoryPool pool, pInterfaceSampleList samples, AlleleListPtr alleles, const Int32ToTeVectorMap& evidence_by_sample)
        : _pool(pool)
        , _evidence_by_sample_index(pool)
        , _filtered_evidence_by_sample_index(pool)
        , _values_by_sample_index(pool)
        , _likelihoods_matrix_evidence_capacity_by_sample_index(pool)
        , _number_of_evidences(pool)
        , _samples(samples)
        , _alleles(alleles)
        , _evidence_index_by_sample_index(pool)
        , _reference_allele_index(alleles->index_of_reference())
        , _sample_matrices(pool)
    {
        CHECK_CONDITION_EXIT(nullptr == samples, "samples is nullptr");
        CHECK_CONDITION_EXIT(nullptr == alleles, "alleles is nullptr");

        size_t sample_count = samples->number_of_samples();
        size_t allele_count = alleles->number_of_alleles();

        _evidence_by_sample_index.resize(sample_count);
        _values_by_sample_index.resize(sample_count);
        _likelihoods_matrix_evidence_capacity_by_sample_index.resize(sample_count);
        _number_of_evidences.resize(sample_count);
        _evidence_index_by_sample_index.resize(sample_count);
        _filtered_evidence_by_sample_index.resize(sample_count);
        _sample_matrices.resize(sample_count);

        for (size_t si = 0; si < sample_count; ++si) {
            const auto& sample_evidences = evidence_by_sample.at(si);

            _number_of_evidences[si] = sample_evidences.empty() ? 0 : (uint32_t)sample_evidences.size();

            TeVector& evidence_at_sample = _evidence_by_sample_index.at(si);
            size_t evidence_at_sample_count = sample_evidences.size();
            evidence_at_sample.reserve(evidence_at_sample_count);
            std::for_each(sample_evidences.begin(), sample_evidences.end(),
                          [&evidence_at_sample](const Te& e) { evidence_at_sample.push_back(e); });

            _likelihoods_matrix_evidence_capacity_by_sample_index[si] = (uint32_t)evidence_at_sample_count;

            _values_by_sample_index[si] = DoubleVector2D(allele_count, DoubleVector(evidence_at_sample_count, NEGATIVE_INFINITY), pool);
        }
    }

    /* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
private:
    // Extends the likelihood arrays-matrices.
    void extends_likelihood_arrays(double initial_likelihood, int32_t sample_index, size_t old_count, size_t new_count)
    {
        DoubleVector2D& sample_values = _values_by_sample_index.at(sample_index);
        int32_t number_of_alleles = _alleles->number_of_alleles();
        ensure_likelihoods_matrix_evidence_capacity(sample_index, (uint32_t)new_count, sample_values, number_of_alleles);
        for (int32_t a = 0; a < number_of_alleles; ++a) {
            DoubleVector& values = sample_values.at(a);
            std::fill(values.begin() + (long)old_count, values.begin() + (long)new_count, initial_likelihood);
        }
    }

    // Resizes the lk value holding arrays to be able to handle at least "x" amount of evidence.
    void ensure_likelihoods_matrix_evidence_capacity(int32_t sample_index, uint32_t x, DoubleVector2D& sample_values,
                                                     int32_t number_of_alleles)
    {
        uint32_t current_capacity = _likelihoods_matrix_evidence_capacity_by_sample_index.at(sample_index);
        if (current_capacity < x) {
            uint32_t new_capacity = std::max(current_capacity, x) << 1;
            for (int32_t a = 0; a < number_of_alleles; ++a) {
                DoubleVector& values = sample_values.at(a);
                values.resize(new_capacity);
                std::fill(values.begin() + (long)current_capacity, values.end(), NAN);
            }
            _likelihoods_matrix_evidence_capacity_by_sample_index.at(sample_index) = new_capacity;
        }
    }

    void append_evidence(const std::pmr::vector<Te>& extra_evidence, int32_t sample_index)
    {
        std::pmr::vector<Te>& sample_evidence = _evidence_by_sample_index.at(sample_index);
        TeToInt32Map& sample_evidence_index = evidence_index_by_sample_index((size_t)sample_index);
        for (const Te& e : extra_evidence) {
            int32_t previous_value = int32_t(sample_evidence.size());
            if (!sample_evidence_index.count(e)) {
                sample_evidence.push_back(e);
                sample_evidence_index.insert({e, previous_value});
            }
            else {
                RovacaLogger::error("logic error");
                exit(EXIT_FAILURE);
            }
        }
        _number_of_evidences.at(sample_index) = (uint32_t)sample_evidence.size();
    }

    pBestAllele search_best_allele(size_t sample_index, size_t evidence_index, bool can_be_reference)
    {
        return search_best_allele(sample_index, evidence_index, can_be_reference, {});
    }

    /*!
     * @brief search the best allele for a unit of evidence
     * @param sample_index including sample index
     * @param evidence_index target evidence index
     * @param can_be_reference
     * @param priorities an array of allele priorities (higher values have higher priority) to be used, if present, to break ties for
     * uninformative likelihoods, in which case the evidence is assigned to the allele with the higher score
     * @return
     */
    pBestAllele search_best_allele(size_t sample_index, size_t evidence_index, bool can_be_reference, const DoubleVector& priorities)
    {
        //    private BestAllele(final int sampleIndex, final int evidenceIndex, final int bestAlleleIndex,
        //                   final double likelihood, final int secondBestAlleleIndex, final double secondBestLikelihood)

        size_t allele_count = _alleles->number_of_alleles();
        if (0 == allele_count || (1 == allele_count && _reference_allele_index == 0 && !can_be_reference)) {
            Te evidence = _evidence_by_sample_index.at(sample_index).at(evidence_index);
            return new ALLOC_TYPE_IN_POOL(_pool, BestAllele)
                BestAllele{(int32_t)sample_index, nullptr, nullptr, evidence, NEGATIVE_INFINITY, NEGATIVE_INFINITY};
        }

        const DoubleVector2D& sample_values = _values_by_sample_index.at(sample_index);
        int32_t best_allele_index = can_be_reference || _reference_allele_index != 0 ? 0 : 1;

        int32_t second_best_index = 0;
        double best_likelihood = sample_values.at(best_allele_index).at(evidence_index);
        double second_best_likelihood = NEGATIVE_INFINITY, candidate_likelihood;
        for (int32_t a = best_allele_index + 1; a < (int32_t)allele_count; a++) {
            if (!can_be_reference && _reference_allele_index == a) {
                continue;
            }
            candidate_likelihood = sample_values.at(a).at(evidence_index);
            if (candidate_likelihood > best_likelihood) {
                second_best_index = best_allele_index;
                best_allele_index = a;
                second_best_likelihood = best_likelihood;
                best_likelihood = candidate_likelihood;
            }
            else if (candidate_likelihood > second_best_likelihood) {
                second_best_index = a;
                second_best_likelihood = candidate_likelihood;
            }
        }

        if (!priorities.empty() && best_likelihood - second_best_likelihood < get_informative_threshold()) {
            double best_priority = priorities.at(best_allele_index);
            double second_best_priority = priorities.at(second_best_index);
            for (int a = 0; a < (int32_t)allele_count; a++) {
                candidate_likelihood = sample_values.at(a).at(evidence_index);
                if (a == best_allele_index || (!can_be_reference && a == _reference_allele_index) ||
                    best_likelihood - candidate_likelihood > get_informative_threshold()) {
                    continue;
                }
                double candidate_priority = priorities[a];

                if (candidate_priority > best_priority) {
                    second_best_index = best_allele_index;
                    best_allele_index = a;
                    second_best_priority = best_priority;
                    best_priority = candidate_priority;
                }
                else if (candidate_priority > second_best_priority) {
                    second_best_index = a;
                    second_best_priority = candidate_priority;
                }
            }
        }

        best_likelihood = sample_values.at(best_allele_index).at(evidence_index);
        second_best_likelihood =
            second_best_index != best_allele_index ? sample_values.at(second_best_index).at(evidence_index) : NEGATIVE_INFINITY;

        Ta best_allele = _alleles->get_allele(best_allele_index);
        Ta second_best_allele = _alleles->get_allele(second_best_index);
        Te evidence = _evidence_by_sample_index.at(sample_index).at(evidence_index);

        return new ALLOC_TYPE_IN_POOL(_pool, BestAllele)
            BestAllele{(int32_t)sample_index, best_allele, second_best_allele, evidence, best_likelihood, second_best_likelihood};
    }

    double get_informative_threshold() const
    {
        return _is_natural_log ? s_natural_log_informative_threshold : s_log_10_informative_threshold;
    }

    void remove_evidence_by_index(size_t sample_index, const std::pmr::vector<size_t>& evidences_to_remove)
    {
        size_t num_to_remove = evidences_to_remove.size();
        if (0 == num_to_remove) {
            return;
        }
        size_t old_evidence_count = _number_of_evidences.at(sample_index);
        size_t new_evidence_count = old_evidence_count - num_to_remove;

        const TeVector& old_evidence = _evidence_by_sample_index.at(sample_index);
        TeVector new_evidence(_pool);
        new_evidence.reserve(new_evidence_count);

        for (size_t i = 0, num_removed = 0; i < old_evidence_count; ++i) {
            if (num_removed < num_to_remove && i == evidences_to_remove.at(num_removed)) {
                ++num_removed;
            }
            else {
                new_evidence.emplace_back(old_evidence.at(i));
                for (auto& allele_values : _values_by_sample_index.at(sample_index)) {
                    allele_values.at(i - num_removed) = allele_values.at(i);
                }
            }
        }

        for (auto& allele_values : _values_by_sample_index.at(sample_index)) {
            std::fill(allele_values.begin() + (long)new_evidence_count, allele_values.end(), NAN);
        }

        _evidence_by_sample_index.at(sample_index) = std::move(new_evidence);
        _number_of_evidences.at(sample_index) = new_evidence_count;
        if (!_evidence_index_by_sample_index.empty()) {
            _evidence_index_by_sample_index.at(sample_index).clear();
        }
    }

    /*! @brief Calculate the marginal likelihoods considering the old -> new allele index mapping */
    DoubleVector3D marginal_likelihoods_direct(size_t new_allele_count, const Int32Vector& old2new_idx_map)
    {
        size_t sample_count = _samples->number_of_samples();
        size_t old_allele_count = old2new_idx_map.size();
        DoubleVector3D result(sample_count, _pool);

        double likelihoods;
        int32_t new_allele_idx;
        size_t si, ai, ri;
        for (si = 0; si < sample_count; ++si) {
            size_t sample_evidence_count = _evidence_by_sample_index.at(si).size();
            const DoubleVector2D& old_sample_values = _values_by_sample_index.at(si);
            DoubleVector2D new_sample_values(new_allele_count, DoubleVector(sample_evidence_count, NEGATIVE_INFINITY), _pool);
            for (ai = 0; ai < old_allele_count; ++ai) {
                new_allele_idx = old2new_idx_map.at(ai);
                if (new_allele_idx == s_missing_index) {
                    continue;
                }
                const DoubleVector& old_lk = old_sample_values.at(ai);
                DoubleVector& new_lk = new_sample_values.at(new_allele_idx);
                for (ri = 0; ri < sample_evidence_count; ++ri) {
                    likelihoods = old_lk.at(ri);
                    if (likelihoods > new_lk.at(ri)) {
                        new_lk[ri] = likelihoods;
                    }
                }
            }
            result[si] = std::move(new_sample_values);
        }
        return result;
    }

    template <typename B>
    Int32Vector old_to_new_allele_index_map(const std::pmr::vector<B>& new_alleles, size_t new_allele_count,
                                            const std::pmr::map<B, std::pmr::vector<Ta>>& new_to_old_allele_map)
    {
        size_t old_allele_count = _alleles->number_of_alleles();
        Int32Vector old2new_idx_map(old_allele_count, s_missing_index, _pool);

        int32_t old_allele_idx;
        for (size_t i = 0; i < new_allele_count; ++i) {
            const std::pmr::vector<Ta>& old_alleles_by_new_idx = new_to_old_allele_map.at(new_alleles.at(i));
            for (const Ta& old_allele : old_alleles_by_new_idx) {
                old_allele_idx = index_of_allele(old_allele);
                CHECK_CONDITION_EXIT(old_allele_idx == s_missing_index, "missing old allele in likelihood collection");
                CHECK_CONDITION_EXIT(old2new_idx_map.at(old_allele_idx) != s_missing_index,
                                     "two new alleles make reference to the same old allele");
                old2new_idx_map.at(old_allele_idx) = int32_t(i);
            }
        }
        return old2new_idx_map;
    }

    TeToInt32Map& evidence_index_by_sample_index(size_t si)
    {
        TeToInt32Map& cache = _evidence_index_by_sample_index.at(si);
        if (!cache.empty()) {
            return cache;
        }
        const TeVector& sample_evidence = _evidence_by_sample_index.at(si);
        size_t i, sample_evidence_count = sample_evidence.size();
        for (i = 0; i < sample_evidence_count; ++i) {
            cache.insert({sample_evidence.at(i), i});
        }
        return cache;
    }
};

}  // namespace rovaca

#endif  // ROVACA_HC_ALLELE_LIKELIHOODS_H_
