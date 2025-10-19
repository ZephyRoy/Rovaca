#ifndef BQSR_REPORTS_H
#define BQSR_REPORTS_H

#include <float.h>

#include <algorithm>
#include <cmath>
#include <map>
#include <string>
#include <vector>

#include "apply_bqsr_utils.h"
struct BQSRReport
{
    int ncols_, nrows_;
    std::string table_name_;
    std::string table_description_;
    std::vector<std::string> format_;
    std::vector<std::string> col_names_;
    std::vector<std::vector<std::string>> underlying_data_;
    std::map<std::string, int> col_names_to_index_;
    std::map<std::string, int> row_names_to_index_;

public:
    BQSRReport(int ncols, int nrows)
        : ncols_(ncols)
        , nrows_(nrows)
    {
        col_names_.resize(ncols);
        underlying_data_.resize(nrows, std::vector<std::string>(ncols, ""));
    }
    ~BQSRReport() {}
};

struct QuantizationInfo
{
    uint8_t* quantized_quals_;
    int64_t* empirical_qual_counts_;
    int quantizationLevels_;

public:
    QuantizationInfo() {}
    ~QuantizationInfo() {}
    inline int get_quantization_levels() { return quantizationLevels_; }
    inline const uint8_t* get_quantized_quals() { return quantized_quals_; }
    inline const int64_t* get_empirical_qual_counts() { return empirical_qual_counts_; }
    void non_quantization()
    {
        quantizationLevels_ = RECALDATUM_MAX_RECALIBRATED_Q_SCORE;
        for (int i = 0; i < quantizationLevels_; i++) {
            quantized_quals_[i] = static_cast<uint8_t>(i);
        }
    }
};

struct RecalDatum
{
    double estimated_qual;
    mutable double empirical_qual;
    double num_mismatches;
    long num_observations;
    bool valid;
    RecalDatum() { valid = false; }
    RecalDatum(long _num_observations, double numMismatches, double reportedQuality)
    {
        // TODO: spdlog error.
        if (_num_observations < 0) throw std::invalid_argument("Invalid num_observations");
        if (numMismatches < 0.0) throw std::invalid_argument("Invalid num_mismatches");
        if (reportedQuality < 0.0) throw std::invalid_argument("Invalid estimated_qual");
        num_observations = _num_observations;
        num_mismatches = numMismatches * RECALDATUM_MULTIPLIER;
        estimated_qual = reportedQuality;
        empirical_qual = RECALDATUM_UNINITIALIZED;
        // calculate_empirical_qual(estimated_qual);
        valid = true;
    }
    double get_empirical_qual(double conditional_prior) const
    {
        if (empirical_qual == RECALDATUM_UNINITIALIZED) {
            calculate_empirical_qual(conditional_prior);
        }
        return empirical_qual;
    }
    double get_num_observations() const { return num_observations; }
    double get_num_mismatches() const { return num_mismatches / RECALDATUM_MULTIPLIER; }

private:
    void calculate_empirical_qual(double conditional_prior) const
    {
        // smoothing is one error and one non-error observation
        long mismatches = static_cast<long>(get_num_mismatches() + 0.5) + RECALDATUM_SMOOTHING_CONSTANT;
        long observations = static_cast<long>(get_num_observations() + RECALDATUM_SMOOTHING_CONSTANT * 2);
        double new_empirical_qual = empirical_quality_bayesian_estimate(observations, mismatches, conditional_prior);
        empirical_qual = std::min(new_empirical_qual, (double)RECALDATUM_MAX_RECALIBRATED_Q_SCORE);
    }
};

struct RecalDatumArray
{
    RecalDatum* data;
    std::vector<int> dims;
    std::vector<int> indices;
    int current_dimention;
    void put_value(const RecalDatum value, const std::vector<int>& key)
    {
        if (key.size() != dims.size()) {
            throw std::invalid_argument("Invalid value dimensions");
        }

        int index = 0;
        for (size_t i = 0; i < key.size(); i++) {
            if (key[i] < 0 || key[i] >= dims[i]) {
                throw std::out_of_range("Value out of range");
            }
            index += key[i] * indices[i];
        }
        // 在对应的索引位置存储值
        data[index] = value;
    }

    const RecalDatum& get_value(const std::vector<int>& key) const
    {
        if (key.size() != dims.size()) {
            throw std::invalid_argument("Invalid key dimensions");
        }

        int index = 0;
        for (size_t i = 0; i < key.size(); i++) {
            if (key[i] < 0 || key[i] >= dims[i]) {
                throw std::out_of_range("Key out of range");
            }
            index += key[i] * indices[i];
        }
        return data[index];
    }
    const RecalDatum& get_value(int key0, int key1) const
    {
        int index = key0 * indices[0] + key1 * indices[1];
        return data[index];
    }
    const RecalDatum& get_value(int key0, int key1, int key2) const
    {
        int index = key0 * indices[0] + key1 * indices[1] + key2 * indices[2];
        return data[index];
    }
    const RecalDatum& get_value(int key0, int key1, int key2, int key3) const
    {
        int index = key0 * indices[0] + key1 * indices[1] + key2 * indices[2] + key3 * indices[3];
        return data[index];
    }

    RecalDatumArray() {}

    ~RecalDatumArray()
    {
        if (data) delete[] data;
    }

    void create(const std::vector<int>& dimention)
    {
        size_t dims_size = dimention.size();
        if (!dims_size) return;
        dims = dimention;
        indices.resize(dims.size());
        int product = 1;
        for (int i = dims.size() - 1; i >= 0; i--) {
            indices[i] = product;
            product *= dims[i];
        }
        // 分配内存
        int total_elements = product;
        if (!total_elements) return;
        data = new RecalDatum[total_elements];
    }
};

struct RecalibrationTables
{
    RecalDatumArray read_group_tables_;
    RecalDatumArray quality_score_tables_;
    RecalDatumArray context_tables_;
    RecalDatumArray cycle_tables_;

public:
    RecalibrationTables() {}
    ~RecalibrationTables() {}
    inline const RecalDatumArray& get_read_group_table() { return read_group_tables_; }
    inline const RecalDatumArray& get_quality_score_table() { return quality_score_tables_; }
    inline const RecalDatumArray& get_context_table() { return context_tables_; }
    inline const RecalDatumArray& get_cycle_table() { return cycle_tables_; }
    RecalDatumArray& get_table_by_name(const std::string& name)
    {
        if (name == "Cycle") {
            return cycle_tables_;
        }
        else
            return context_tables_;
    }
};

#endif  // BQSR_REPORTS_H