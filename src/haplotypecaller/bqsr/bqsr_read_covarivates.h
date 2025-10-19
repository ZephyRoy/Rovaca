#ifndef BQSR_READ_COVARIVATES_H
#define BQSR_READ_COVARIVATES_H

#include <htslib/sam.h>

#include <functional>
#include <map>
#include <string>
#include <vector>

#include <stdexcept>

#include "apply_bqsr_utils.h"
struct ReadCovariates
{
    // int*** keys;
    // int keys[EVENT_TYPE_MAX][BQSR_MAX_QUAL_LENGTH][COVARIATE_TABLE_COUNT];
    int keys[BQSR_MAX_QUAL_LENGTH][COVARIATE_TABLE_COUNT];
    int current_covariate_index = 0;
    bool negative_strand;
    uint32_t read_len;
    bam1_t* read;
    uint32_t read_len_after_clip;
    uint8_t* read_base_quality;
    uint8_t* read_insert_quality;
    uint8_t* read_deletion_quality;
    uint8_t* stranded_clipped_bases;
    std::vector<int> mismatch_keys;
    void record_rg_covariate_value();
    void record_qual_score_covariate_value();
    void record_context_covariate_value();
    void record_cycle_covariate_value();
    void record_all_covariate_values();
    void add_covariate(int mismatch, int insertion, int deletion, int read_offset);
    void process_covariates(const std::function<int(uint32_t)>& covariate_func, int index, bool need_clipping = false);
    ReadCovariates(bam1_t* read)
    {
        read_insert_quality = nullptr;
        read_deletion_quality = nullptr;
        read_len = read->core.l_qseq;
        stranded_clipped_bases = new uint8_t[read_len];
        mismatch_keys.resize(read_len);
    }
    ~ReadCovariates() { delete[] stranded_clipped_bases; }
};
struct ReadGroupCovariate
{
    std::map<std::string, int32_t> rg_table_;
    std::map<int32_t, std::string> rg_reverse_table_;
};
struct QualityScoreCovariate
{
    int key_from_value(const std::string& value)
    {
        int i = parseInt(value, 10);
        if (i >= -128 && i <= 127) return i;
        return INT32_MAX;
    }
};
struct ContextCovariate
{
    int mismatchesContextSize;
    int indelsContextSize;
    int mismatchesKeyMask;
    int indelsKeyMask;

    int LENGTH_MASK = 15;
    // int MAX_DNA_CONTEXT = 13;
    uint8_t low_qual_tail;
    int max_key_value()
    {
        int length = std::max(mismatchesContextSize, indelsContextSize);
        int key = length;
        int bit_offset = LENGTH_BITS;
        for (int i = 0; i < length; ++i) {
            key |= (3 << bit_offset);
            bit_offset += 2;
        }
        return key;
    }
    int key_from_value(const std::string& value) { return bqsr_covariate_key_from_context((uint8_t*)value.data(), 0, value.length()); }
};
struct CycleCovariates
{
    int MAXIMUM_CYCLE_VALUE;
    int CUSHION_FOR_INDELS = 4;
    int max_key_value() { return (MAXIMUM_CYCLE_VALUE << 1) + 1; }
    int key_from_value(const std::string& value)
    {
        int i = parseInt(value, 10);
        int result = std::abs(i);
        if (result > MAXIMUM_CYCLE_VALUE) return INT32_MAX;
        result <<= 1;
        if (i < 0) result++;
        return result;
    }
};
struct StandardCovariates
{
    ReadGroupCovariate rg_covariate_;
    QualityScoreCovariate quality_score_covariate_;
    ContextCovariate context_covariates_;
    CycleCovariates cycle_covariates_;
};

struct RecalArguments
{
    bool DO_NOT_USE_STANDARD_COVARIATES = false;
    std::string SOLID_RECAL_MODE = "SET_Q_ZERO";
    std::string SOLID_NOCALL_STRATEGY = "THROW_EXCEPTION";
    bool UN_WITHOUT_DBSNP = false;
    int MISMATCHES_CONTEXT_SIZE = 2;
    int INDELS_CONTEXT_SIZE = 3;
    int MAXIMUM_CYCLE_VALUE = 500;
    uint8_t MISMATCHES_DEFAULT_QUALITY = -1;
    uint8_t INSERTIONS_DEFAULT_QUALITY = 45;
    uint8_t DELETIONS_DEFAULT_QUALITY = 45;
    int LOW_QUAL_TAIL = 2;
    double BAQGOP = 40;
    int QUANTIZING_LEVELS = 16;
    bool enableBAQ = false;
    bool computeIndelBQSRTables = false;
    bool useOriginalBaseQualities = false;
    uint8_t defaultBaseQualities = -1;
};
int get_stranded_clipped_bases(bam1_t* read, uint8_t* clipped_base, int low_qual_tail);
void context_with(std::vector<int>& mismatch_keys, int read_len, uint8_t* bases, int context_size, int mask);
int cycle_key(int base_number, bam1_t* read, bool indel_in_count, int max_cycle);
int key_from_cycle(int cycle, int max_cycle);

#endif  // BQSR_READ_COVARIVATES_H