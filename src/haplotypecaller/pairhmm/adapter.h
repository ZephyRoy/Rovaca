#ifndef __ADAPTER_H__
#define __ADAPTER_H__
#include <cstdint>

constexpr int MIN_QUALITY = 6;
constexpr int MAX_REPEAT_LENGTH = 20;
constexpr int MAX_STR_UNIT_LENGTH = 8;
constexpr int MIN_ADJUSTED_QSCORE = 10;
constexpr int MIN_QUALITY_THRESHOLD = 18;
constexpr int DYNAMIC_QUAL_THRESHOLD_TABLE_ENTRY_LENGTH = 3;
constexpr int MAXIMUM_DYNAMIC_QUAL_THRESHOLD_ENTRY_BASEQ = 40;
constexpr int DYNAMIC_QUAL_THRESHOLD_TABLE_ENTRY_MEAN_OFFSET = 1;
constexpr int PCR_INDEL_MODEL_CACHE_HOSTILE[21] = {40, 40, 39, 38, 37, 36, 34, 32, 28, 23, 17, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10};
constexpr int PCR_INDEL_MODEL_CACHE_AGGRESSIVE[21] = {40, 40, 40, 39, 39, 39, 38, 38, 37, 37, 36, 35, 34, 33, 32, 30, 28, 26, 23, 20, 17};
constexpr int PCR_INDEL_MODEL_CACHE_CONSERVATIVE[21] = {40, 40, 40, 40, 39, 39, 39, 39, 39, 38, 38, 38, 37, 37, 37, 36, 36, 35, 34, 33, 33};

constexpr double ROUND_UP = 0.5;
constexpr double INITIAL_QSCORE = 40.0;
constexpr double TRISTATE_CORRECTION = 3.0;
constexpr double LOG10_QUALITY_PER_BASE = -4.0;
constexpr double EXPECTED_ERROR_RATE_PER_BASE = 0.02;
constexpr double MAXIMUM_EXPECTED_ERROR_PER_READ = 2.0;
constexpr double MAXIMUM_BEST_ALT_LIKELIHOOD_DIFFERENCE = -4.5;

constexpr float LOG10_INITIAL_CONSTANT_F = 36.1236000061;
constexpr double LOG10_INITIAL_CONSTANT_D = 307.050595577260822;

constexpr bool CAP_MIN_LIKELIHOOD = true;
constexpr bool DYNAMIC_DISQUALIFICATION = false;

constexpr size_t k_alignment = 64;
constexpr double k_length_ratio = 1.2;
constexpr uint32_t k_min_diff_step = 15;
constexpr uint32_t k_sum_count_len = 5;
constexpr uint32_t k_read_group_count = 8;
constexpr uint32_t k_avx512_float_hap_concurrency = 16;
constexpr uint32_t k_avx512_float_read_concurrency = 16;

#endif  // __ADAPTER_H__