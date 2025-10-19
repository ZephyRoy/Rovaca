#include "pairhmm_internal.h"

#include <memory_resource>

#include "../common/utils/rovaca_memory_pool.h"
#include "../genotype/haplotype.h"
#include "../genotype/read_record.h"
#include "adapter.h"
#include "core/avx.h"
#include "core/avx_impl.h"
#include "rovaca/avx_512_float.h"
#include "rovaca/common.h"
#include "rovaca/context.h"

extern float (*s_compute_full_prob_float_old_pairhmm)(const TestCase& tc);
extern double (*s_compute_full_prob_double_old_pairhmm)(const TestCase& tc);

namespace rovaca
{

inline void init_native()
{
    // enable FTZ
    _MM_SET_FLUSH_ZERO_MODE(_MM_FLUSH_ZERO_ON);
}

/*
 * 'A' : 0x1
 * 'C' : 0x2
 * 'T' : 0x4
 * 'G' : 0x8
 * 'N' : 0xF
 */
static uint32_t convert_arr[256] = {
    0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,   0, 0,   0, 0, 0, 0,   0, 0,
    0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0x1, 0, 0x2, 0, 0, 0, 0x8, 0, 0,
    0, 0, 0, 0, 0xF, 0, 0, 0, 0, 0, 0x4, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,   0, 0,   0, 0, 0, 0,   0, 0,
    0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,   0, 0,   0, 0, 0, 0,   0, 0,
    0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,   0, 0,   0, 0, 0, 0,   0, 0,
    0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,   0, 0,   0, 0, 0, 0,   0, 0,
    0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,   0, 0,   0, 0, 0};

bool string_equal(const uint8_t* str1, const uint8_t* str2, int32_t len)
{
    for (int32_t i{0}; i < len; ++i) {
        if (str1[i] != str2[i]) {
            return false;
        }
    }
    return true;
}

void normalize_likelihoods(DoubleVector2D& log_likelihoods)
{
    double best, cmp;
    for (DoubleVector& likelihoods : log_likelihoods) {
        best = std::max_element(likelihoods.begin(), likelihoods.end()).operator*();
        cmp = best + MAXIMUM_BEST_ALT_LIKELIHOOD_DIFFERENCE;
        for (double& val : likelihoods) {
            val = val < cmp ? cmp : val;
        }
    }
}

int32_t find_tandem_repeat_units(pBases bases, uint32_t offset1)
{
    int32_t offset = static_cast<int32_t>(offset1);
    int32_t base_num = static_cast<int32_t>(bases->num);

    int32_t max_bw = 0;
    std::pair<int32_t, int32_t> best_bw_repeat_unit(offset, 1);
    for (int32_t str{1}; str <= MAX_STR_UNIT_LENGTH; str++) {
        if (offset + 1 - str < 0) {
            break;
        }
        max_bw = find_number_of_repetitions(bases, offset - str + 1, str, bases, 0, offset + 1, false);
        if (max_bw > 1) {
            best_bw_repeat_unit.first = offset - str - 1;
            best_bw_repeat_unit.second = str;
            break;
        }
    }

    int32_t max_rl = max_bw;
    std::pair<int32_t, int32_t> best_repeat_unit = best_bw_repeat_unit;

    if (offset < base_num - 1) {
        std::pair<int32_t, int32_t> best_FW_repeat_unit(offset + 1, 1);
        int32_t max_fw = 0;

        for (int32_t str{1}; str <= MAX_STR_UNIT_LENGTH; str++) {
            if (offset + str + 1 > base_num) {
                break;
            }

            max_fw = find_number_of_repetitions(bases, offset + 1, str, bases, offset + 1, bases->num - offset - 1, true);
            if (max_fw > 1) {
                best_FW_repeat_unit.first = offset + 1;
                best_FW_repeat_unit.second = str;
                break;
            }
        }
        if (best_FW_repeat_unit == best_bw_repeat_unit) {
            max_rl = max_bw + max_fw;
            best_repeat_unit = best_FW_repeat_unit;
        }
        else {
            max_bw = find_number_of_repetitions(bases, best_FW_repeat_unit.first, best_FW_repeat_unit.second, bases, 0, offset + 1, false);
            max_rl = max_bw + max_fw;
            best_repeat_unit = best_FW_repeat_unit;
        }
    }

    if (max_rl > MAX_REPEAT_LENGTH) {
        max_rl = MAX_REPEAT_LENGTH;
    }

    return max_rl;
}

void apply_pcr_error_model(pBases bases, pBases gap_qual, PcrIndelModel pcr_option)
{
    int32_t repeat_length = 0;
    switch (pcr_option) {
        case PcrIndelModel::HOSTILE: {
            for (uint32_t i = 1; i < bases->num; i++) {
                repeat_length = find_tandem_repeat_units(bases, i - 1);
                gap_qual->data[i - 1] = gap_qual->data[i - 1] < PCR_INDEL_MODEL_CACHE_HOSTILE[repeat_length]
                                            ? gap_qual->data[i - 1]
                                            : PCR_INDEL_MODEL_CACHE_HOSTILE[repeat_length];
            }
            break;
        }

        case PcrIndelModel::AGGRESSIVE: {
            for (uint32_t i = 1; i < bases->num; i++) {
                repeat_length = find_tandem_repeat_units(bases, i - 1);
                gap_qual->data[i - 1] = gap_qual->data[i - 1] < PCR_INDEL_MODEL_CACHE_AGGRESSIVE[repeat_length]
                                            ? gap_qual->data[i - 1]
                                            : PCR_INDEL_MODEL_CACHE_AGGRESSIVE[repeat_length];
            }
            break;
        }

        case PcrIndelModel::CONSERVATIVE: {
            for (uint32_t i = 1; i < bases->num; i++) {
                repeat_length = find_tandem_repeat_units(bases, i - 1);
                gap_qual->data[i - 1] = gap_qual->data[i - 1] < PCR_INDEL_MODEL_CACHE_CONSERVATIVE[repeat_length]
                                            ? gap_qual->data[i - 1]
                                            : PCR_INDEL_MODEL_CACHE_CONSERVATIVE[repeat_length];
            }
            break;
        }

        default: {
            break;
        }
    }
}

void filter_poorly_modelled_evidence(ReadVector& reads, DoubleVector2D& likelihoods, pMemoryPool pool)
{
    pReadRecord read_i;
    double best_likelihoods, log10_min_likelihoods;
    size_t number_of_evidence = reads.size();
    std::pmr::vector<size_t> remove_idx{pool};
    for (size_t i{0}; i < number_of_evidence; ++i) {
        read_i = reads[i];
        const auto& likelihoods_i = likelihoods[i];
        best_likelihoods = std::max_element(likelihoods_i.begin(), likelihoods_i.end()).operator*();
        log10_min_likelihoods = std::min(2.0, std::ceil(read_i->seq_length() * EXPECTED_ERROR_RATE_PER_BASE)) * LOG10_QUALITY_PER_BASE;

        if (best_likelihoods < log10_min_likelihoods) {
            remove_idx.push_back(i);
        }
    }

    for (auto it = remove_idx.rbegin(); it != remove_idx.rend(); ++it) {
        reads.erase(reads.begin() + *it);
        likelihoods.erase(likelihoods.begin() + *it);
    }
}

void transpose_likelihood_matrix(const DoubleVector2D& source, DoubleVector2D& target)
{
    size_t raw = source.size(), col = target.size();
    for (size_t i{0}; i < raw; ++i) {
        const DoubleVector& source_i = source[i];
        for (size_t j{0}; j < col; ++j) {
            target[j][i] = source_i[j];
        }
    }
}

int32_t find_number_of_repetitions(pBases repeat_unit_full, int32_t offset_in_repeat_unit_full, int32_t repeat_unit_length,
                                   pBases test_string_full, int32_t offset_in_test_string_full, int32_t test_string_length,
                                   bool leading_repeats)
{
    if (test_string_length == 0) {
        return 0;
    }

    int32_t num_repeats = 0;
    int32_t length_diff = test_string_length - repeat_unit_length;

    if (leading_repeats) {
        for (int32_t start{0}; start <= length_diff; start += repeat_unit_length) {
            if (string_equal(test_string_full->data + start + offset_in_test_string_full,
                             repeat_unit_full->data + offset_in_repeat_unit_full, repeat_unit_length)) {
                ++num_repeats;
            }
            else {
                return num_repeats;
            }
        }
    }
    else {
        for (int32_t start = length_diff; start >= 0; start -= repeat_unit_length) {
            if (string_equal(test_string_full->data + start + offset_in_test_string_full,
                             repeat_unit_full->data + offset_in_repeat_unit_full, repeat_unit_length)) {
                ++num_repeats;
            }
            else {
                return num_repeats;
            }
        }
    }
    return num_repeats;
}

std::pmr::vector<ReadLenInfo> generate_histogram(const ReadVector& rs, uint32_t read_num, uint32_t max_read_len, pMemoryPool pool)
{
    std::pmr::vector<ReadLenInfo> read_len_histogram{max_read_len + 1, {0, 0}, pool};

    uint32_t read_len;
    for (uint32_t i{0}; i < read_num; ++i) {
        read_len = rs.at(i)->seq_length();
        ReadLenInfo& info = read_len_histogram.at(read_len);
        if (0 == info.count) {
            info.idx = i;
        }
        ++info.count;
    }

    return read_len_histogram;
}

void accumulation(ReadLenInfo& info, ReadGroup& acc, uint32_t& read_count, uint32_t& idx)
{
    uint32_t count = info.count;
    for (uint32_t k{0}; k < count; ++k) {
        acc.idxs[idx++] = info.idx + k;
    }
    read_count += count;
    info.count = 0;
}

void debug_print(uint32_t* arr, uint32_t len = k_avx512_float_read_concurrency)
{
    for (uint32_t i{0}; i < len; ++i) {
        printf("%d\t", arr[i]);
    }

    printf("\n");
    fflush(stdout);
}

std::pmr::vector<ReadGroup> generate_new_pairhmm_idx(std::pmr::vector<ReadLenInfo>& histogram, pMemoryPool pool)
{
    ReadGroup temp{};
    std::pmr::vector<ReadGroup> result{pool};
    uint32_t max_read_len = histogram.size();

    // First calculate reads with count ≥16 and same length
    for (uint32_t i{0}; i < max_read_len; ++i) {
        ReadLenInfo& info = histogram.at(i);
        while (info.count >= k_avx512_float_read_concurrency) {
            for (uint32_t k{0}; k < k_avx512_float_read_concurrency; ++k) {
                temp.idxs[k] = info.idx + k;
            }
            result.push_back(temp);

            info.idx += k_avx512_float_read_concurrency;
            info.count -= k_avx512_float_read_concurrency;
        }
    }

    uint32_t count;
    uint32_t read_count{0};
    uint32_t sub1, sub2, sub3, sub4, sub5;
    uint32_t add1, add2, add3, add4, add5;
    // Calculate reads with count ≥k_read_group_count and same length, and group with similar length reads into groups of 16
    for (uint32_t i{max_read_len - 1}; i < max_read_len; --i) {
        ReadLenInfo& info_out = histogram.at(i);
        if (info_out.count < k_read_group_count) {
            continue;
        }

        read_count = 0;

    // Prefer to find shorter ones first
        sub1 = i >= 1 ? histogram.at(i - 1).count : 0;
        sub2 = i >= 2 ? histogram.at(i - 2).count : 0;
        sub3 = i >= 3 ? histogram.at(i - 3).count : 0;
        sub4 = i >= 4 ? histogram.at(i - 4).count : 0;
        sub5 = i >= 5 ? histogram.at(i - 5).count : 0;

    // The sum on the left meets 16 reads, no need to verify the index below, it must exist
        if (sub1 + sub2 + sub3 + sub4 + sub5 + info_out.count >= k_avx512_float_read_concurrency) {
            for (uint32_t k{0}; k <= k_sum_count_len; ++k) {
                ReadLenInfo& info = histogram.at(i - k);
                count = info.count;

                // Greater than or equal to 16 reads
                if (read_count + count >= k_avx512_float_read_concurrency) {
                    uint32_t required_count = k_avx512_float_read_concurrency - read_count;
                    uint32_t offset = read_count + count - k_avx512_float_read_concurrency;
                    for (uint32_t j = 0; j < required_count; ++j) {
                        temp.idxs[j] = info.idx + j + offset;
                    }
                    info.count -= required_count;
                    break;
                }

                read_count += count;
                for (uint32_t j = 0, idx = k_avx512_float_read_concurrency - read_count; j < count; ++j, ++idx) {
                    temp.idxs[idx] = info.idx + j;
                }
                info.count = 0;
            }

            result.push_back(temp);
            continue;
        }

        add1 = i < (max_read_len - 1) ? histogram.at(i + 1).count : 0;
        add2 = i < (max_read_len - 2) ? histogram.at(i + 2).count : 0;
        add3 = i < (max_read_len - 3) ? histogram.at(i + 3).count : 0;
        add4 = i < (max_read_len - 4) ? histogram.at(i + 4).count : 0;
        add5 = i < (max_read_len - 5) ? histogram.at(i + 5).count : 0;

    // The sum of left and right meets 16 reads, no need to verify the index below, it must exist
        if (sub1 + sub2 + sub3 + sub4 + sub5 + add1 + add2 + add3 + add4 + add5 + info_out.count >= k_avx512_float_read_concurrency) {
            uint32_t idx = 0;

            if (sub1) {
                accumulation(histogram.at(i - 1), temp, read_count, idx);
            }
            if (sub2) {
                accumulation(histogram.at(i - 2), temp, read_count, idx);
            }
            if (sub3) {
                accumulation(histogram.at(i - 3), temp, read_count, idx);
            }
            if (sub4) {
                accumulation(histogram.at(i - 4), temp, read_count, idx);
            }
            if (sub5) {
                accumulation(histogram.at(i - 5), temp, read_count, idx);
            }

            accumulation(histogram.at(i), temp, read_count, idx);

            for (uint32_t k{1}; k <= k_sum_count_len; ++k) {
                ReadLenInfo& info = histogram.at(i + k);
                count = info.count;

                // Greater than or equal to 16 reads
                if (read_count + count >= k_avx512_float_read_concurrency) {
                    uint32_t required_count = k_avx512_float_read_concurrency - read_count;
                    for (uint32_t j = 0; j < required_count; ++j) {
                        temp.idxs[idx++] = info.idx + j;
                    }
                    info.idx += required_count;
                    info.count -= required_count;
                    break;
                }

                read_count += count;
                for (uint32_t j{0}; j < count; ++j) {
                    temp.idxs[idx++] = info.idx + j;
                }
                info.count = 0;
            }

            result.push_back(temp);
        }
    }

    return result;
}

void init_test_case(TestCase& tc, uint32_t max_read_len, uint32_t max_hap_len, pMemoryPool pool)
{
    tc.hap_num = k_avx512_float_hap_concurrency;
    tc.read_num = k_avx512_float_read_concurrency;
    tc.mm = (float*)pool->allocate(sizeof(int32_t) * k_avx512_float_read_concurrency * max_read_len, k_alignment);
    tc.mi = (float*)pool->allocate(sizeof(int32_t) * k_avx512_float_read_concurrency * max_read_len, k_alignment);
    tc.ii = (float*)pool->allocate(sizeof(int32_t) * k_avx512_float_read_concurrency * max_read_len, k_alignment);
    tc.md = (float*)pool->allocate(sizeof(int32_t) * k_avx512_float_read_concurrency * max_read_len, k_alignment);
    tc.dd = (float*)pool->allocate(sizeof(int32_t) * k_avx512_float_read_concurrency * max_read_len, k_alignment);
    tc.gapm = (float*)pool->allocate(sizeof(int32_t) * k_avx512_float_read_concurrency * max_read_len, k_alignment);
    tc.distm = (float*)pool->allocate(sizeof(int32_t) * k_avx512_float_read_concurrency * max_read_len, k_alignment);
    tc._1_distm = (float*)pool->allocate(sizeof(int32_t) * k_avx512_float_read_concurrency * max_read_len, k_alignment);

    tc.hap_len_arr = (uint32_t*)pool->allocate(sizeof(uint32_t) * k_avx512_float_read_concurrency, k_alignment);
    tc.read_len_arr = (uint32_t*)pool->allocate(sizeof(uint32_t) * k_avx512_float_read_concurrency, k_alignment);

    tc.haps = (int32_t*)pool->allocate(sizeof(int32_t) * max_hap_len, k_alignment);
    tc.reads = (int32_t*)pool->allocate(sizeof(int32_t) * tc.read_num * max_read_len, k_alignment);
}

void decompression_data(ReadVector& rs, TestCase& tc, rovaca::pBases* read_base, rovaca::pBases* read_qual, rovaca::pBases* ins_gops,
                        rovaca::pBases* del_gops, rovaca::pBases* gap_conts, const uint32_t* read_idxs, PcrIndelModel pcr_option,
                        int32_t min_quality_threshold)
{
    uint32_t read_idx, read_len;
    uint32_t min_len = UINT32_MAX, max_len = 0;
    for (uint32_t k{0}; k < k_avx512_float_read_concurrency; ++k) {
        read_idx = read_idxs[k];

        rs[read_idx]->decode_base(read_base[k]);
        rs[read_idx]->decode_qual(read_qual[k]);
        rs[read_idx]->ins_gops(ins_gops[k]);
        del_gops[k] = ins_gops[k];
        // rs[read_idx]->del_gops(del_gops[k]);
        rs[read_idx]->gap_conts(gap_conts[k]);

        read_len = read_base[k]->num;
        tc.read_len_arr[k] = read_len;
        min_len = std::min(min_len, read_len);
        max_len = std::max(max_len, read_len);

        // modify_read_qualities
        uint8_t mq = rs[read_idx]->mapping_quality();
        for (uint32_t x = 0; x < read_qual[k]->num; ++x) {
            read_qual[k]->data[x] = std::min(read_qual[k]->data[x], mq);
            read_qual[k]->data[x] = read_qual[k]->data[x] < min_quality_threshold ? MIN_QUALITY : read_qual[k]->data[x];
        }

        if (pcr_option != PcrIndelModel::NONE) {
            apply_pcr_error_model(read_base[k], ins_gops[k], pcr_option);
        }
    }

    tc.max_read_len = max_len;
    tc.min_read_len = min_len;
}

void prepare_data(TestCase& tc, rovaca::Context<float>& ctx, rovaca::pBases* read_base, rovaca::pBases* read_qual, rovaca::pBases* ins_gops,
                  rovaca::pBases* del_gops, rovaca::pBases* gap_conts)
{
    rovaca::pBases b_, q_, i_, d_, c_;
    int _i, _d, _c, _q;

    for (uint32_t i{0}; i < k_avx512_float_read_concurrency; ++i) {
        b_ = read_base[i];
        q_ = read_qual[i];
        i_ = ins_gops[i];
        d_ = del_gops[i];
        c_ = gap_conts[i];

        for (uint32_t j{0}; j < b_->num; ++j) {
            tc.reads[i + j * k_avx512_float_read_concurrency] = convert_arr[read_base[i]->data[j]];

            _i = i_->data[j] & 127;
            _d = d_->data[j] & 127;
            _c = c_->data[j] & 127;
            _q = q_->data[j] & 127;

            tc.mm[i + j * k_avx512_float_read_concurrency] = ctx.set_mm_prob(_i, _d);
            tc.gapm[i + j * k_avx512_float_read_concurrency] = rovaca::Context<float>::_(1.0) - rovaca::Context<float>::ph2pr[_c];
            tc.mi[i + j * k_avx512_float_read_concurrency] = rovaca::Context<float>::ph2pr[_i];
            tc.ii[i + j * k_avx512_float_read_concurrency] = rovaca::Context<float>::ph2pr[_c];
            tc.md[i + j * k_avx512_float_read_concurrency] = rovaca::Context<float>::ph2pr[_d];
            tc.dd[i + j * k_avx512_float_read_concurrency] = rovaca::Context<float>::ph2pr[_c];

            float dist = rovaca::Context<float>::ph2pr[_q];
            tc._1_distm[i + j * k_avx512_float_read_concurrency] = 1.0f - dist;
            tc.distm[i + j * k_avx512_float_read_concurrency] = dist / 3.0f;
        }
    }
}

void debug_tools_print_arr(const uint8_t* arr, int32_t len, uint8_t sub_num = 0)
{
    for (int32_t i{0}; i < len; ++i) {
        printf("%c", arr[i] + sub_num);
    }
    printf("\n");
}

void debug_tools_print_temp(const ::TestCase& temp)
{
    printf("%d \t %d\n", temp.rslen, temp.haplen);
    debug_tools_print_arr(temp.hap, temp.haplen);
    debug_tools_print_arr(temp.rs, temp.rslen);
    debug_tools_print_arr(temp.q, temp.rslen, 33);
    debug_tools_print_arr(temp.i, temp.rslen, 33);
    debug_tools_print_arr(temp.d, temp.rslen, 33);
    debug_tools_print_arr(temp.c, temp.rslen, 33);
    fflush(stdout);
}

DoubleVector2D call_old_pairhmm(const HaplotypeVector& hs, ReadVector& rs, int min_quality_threshold, PcrIndelModel pcr_option,
                                pMemoryPool pool)
{
    if (hs.empty() || rs.empty()) {
        return {};
    }

    init_native();

    size_t reads_num = rs.size();
    size_t haplotypes_num = hs.size();
    DoubleVector2D final_likelihood{haplotypes_num, DoubleVector{reads_num, NAN, pool}, pool};
    auto lp = static_cast<rovaca::RovacaMemoryPool*>(pool);
    rovaca::MemoryPoolGuard pairhmm_guard{lp};

    int32_t max_len = INT32_MIN;
    std::for_each(rs.begin(), rs.end(), [&](pReadRecord r) { max_len = std::max(max_len, r->seq_length()); });

    pBases read_qual = new ALLOC_FLEXIBLE_IN_POOL(pool, Bases, max_len, uint8_t) Bases{(uint32_t)max_len};
    pBases ins_gops = new ALLOC_FLEXIBLE_IN_POOL(pool, Bases, max_len, uint8_t) Bases{(uint32_t)max_len};
    pBases del_gops;
    pBases gap_conts = new ALLOC_FLEXIBLE_IN_POOL(pool, Bases, max_len, uint8_t) Bases{(uint32_t)max_len};
    pBases read_bases = new ALLOC_FLEXIBLE_IN_POOL(pool, Bases, max_len, uint8_t) Bases{(uint32_t)max_len};

    DoubleVector2D likelihood_array(reads_num, DoubleVector{haplotypes_num, NAN, pool}, pool);

    uint8_t mq;
    ::TestCase temp{};
    pReadRecord read;
    pBases hap_bases;
    double result_final, result_double;
    float result_float;

    for (size_t i = 0; i < reads_num; ++i) {
        read = rs.at(i);
        read->decode_base(read_bases);
        read_qual->num = read_bases->num;
        read_qual->data[read_qual->num] = '\0';
        memcpy(read_qual->data, read->qual(), read_qual->num * sizeof(uint8_t));

        // modify_read_qualities
        mq = read->mapping_quality();
        for (uint32_t k = 0; k < read_qual->num; ++k) {
            read_qual->data[k] = std::min(read_qual->data[k], mq);
            read_qual->data[k] = read_qual->data[k] < min_quality_threshold ? MIN_QUALITY : read_qual->data[k];
        }

        ins_gops->num = read_bases->num;
        ins_gops->data[ins_gops->num] = '\0';
        memset(ins_gops->data, '-', ins_gops->num);
        if (pcr_option != PcrIndelModel::NONE) {
            apply_pcr_error_model(read_bases, ins_gops, pcr_option);
        }

        del_gops = ins_gops;

        gap_conts->num = read_bases->num;
        gap_conts->data[gap_conts->num] = '\0';
        memset(gap_conts->data, '+' - 33, gap_conts->num);

        std::pmr::vector<double>& likelihood_array_at_i = likelihood_array[i];
        for (size_t j = 0; j < haplotypes_num; ++j) {
            hap_bases = hs.at(j)->get_display_string();

            temp.rslen = read_bases->num;
            temp.haplen = hap_bases->num;
            temp.q = read_qual->data;
            temp.i = ins_gops->data;
            temp.d = del_gops->data;
            temp.c = gap_conts->data;
            temp.hap = hap_bases->data;
            temp.rs = read_bases->data;
            result_float = s_compute_full_prob_float_old_pairhmm(temp);
            if (result_float < MIN_ACCEPTED) {
                result_double = s_compute_full_prob_double_old_pairhmm(temp);
                result_final = log10(result_double) - LOG10_INITIAL_CONSTANT_D;
            }
            else {
                result_final = (double)(log10f(result_float) - LOG10_INITIAL_CONSTANT_F);
            }
            likelihood_array_at_i[j] = result_final;
        }
    }

    normalize_likelihoods(likelihood_array);
    filter_poorly_modelled_evidence(rs, likelihood_array, lp);
    transpose_likelihood_matrix(likelihood_array, final_likelihood);

    return final_likelihood;
}

DoubleVector2D call_new_pairhmm002(const HaplotypeVector& hs, ReadVector& rs, int32_t min_quality_threshold, PcrIndelModel pcr_option,
                                   pMemoryPool pool)
{
    if (rs.size() <= 2 * k_avx512_float_read_concurrency) {
        return call_old_pairhmm(hs, rs, min_quality_threshold, pcr_option, pool);
    }
    assert(rs.back()->seq_length() >= rs.front()->seq_length());

    init_native();

    size_t hap_num = hs.size();
    size_t read_num = rs.size();
    DoubleVector2D result(hap_num, DoubleVector(read_num, NAN), pool);
    auto lp = static_cast<rovaca::RovacaMemoryPool*>(pool);
    rovaca::MemoryPoolGuard pairhmm_guard{lp};

    uint32_t max_hap_len = 0;
    uint32_t max_read_len = static_cast<uint32_t>(rs.back()->seq_length());
    for (auto h : hs) {
        max_hap_len = std::max(max_hap_len, h->length());
    }

    DoubleVector2D temp_result(read_num, DoubleVector(hap_num, NAN), pool);

    ::TestCase temp{};
    uint32_t hap_len;
    rovaca::pHaplotype h;
    float result_float;
    rovaca::TestCase tc{};
    rovaca::pBases hap_base;
    rovaca::pReadRecord read;
    rovaca::Context<float> ctx{};
    double result_double, result_final;
    float result_per_loop[k_avx512_float_read_concurrency]{};

    init_test_case(tc, max_read_len, max_hap_len, pool);

    rovaca::pBases read_base[k_avx512_float_read_concurrency]{};
    rovaca::pBases read_qual[k_avx512_float_read_concurrency]{};
    rovaca::pBases del_gops[k_avx512_float_read_concurrency]{};
    rovaca::pBases ins_gops[k_avx512_float_read_concurrency]{};
    rovaca::pBases gap_conts[k_avx512_float_read_concurrency]{};
    for (uint32_t k{0}; k < k_avx512_float_read_concurrency; ++k) {
        read_base[k] = new ALLOC_FLEXIBLE_IN_POOL(pool, Bases, max_read_len, uint8_t) Bases{max_read_len};
        read_qual[k] = new ALLOC_FLEXIBLE_IN_POOL(pool, Bases, max_read_len, uint8_t) Bases{max_read_len};
        ins_gops[k] = new ALLOC_FLEXIBLE_IN_POOL(pool, Bases, max_read_len, uint8_t) Bases{max_read_len};
        gap_conts[k] = new ALLOC_FLEXIBLE_IN_POOL(pool, Bases, max_read_len, uint8_t) Bases{max_read_len};
    }

    std::pmr::vector<ReadLenInfo> histogram = generate_histogram(rs, read_num, max_read_len, pool);

    std::pmr::vector<ReadGroup> group = generate_new_pairhmm_idx(histogram, pool);

    for (const ReadGroup& rg : group) {
    // Decompress read
        decompression_data(rs, tc, read_base, read_qual, ins_gops, del_gops, gap_conts, rg.idxs, pcr_option, min_quality_threshold);

    // Arrange read
        prepare_data(tc, ctx, read_base, read_qual, ins_gops, del_gops, gap_conts);

    // Loop through haplotypes and calculate likelihood
        for (uint32_t hap_idx{0}; hap_idx < hap_num; ++hap_idx) {
            h = hs.at(hap_idx);
            hap_base = h->get_display_string();

            hap_len = h->length();
            tc.min_hap_len = tc.max_hap_len = hap_len;

            // Initialize haps, here hap only takes itself
            for (uint32_t hap_base_idx{0}; hap_base_idx < hap_len; ++hap_base_idx) {
                tc.haps[hap_base_idx] = convert_arr[hap_base->data[hap_base_idx]];
            }

            rovaca::compute_full_prob_avx512_float_1xn(tc, result_per_loop);

            for (uint32_t kk{0}; kk < k_avx512_float_read_concurrency; ++kk) {
                result_float = result_per_loop[kk];

                if (result_float < MIN_ACCEPTED) {
                    temp.rslen = read_base[kk]->num;
                    temp.haplen = hap_base->num;
                    temp.q = read_qual[kk]->data;
                    temp.i = ins_gops[kk]->data;
                    temp.d = del_gops[kk]->data;
                    temp.c = gap_conts[kk]->data;
                    temp.hap = hap_base->data;
                    temp.rs = read_base[kk]->data;

                    result_double = s_compute_full_prob_double_old_pairhmm(temp);
                    result_final = log10(result_double) - LOG10_INITIAL_CONSTANT_D;
                }
                else {
                    result_final = (double)(log10f(result_float) - LOG10_INITIAL_CONSTANT_F);
                }
                temp_result[rg.idxs[kk]][hap_idx] = result_final;
            }
        }
    }

    for (ReadLenInfo& info : histogram) {
        if (0 == info.count) {
            continue;
        }

        for (uint32_t offset{0}; offset < info.count; ++offset) {
            uint32_t read_idx = info.idx + offset;
            read = rs.at(read_idx);
            read->decode_base(read_base[0]);
            read->decode_qual(read_qual[0]);
            read->ins_gops(ins_gops[0]);
            del_gops[0] = ins_gops[0];
            // rs[read_idx]->del_gops(del_gops[k]);
            read->gap_conts(gap_conts[0]);

            // modify_read_qualities
            uint8_t mq = rs[read_idx]->mapping_quality();
            for (uint32_t x = 0; x < read_qual[0]->num; ++x) {
                read_qual[0]->data[x] = std::min(read_qual[0]->data[x], mq);
                read_qual[0]->data[x] = read_qual[0]->data[x] < min_quality_threshold ? MIN_QUALITY : read_qual[0]->data[x];
            }

            if (pcr_option != PcrIndelModel::NONE) {
                apply_pcr_error_model(read_base[0], ins_gops[0], pcr_option);
            }

            DoubleVector& likelihoods_at_read = temp_result.at(read_idx);

            for (uint32_t hap_idx{0}; hap_idx < hap_num; ++hap_idx) {
                hap_base = hs.at(hap_idx)->get_display_string();

                temp.rslen = read_base[0]->num;
                temp.haplen = hap_base->num;
                temp.q = read_qual[0]->data;
                temp.i = ins_gops[0]->data;
                temp.d = del_gops[0]->data;
                temp.c = gap_conts[0]->data;
                temp.hap = hap_base->data;
                temp.rs = read_base[0]->data;

                result_float = s_compute_full_prob_float_old_pairhmm(temp);
                if (result_float < MIN_ACCEPTED) {
                    result_double = s_compute_full_prob_double_old_pairhmm(temp);
                    result_final = log10(result_double) - LOG10_INITIAL_CONSTANT_D;
                }
                else {
                    result_final = (double)(log10f(result_float) - LOG10_INITIAL_CONSTANT_F);
                }
                likelihoods_at_read[hap_idx] = result_final;
            }
        }
    }

    normalize_likelihoods(temp_result);
    filter_poorly_modelled_evidence(rs, temp_result, lp);
    transpose_likelihood_matrix(temp_result, result);

    return result;
}

}  // namespace rovaca