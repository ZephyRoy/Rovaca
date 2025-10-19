#ifndef DOWNSAMPLER_HC_H
#define DOWNSAMPLER_HC_H

#include <list>
#include <vector>

#include "../rovaca_logger/rovaca_logger.h"
#include "downsampler.h"

#define HC_RANDOM_MULTIPLIER       (0x5DEECE66DL)
#define HC_RANDOM_ADDEND           (0xBL)
#define HC_RANDOM_MASK             ((1L << 48) - 1)
#define HC_DOWN_SAMPLE_RANDOM_SEED (25170011922L)

class HCDownsampler : public Downsampler
{
private:
    uint32_t target_coverage;
    std::list<bam1_t*> finalized_items;
    std::list<bam1_t*> reservoir;
    bam1_t* previous_reads;
    bool eof_input_stream;
    uint32_t total_reads_seen;
    uint64_t hc_random_seed;

public:
    HCDownsampler(uint32_t target_sample_size);
    ~HCDownsampler() override;

    /*!
     * @brief Add a read item to downsample reservoir.
     * @param bam1_t* item
     * @return read to be abandoned.
     */
    bam1_t* submit(bam1_t* item) override;
    bool has_finalized_items() override;
    void input_end_signal() override;
    int consume_finalized_items(std::list<bam1_t*>& cache) override;

private:
    void handle_positional_change(bam1_t* read);
    void finalize_reservoir(bool expect);
    uint32_t apply_random_method(uint32_t bound);
    uint32_t apply_random_method_next(uint32_t bits);
};

static bool read_has_no_assigned_position(bam1_t* read) { return read->core.tid < 0 || read->core.pos < 0; }

static uint32_t compair_reads_coordinates(bam1_t* one, bam1_t* other)
{
    int tid_one = one->core.tid;
    int tid_other = other->core.tid;

    if (tid_one == -1) {
        return (tid_other == -1 ? 0 : 1);
    }
    else if (tid_other == -1) {
        return -1;
    }

    int index_diff = tid_one - tid_other;
    if (index_diff != 0) {
        return index_diff;
    }

    if (one->core.pos > other->core.pos) {
        return 1;
    }
    else if (one->core.pos < other->core.pos) {
        return -1;
    }
    else {
        return 0;
    }
    return 0;
}

uint32_t HCDownsampler::apply_random_method(uint32_t bound)
{
    int32_t r;
    int32_t m;

    if (bound <= 0) {
        return -1;
    }

    r = apply_random_method_next(31);
    m = bound - 1;

    if ((bound & m) == 0) {  // i.e., bound is a power of 2
        r = (int32_t)((bound * (int64_t)r) >> 31);
    }
    else {
        int32_t u = r;
        for (; u - (r = u % bound) + m < 0; u = apply_random_method_next(31)) {}
    }
    return r;
}

uint32_t HCDownsampler::apply_random_method_next(uint32_t bits)
{
    uint64_t nextseed;
    uint64_t seed = hc_random_seed;
    uint32_t ret_u;
    int32_t ret;

    nextseed = (seed * HC_RANDOM_MULTIPLIER + HC_RANDOM_ADDEND) & HC_RANDOM_MASK;
    seed = nextseed;

    ret_u = (int32_t)(nextseed >> (48 - bits));
    ret = ret_u;
    hc_random_seed = nextseed;

    return ret;
}

HCDownsampler::HCDownsampler(uint32_t target_sample_size)
    : target_coverage(target_sample_size)
    , finalized_items({})
    , reservoir({})
    , previous_reads(nullptr)
    , eof_input_stream(false)
    , total_reads_seen(0)
    , hc_random_seed(HC_DOWN_SAMPLE_RANDOM_SEED)
{
    // finalized_items.reserve(50);
}

HCDownsampler::~HCDownsampler() {}

bam1_t* HCDownsampler::submit(bam1_t* item)
{
    bam1_t* ret;
    if (item == nullptr) {
        RovacaLogger::warn("an empty read");
        return nullptr;
    }
    if (eof_input_stream) {
        RovacaLogger::warn("attempt to call submit after EOI");
        return nullptr;
    }
    handle_positional_change(item);
    // Ensure a normal read here.
    if (!read_has_no_assigned_position(item)) {
        total_reads_seen++;
        if (total_reads_seen <= target_coverage) {
            reservoir.emplace_back(item);
            previous_reads = item;
        }
        else {
            increment_number_of_discarded_items();
            uint32_t random_slot = apply_random_method(total_reads_seen);
            if (random_slot < target_coverage) {
                auto it = std::next(reservoir.begin(), random_slot);
                ret = *it;
                *it = item;
                previous_reads = item;
            }
            else {
                ret = item;
            }
            return ret;
        }
    }
    else {
        finalized_items.emplace_back(item);
    }
    return nullptr;
}

bool HCDownsampler::has_finalized_items() { return !finalized_items.empty(); }

int HCDownsampler::consume_finalized_items(std::list<bam1_t*>& cache)
{
    cache = std::move(finalized_items);
    return (int)cache.size();
}

void HCDownsampler::input_end_signal() { finalize_reservoir(false); }

void HCDownsampler::handle_positional_change(bam1_t* read)
{
    if (previous_reads != nullptr) {
        int read_cmp = compair_reads_coordinates(previous_reads, read);
        if (read_cmp == 1) {
            RovacaLogger::warn("reads must be coordinate sorted");
        }
        else if (read_cmp != 0) {
            finalize_reservoir(true);
        }
    }
}

void HCDownsampler::finalize_reservoir(bool expect)
{
    eof_input_stream = true;

    if (expect && reservoir.empty()) {
        RovacaLogger::warn("expected downsampled items to be present when none are");
    }

    if (!reservoir.empty()) {
        finalized_items.splice(finalized_items.end(), std::move(reservoir));
        total_reads_seen = 0;
    }

    reservoir.clear();
    eof_input_stream = !reservoir.empty();

    previous_reads = nullptr;
}

#endif  // DOWNSAMPLER_HC_H