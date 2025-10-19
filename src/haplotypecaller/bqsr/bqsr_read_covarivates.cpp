#include "bqsr_read_covarivates.h"

static int simple_base_to_base_index(uint8_t base) { return g_bqsr_seq_table[base]; }

void ReadCovariates::record_rg_covariate_value()
{
    int key = 0;
    current_covariate_index = 0;
    for (uint32_t i = 0; i < read_len; ++i) {
        add_covariate(key, key, key, i);
    }
}
void ReadCovariates::record_qual_score_covariate_value()
{
    current_covariate_index = 1;
    for (uint32_t i = 0; i < read_len; ++i) {
        add_covariate(read_base_quality[i], 0, 0, i);
    }
}
void ReadCovariates::record_context_covariate_value()
{
    current_covariate_index = 2;
    for (uint32_t i = 0; i < read_len_after_clip; ++i) {
        int read_offset = negative_strand ? (read_len_after_clip - i - 1) : i;
        add_covariate(i >= mismatch_keys.size() ? 0 : mismatch_keys[i], 0, 0, read_offset);
    }
}
void ReadCovariates::record_cycle_covariate_value()
{
    current_covariate_index = 3;
    for (uint32_t i = 0; i < read_len; ++i) {
        add_covariate(cycle_key(i, read, false, BQSR_MAXIMUM_CYCLE_VALUE), 0, 0, i);
    }
}

void ReadCovariates::add_covariate(int mismatch, __attribute((unused)) int insertion, __attribute((unused)) int deletion, int read_offset)
{
    auto& match_keys = keys[read_offset];
    match_keys[current_covariate_index] = mismatch;

    //     #ifdef __COUNT_INDEL__
    //     auto& insert_keys = keys[BQSR_INSERT][read_offset];
    //     auto& deletion_keys = keys[BQSR_DELETION][read_offset];
    //     insert_keys[current_covariate_index] = insertion;
    //     deletion_keys[current_covariate_index] = deletion;
    // #else
    //     (void)insertion;
    //     (void)deletion;
    // #endif
}

int get_stranded_clipped_bases(bam1_t* read, uint8_t* clipped_base, int low_qual_tail)
{
    int read_len = read->core.l_qseq;
    int left_clip_index = 0;
    int right_clip_index = read_len - 1;
    uint8_t* qual = bam_get_qual(read);

    while (right_clip_index >= 0 && qual[right_clip_index] <= low_qual_tail) {
        right_clip_index--;
    }
    while (left_clip_index < read_len && qual[left_clip_index] <= low_qual_tail) {
        left_clip_index++;
    }

    if (left_clip_index > right_clip_index) return -1;
    if (right_clip_index < read_len - 1) {
        clip_low_qual_ends_in_method_WRITE_NS(read, clipped_base, right_clip_index + 1, read_len - 1);
    }
    if (left_clip_index > 0) {
        clip_low_qual_ends_in_method_WRITE_NS(read, clipped_base, 0, left_clip_index - 1);
    }
    if (bam_is_rev(read)) reverse_bam_seq(read, clipped_base);
    return read_len;
}

void context_with(std::vector<int>& mismatch_keys, int read_len, uint8_t* bases, int context_size, int mask)
{
    for (int i = 1; i < context_size && i <= read_len; ++i) {
        mismatch_keys[i - 1] = -1;
    }
    if (read_len < context_size) return;
    int new_base_offset = 2 * (context_size - 1) + LENGTH_BITS;
    int current_key = bqsr_covariate_key_from_context(bases, 0, context_size);
    mismatch_keys[context_size - 1] = current_key;

    int current_penalty = 0;
    if (current_key == -1) {
        current_key = 0;
        current_penalty = context_size - 1;
        int offset = new_base_offset;
        int base_index;
        while ((base_index = simple_base_to_base_index(bases[current_penalty])) != -1) {
            current_key |= (base_index << offset);
            offset -= 2;
            current_penalty--;
        }
    }
    for (int current_index = context_size; current_index < read_len; current_index++) {
        int base_index = simple_base_to_base_index(bases[current_index]);
        if (base_index == -1) {
            current_penalty = context_size;
            current_key = 0;
        }
        else {
            current_key = (current_key >> 2) & mask;
            current_key |= (base_index << new_base_offset);
            current_key |= context_size;
        }

        if (current_penalty == 0) {
            mismatch_keys[current_index] = current_key;
        }
        else {
            current_penalty--;
            mismatch_keys[current_index] = -1;
        }
    }
    return;
}

int cycle_key(int base_number, bam1_t* read, bool indel_in_count, int max_cycle)
{
    bool is_negstrand = bam_is_rev(read);
    bool is_second_in_pair = ((read->core.flag & BAM_FPAIRED) != 0) && ((read->core.flag & BAM_FREAD2) != 0);
    int read_order_factor = is_second_in_pair ? -1 : 1;
    int increment;
    int cycle;
    if (is_negstrand) {
        cycle = read->core.l_qseq * read_order_factor;
        increment = -1 * read_order_factor;
    }
    else {
        cycle = read_order_factor;
        increment = read_order_factor;
    }
    cycle += base_number * increment;
    if (!indel_in_count) {
        return key_from_cycle(cycle, max_cycle);
    }
    // TODO: if indel are in count
    return -1;
}

// Encodes the cycle number as a key.
int key_from_cycle(int cycle, int max_cycle)
{
    int ret = std::abs(cycle);
    if (ret > max_cycle) {
        throw std::out_of_range("abs of cycle cannot beyond max_cycle.");
    }
    ret <<= 1;  // shift so we can add the "sign" bit
    if (cycle < 0) {
        ret++;  // negative cycles get the lower-most bit set
    }
    return ret;
}
