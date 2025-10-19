#include "alignment_utils.h"

#include <algorithm>
#include <memory_resource>

#include "cigar_builder.h"
#include "genotype_macors.h"
#include "haplotype.h"
#include "rovaca_logger.h"
#include "read_clipper.h"
#include "read_record.h"
#include "simple_interval.h"
#include "smithwaterman_common.h"
#include "utils/base_utils.h"
#include "utils/cigar_utils.h"

namespace rovaca
{

#define PADDED_SIZE (1000)

namespace AlignmentUtils
{

class CigarPairTransform
{
    /* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
private:
    std::pmr::set<uint32_t> _op12, _op23;
    uint32_t _op13, _advance12, _advance23;

    /* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
public:
    static pCigarPairTransform create(uint32_t op12, uint32_t op23, uint32_t op13, uint32_t advance12, uint32_t advance23, pMemoryPool pool)
    {
        return new ALLOC_TYPE_IN_POOL(pool, CigarPairTransform) CigarPairTransform{op12, op23, op13, advance12, advance23, pool};
    }

    uint32_t op13() const { return _op13; }
    uint32_t advance12() const { return _advance12; }
    uint32_t advance23() const { return _advance23; }
    const std::pmr::set<uint32_t>& op12() const { return _op12; }
    const std::pmr::set<uint32_t>& op23() const { return _op23; }

    /* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
private:
    CigarPairTransform(uint32_t op12, uint32_t op23, uint32_t op13, uint32_t advance12, uint32_t advance23, pMemoryPool pool)
        : _op12(pool)
        , _op23(pool)
        , _op13(op13)
        , _advance12(advance12)
        , _advance23(advance23)
    {
        _op12 = get_cigar_set(op12, pool);
        _op23 = get_cigar_set(op23, pool);
    }

    static std::pmr::set<uint32_t> get_cigar_set(uint32_t master_op, pMemoryPool pool)
    {
        switch (master_op) {
            case BAM_CMATCH: {
                return {{BAM_CMATCH, BAM_CEQUAL, BAM_CDIFF}, pool};
            }
            case BAM_CINS: {
                return {{BAM_CINS, BAM_CSOFT_CLIP}, pool};
            }
            case BAM_CDEL: {
                return {{BAM_CDEL}, pool};
            }
            default: {
                RovacaLogger::error("unexpected state: {}", master_op);
                exit(EXIT_FAILURE);
            }
        }
    }
};

class CigarPairTransformManager
{
    static constexpr size_t s_buffer_size = 4096;

    /* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
private:
    uint8_t _data[s_buffer_size];
    std::pmr::monotonic_buffer_resource _pool;

    std::pmr::vector<pCigarPairTransform> _transformers;

    /* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
public:
    DISALLOW_COPY_AND_ASSIGN(CigarPairTransformManager);

    static pCigarPairTransformManager get_instance()
    {
        static std::unique_ptr<CigarPairTransformManager> m(new CigarPairTransformManager{});
        return m.get();
    }

    const std::pmr::vector<pCigarPairTransform>& cigar_pair_transformers() const { return _transformers; }

    /* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
private:
    CigarPairTransformManager()
        : _data()
        , _pool(_data, s_buffer_size, std::pmr::null_memory_resource())
        , _transformers(&_pool)
    {
        _transformers.reserve(9);

        //
        // op12 is a match
        //
        // 3: xxx B yyy
        // ^^^^^^^^^^^^
        // 2: xxx M yyy
        // 1: xxx M yyy
        _transformers.emplace_back(CigarPairTransform::create(BAM_CMATCH, BAM_CMATCH, BAM_CMATCH, 1, 1, &_pool));
        // 3: xxx I yyy
        // ^^^^^^^^^^^^
        // 2: xxx I yyy
        // 1: xxx M yyy
        _transformers.emplace_back(CigarPairTransform::create(BAM_CMATCH, BAM_CINS, BAM_CINS, 1, 1, &_pool));
        // 3: xxx D yyy
        // ^^^^^^^^^^^^
        // 2: xxx D yyy
        // 1: xxx M yyy
        _transformers.emplace_back(CigarPairTransform::create(BAM_CMATCH, BAM_CDEL, BAM_CDEL, 0, 1, &_pool));

        //
        // op12 is a deletion
        //
        // 3: xxx D M yyy
        // ^^^^^^^^^^^^
        // 2: xxx M yyy
        // 1: xxx D yyy
        _transformers.emplace_back(CigarPairTransform::create(BAM_CDEL, BAM_CMATCH, BAM_CDEL, 1, 1, &_pool));
        // 3: xxx D2 D1 yyy
        // ^^^^^^^^^^^^
        // 2: xxx D2 yyy
        // 1: xxx D1 yyy
        _transformers.emplace_back(CigarPairTransform::create(BAM_CDEL, BAM_CDEL, BAM_CDEL, 0, 1, &_pool));
        // 3: xxx X yyy => no-op, we skip emitting anything here
        // ^^^^^^^^^^^^
        // 2: xxx I yyy
        // 1: xxx D yyy
        _transformers.emplace_back(CigarPairTransform::create(BAM_CDEL, BAM_CINS, BAM_CUNINITIALIZE, 1, 1, &_pool));

        //
        // op12 is a insertion
        //
        // 3: xxx I M yyy
        // ^^^^^^^^^^^^
        // 2: xxx M yyy
        // 1: xxx I yyy
        _transformers.emplace_back(CigarPairTransform::create(BAM_CINS, BAM_CMATCH, BAM_CINS, 1, 0, &_pool));
        // 3: xxx I D yyy
        // ^^^^^^^^^^^^
        // 2: xxx D yyy
        // 1: xxx I yyy
        _transformers.emplace_back(CigarPairTransform::create(BAM_CINS, BAM_CDEL, BAM_CINS, 1, 0, &_pool));
        // 3: xxx I1 I2 yyy
        // ^^^^^^^^^^^^
        // 2: xxx I2 yyy
        // 1: xxx I1 yyy
        _transformers.emplace_back(CigarPairTransform::create(BAM_CINS, BAM_CINS, BAM_CINS, 1, 0, &_pool));
    }
};

}  // namespace AlignmentUtils

uint32_t AlignmentUtils::length_on_read(uint32_t cigar_elenment)
{
    return consumes_read_bases(bam_cigar_op(cigar_elenment)) ? bam_cigar_oplen(cigar_elenment) : 0;
}

uint32_t AlignmentUtils::length_on_reference(uint32_t cigar_elenment)
{
    return consumes_ref_bases(bam_cigar_op(cigar_elenment)) ? bam_cigar_oplen(cigar_elenment) : 0;
}

/*!
 * @brief given a read's first aligned base on an alt haplotype, find the first aligned base in the reference haplotype.  this method
 * assumes that the alt haplotype and reference haplotype start at the same place.  that is, the alt haplotype starts at index 0 within the
 * reference base array.
 * @param haplotype_vs_ref_cigar
 * @param read_start_on_haplotype
 * @return
 */
int64_t AlignmentUtils::read_start_on_reference_haplotype(pCigar haplotype_vs_ref_cigar, int64_t read_start_on_haplotype)
{
    if (0 == read_start_on_haplotype) {
        return 0;
    }

    // move forward in the haplotype vs ref cigar until we have consumed enough haplotype bases to reach the read start the number of
    // reference bases consumed during this traversal gives us the reference start
    int64_t haplotype_bases_consumed = 0;
    int64_t ref_bases_consumed_before_start = 0;

    for (uint32_t i = 0; i < haplotype_vs_ref_cigar->num; i++) {
        uint32_t cigar_element = haplotype_vs_ref_cigar->data[i];

        ref_bases_consumed_before_start += (int64_t)length_on_reference(cigar_element);
        haplotype_bases_consumed += (int64_t)length_on_read(cigar_element);

        if (haplotype_bases_consumed >= read_start_on_haplotype) {
            int64_t excess = consumes_ref_bases(bam_cigar_op(cigar_element)) ? haplotype_bases_consumed - read_start_on_haplotype : 0;
            return ref_bases_consumed_before_start - excess;
        }
    }

    CHECK_CONDITION_EXIT(true, "Cigar doesn't reach the read start");
}

/**
 * Trim cigar down to one that starts at start base in the cigar and extends to (inclusive) end base
 *
 * @param cigar a non-null Cigar to trim down
 * @param start Where should we start keeping bases in the cigar (inclusive)?  The first position is 0
 * @param end Where should we stop keeping bases in the cigar (inclusive)?  The maximum value is cigar.getLength() - 1
 * @return a new Cigar containing == start - end + 1 reads
 */
pCigar AlignmentUtils::trim_cigar_by_bases(pCigar cigar, int64_t start, int64_t end, pMemoryPool pool)
{
    return trim_cigar(cigar, start, end, false, pool);
}

pCigar AlignmentUtils::trim_cigar_by_reference(pCigar cigar, int64_t start, int64_t end, pMemoryPool pool)
{
    return trim_cigar(cigar, start, end, true, pool);
}

/**
 * Workhorse for trimCigarByBases and trimCigarByReference
 *
 * @param cigar a non-null Cigar to trim down
 * @param start Where should we start keeping bases in the cigar (inclusive)?  The first position is 0
 * @param end Where should we stop keeping bases in the cigar (inclusive)?  The maximum value is cigar.getLength() - 1
 * @param byReference should start and end be interpreted as position in the reference or the read to trim to/from?
 * @return a non-null cigar
 */
pCigar AlignmentUtils::trim_cigar(pCigar cigar, int64_t start, int64_t end, bool by_reference, pMemoryPool pool)
{
    CHECK_CONDITION_EXIT(start < 0, "start position can't be negative");
    CHECK_CONDITION_EXIT(end < start, "end is smaller than start");

    pCigarBuilder builder = CigarBuilder::create(pool);

    // these variables track the inclusive start and exclusive end of the current cigar element in reference (if byReference) or read
    // (otherwise) coordinates
    int64_t element_start;    // inclusive
    int64_t element_end = 0;  // exclusive -- start of next element
    int64_t overlap_length;
    uint32_t element;
    for (uint32_t i = 0; i < cigar->num; i++) {
        element = cigar->data[i];
        element_start = element_end;
        element_end = element_start + (by_reference ? length_on_reference(element) : length_on_read(element));

        // we are careful to include zero-length elements at both ends, that is, elements with elementStart == elementEnd == start and
        // elementStart == elementEnd == end + 1
        if (element_end < start || (element_end == start && element_start < start)) {
            continue;
        }
        else if (element_start > end && element_end > end + 1) {
            break;
        }

        overlap_length =
            element_end == element_start ? bam_cigar_oplen(element) : std::min(end + 1, element_end) - std::max(start, element_start);

        builder->add(overlap_length << BAM_CIGAR_SHIFT | bam_cigar_op(element));
    }

    CHECK_CONDITION_EXIT(element_end <= end, "cigar elements don't reach end position (inclusive) %d", end);

    return builder->make_and_record_deletions_removed_result().cigar;
}

pCigar AlignmentUtils::apply_cigar_to_cigar(pCigar first_to_second, pCigar second_to_third, pMemoryPool pool)
{
    pCigarBuilder new_elements = CigarBuilder::create(pool);
    uint32_t num_elements_12 = first_to_second->num;
    uint32_t num_elements_23 = second_to_third->num;

    uint32_t cigar12_I = 0, cigar23_I = 0, elt12_I = 0, elt23_I = 0;
    uint32_t elt12, elt23;

    while (cigar12_I < num_elements_12 && cigar23_I < num_elements_23) {
        elt12 = first_to_second->data[cigar12_I];  // elt12 and elt23 are cigar elements
        elt23 = second_to_third->data[cigar23_I];

        pCigarPairTransform transformer = get_transformer(bam_cigar_op(elt12), bam_cigar_op(elt23));

        if (transformer->op13() != BAM_CUNINITIALIZE) {  // skip no ops
            new_elements->add((1 << BAM_CIGAR_SHIFT | transformer->op13()));
        }

        elt12_I += transformer->advance12();
        elt23_I += transformer->advance23();

        // if have exhausted our current element, advance to the next one
        if (elt12_I == bam_cigar_oplen(elt12)) {
            cigar12_I++;
            elt12_I = 0;
        }

        if (elt23_I == bam_cigar_oplen(elt23)) {
            cigar23_I++;
            elt23_I = 0;
        }
    }

    return new_elements->make();
}

AlignmentUtils::pCigarPairTransform AlignmentUtils::get_transformer(uint32_t op12, uint32_t op23)
{
    const auto& cigar_pair_transformers = CigarPairTransformManager::get_instance()->cigar_pair_transformers();
    for (pCigarPairTransform transform : cigar_pair_transformers) {
        if (transform->op12().count(op12) && transform->op23().count(op23)) {
            return transform;
        }
    }
    CHECK_CONDITION_EXIT(true, "No transformer for operators %d and %d", op12, op23);
}

bool AlignmentUtils::last_base_on_right_is_same(const std::pmr::vector<uint8_t*>& sequences, std::pmr::vector<IndexRange>& bounds)
{
    uint8_t* seq;
    int32_t idx = bounds.at(0).get_end() - 1;
    uint8_t last_base_on_right = sequences.at(0)[idx];
    std::size_t n, sequences_size = sequences.size();
    for (n = 1; n < sequences_size; ++n) {
        seq = sequences.at(n);
        idx = bounds.at(n).get_end() - 1;
        if (seq[idx] != last_base_on_right) {
            return false;
        }
    }
    return true;
}

bool AlignmentUtils::first_base_on_left_is_same(const std::pmr::vector<uint8_t*>& sequences, std::pmr::vector<IndexRange>& bounds)
{
    uint8_t* seq;
    int32_t idx = bounds.at(0).get_start();
    uint8_t first_base_on_left = sequences.at(0)[idx];
    std::size_t n, sequences_size = sequences.size();
    for (n = 1; n < sequences_size; ++n) {
        seq = sequences.at(n);
        idx = bounds.at(n).get_start();
        if (seq[idx] != first_base_on_left) {
            return false;
        }
    }
    return true;
}

bool AlignmentUtils::next_base_on_left_is_same(const std::pmr::vector<uint8_t*>& sequences, std::pmr::vector<IndexRange>& bounds)
{
    uint8_t* seq;
    int32_t idx = bounds.at(0).get_start() - 1;
    uint8_t next_base_on_left = sequences.at(0)[idx];
    std::size_t n, sequences_size = sequences.size();
    for (n = 1; n < sequences_size; ++n) {
        seq = sequences.at(n);
        idx = bounds.at(n).get_start() - 1;
        if (seq[idx] != next_base_on_left) {
            return false;
        }
    }
    return true;
}

std::pair<int32_t, int32_t> AlignmentUtils::normalize_alleles(const std::pmr::vector<uint8_t*>& sequences,
                                                              std::pmr::vector<IndexRange>& bounds, int32_t max_shift, bool trim)
{
    CHECK_CONDITION_EXIT(sequences.empty(), "sequences is empty");
    CHECK_CONDITION_EXIT(sequences.size() != bounds.size(), "must have one initial allele range per sequence");
    for (const IndexRange& bound : bounds) {
        CHECK_CONDITION_EXIT(bound.get_start() < (int32_t)max_shift, "max_shift goes past the start of a sequence");
    }

    int32_t start_shift = 0, end_shift = 0;

    // consume any redundant shared bases at the end of the alleles
    int32_t min_size = INT32_MAX;
    std::for_each(bounds.begin(), bounds.end(), [&min_size](const IndexRange& bound) { min_size = std::min(min_size, bound.size()); });

    while (trim && min_size > 0 && last_base_on_right_is_same(sequences, bounds)) {
        std::for_each(bounds.begin(), bounds.end(), [](IndexRange& bound) { bound.shift_end_left(1); });
        --min_size;
        ++end_shift;
    }

    while (trim && min_size > 0 && first_base_on_left_is_same(sequences, bounds)) {
        std::for_each(bounds.begin(), bounds.end(), [](IndexRange& bound) { bound.shift_start(1); });
        --min_size;
        --start_shift;
    }

    // we shift left as long as the last bases on the right are equal among all sequences and the next bases on the left are all
    // equal. if a sequence is empty (eg the reference relative to an insertion alt allele) the last base on the right is the next
    // base on the left
    while (start_shift < max_shift && next_base_on_left_is_same(sequences, bounds) && last_base_on_right_is_same(sequences, bounds)) {
        std::for_each(bounds.begin(), bounds.end(), [](IndexRange& bound) { bound.shift_left(1); });
        ++start_shift;
        ++end_shift;
    }

    return {start_shift, end_shift};
}

CigarBuilder::Result AlignmentUtils::left_align_indels(pCigar cigar, pBases ref_bases, pBases read_bases, int64_t read_start,
                                                       pMemoryPool pool)
{
    CHECK_CONDITION_EXIT(read_start < 0, "read start within reference base array must be non-negative");

    bool ret = false;
    // Find if there is INDEL existing in CIGAR; If no return original CIGAR
    for (uint32_t i = 0; i < cigar->num; i++) {
        if (cigar_op_is_indel(bam_cigar_op(cigar->data[i]))) {
            ret = true;
            break;
        }
    }
    if (!ret) {
        return {cigar, 0, 0};
    }

    // we need reference bases from the start of the read to the rightmost indel
    int64_t necessary_ref_length = read_start;
    for (int32_t i = (int32_t)cigar->num - 1; i >= 0; --i) {
        if (cigar_op_is_indel(bam_cigar_op(cigar->data[i]))) {
            break;
        }
        necessary_ref_length += length_on_reference(cigar->data[i]);
    }
    CHECK_CONDITION_EXIT(necessary_ref_length > (int64_t)ref_bases->num, "read goes past end of reference");

    // at this point, we are one base past the end of the read.  Now we traverse the cigar from right to left
    int32_t cigar_num_i32 = (int32_t)cigar->num;
    int32_t ref_length = (int32_t)bam_cigar2rlen((int32_t)cigar->num, cigar->data);
    int32_t read_length = (int32_t)read_bases->num;
    Uint32Vector result_right_to_left(pool);
    result_right_to_left.reserve(cigar_num_i32 * 3);

    std::pmr::vector<IndexRange> ranges(pool);
    ranges.reserve(2);
    ranges.emplace_back(read_start + ref_length, read_start + ref_length);
    ranges.emplace_back(read_length, read_length);
    IndexRange& ref_indel_range = ranges.at(0);
    IndexRange& read_indel_range = ranges.at(1);
    std::pmr::vector<uint8_t*> sequences{{ref_bases->data, read_bases->data}, pool};
    uint32_t element, element_op, element_len;
    for (int32_t n = cigar_num_i32 - 1; n >= 0; --n) {
        element = cigar->data[n];
        element_op = bam_cigar_op(element);
        element_len = bam_cigar_oplen(element);
        // if it's an indel, just accumulate the read and ref bases consumed.  We won't shift the indel until we hit an alignment
        // block or the read start.
        if (cigar_op_is_indel(element_op)) {
            ref_indel_range.shift_start_left((int32_t)length_on_reference(element));
            read_indel_range.shift_start_left((int32_t)length_on_read(element));
        }
        else if (ref_indel_range.size() == 0 && read_indel_range.size() == 0) {
            // no indel, just add the cigar element to the result
            result_right_to_left.emplace_back(element);
            ref_indel_range.shift_left((int32_t)length_on_reference(element));
            read_indel_range.shift_left((int32_t)length_on_read(element));
        }
        else {
            // there's an indel that we need to left-align
            // we can left-align into match cigar elements but not into clips
            int32_t max_shift = cigar_op_is_alignment(element_op) ? (int32_t)element_len : 0;

            std::pair<int32_t, int32_t> shifts = normalize_alleles(sequences, ranges, max_shift, true);

            // account for new match alignments on the right due to left-alignment
            result_right_to_left.emplace_back((shifts.second << BAM_CIGAR_SHIFT | BAM_CMATCH));

            // emit if we didn't go all the way to the start of an alignment block OR we have reached clips OR we have reached the
            // start of the cigar
            bool emit_indel = (n == 0) || (shifts.first < max_shift) || (!cigar_op_is_alignment(element_op));
            int32_t new_match_on_left_due_to_trimming = shifts.first < 0 ? -shifts.first : 0;
            int32_t remaining_bases_on_left = shifts.first < 0 ? (int32_t)element_len : ((int32_t)element_len - shifts.first);

            if (emit_indel) {
                // some of this alignment block remains after left-alignment -- emit the indel
                result_right_to_left.emplace_back(ref_indel_range.size() << BAM_CIGAR_SHIFT | BAM_CDEL);
                result_right_to_left.emplace_back(read_indel_range.size() << BAM_CIGAR_SHIFT | BAM_CINS);

                // ref indel range is now empty and points to start of left-aligned indel
                ref_indel_range.shift_end_left(ref_indel_range.size());
                // read indel range is now empty and points to start of left-aligned indel
                read_indel_range.shift_end_left(read_indel_range.size());

                ref_indel_range.shift_left(new_match_on_left_due_to_trimming +
                                           (consumes_ref_bases(element_op) ? remaining_bases_on_left : 0));
                read_indel_range.shift_left(new_match_on_left_due_to_trimming +
                                            (consumes_read_bases(element_op) ? remaining_bases_on_left : 0));
                // now read and ref indel ranges are empty and point to end of element preceding this block
            }
            result_right_to_left.emplace_back(new_match_on_left_due_to_trimming << BAM_CIGAR_SHIFT | BAM_CMATCH);
            result_right_to_left.emplace_back(remaining_bases_on_left << BAM_CIGAR_SHIFT | element_op);
        }
    }

    // account for any indels at the start of the cigar that weren't processed because they have no adjacent non-indel element to
    // the left
    result_right_to_left.emplace_back(ref_indel_range.size() << BAM_CIGAR_SHIFT | BAM_CDEL);
    result_right_to_left.emplace_back(read_indel_range.size() << BAM_CIGAR_SHIFT | BAM_CINS);

    if (read_indel_range.get_start() != 0) {
        printf("failed\n");
    }
    CHECK_CONDITION_EXIT(read_indel_range.get_start() != 0, "given cigar does not account for all bases of the read");

    pCigarBuilder builder = CigarBuilder::create(pool);
    std::reverse(result_right_to_left.begin(), result_right_to_left.end());
    return builder->add_all(result_right_to_left.data(), result_right_to_left.size())->make_and_record_deletions_removed_result();
}

pCigar AlignmentUtils::append_clipped_elements_from_cigar_to_cigar(pCigar cigar_to_have_clipped_elements_added,
                                                                   const uint32_t* original_clipped_cigar, uint32_t cigar_num,
                                                                   pMemoryPool pool)
{
    uint32_t first_index = 0, last_index = cigar_num - 1;
    uint32_t first_element = original_clipped_cigar[first_index];
    uint32_t last_element = original_clipped_cigar[last_index];

    std::pmr::vector<uint32_t> read_to_ref_cigar_elements_with_hard_clips(pool);
    read_to_ref_cigar_elements_with_hard_clips.reserve(3 * cigar_num);

    while (cigar_op_is_clipping(bam_cigar_op(first_element)) && first_index != last_index) {
        read_to_ref_cigar_elements_with_hard_clips.emplace_back(first_element);
        first_element = original_clipped_cigar[++first_index];
    }

    for (uint32_t i = 0; i < cigar_to_have_clipped_elements_added->num; i++) {
        read_to_ref_cigar_elements_with_hard_clips.emplace_back(cigar_to_have_clipped_elements_added->data[i]);
    }

    std::pmr::vector<uint32_t> end_cigar_elements_to_reverse(pool);
    end_cigar_elements_to_reverse.reserve(cigar_num);

    while (cigar_op_is_clipping(bam_cigar_op(last_element)) && first_index != last_index) {
        end_cigar_elements_to_reverse.emplace_back(last_element);
        last_element = original_clipped_cigar[--last_index];
    }

    std::reverse(end_cigar_elements_to_reverse.begin(), end_cigar_elements_to_reverse.end());
    std::copy(end_cigar_elements_to_reverse.begin(), end_cigar_elements_to_reverse.end(),
              std::back_inserter(read_to_ref_cigar_elements_with_hard_clips));

    uint32_t num = (uint32_t)read_to_ref_cigar_elements_with_hard_clips.size();
    pCigar result = new ALLOC_FLEXIBLE_IN_POOL(pool, Cigar, num, uint32_t) Cigar{num};
    memcpy(result->data, read_to_ref_cigar_elements_with_hard_clips.data(), sizeof(uint32_t) * num);
    return result;
}

void AlignmentUtils::create_read_aligned_to_ref(pReadRecord original_read, pHaplotype haplotype, pHaplotype ref_haplotype,
                                                int64_t reference_start, bool is_informative, p_lib_sw_avx sw, pMemoryPool pool,
                                                pBamDataPool bam_pool)
{
    un_used(is_informative);
    CHECK_CONDITION_EXIT(reference_start < 1, "reference start much be >= 1 but got {}", reference_start);

    ReadClipper rc(original_read, pool, bam_pool);
    pReadRecord read_minus_soft_clips = rc.hard_clip_soft_clipped_bases();
    int32_t soft_clipped_bases = original_read->seq_length() - read_minus_soft_clips->seq_length();

    sw->len1 = (int16_t)haplotype->length();
    memcpy(sw->seq1, haplotype->get_display_string()->data, sw->len1);
    sw->len2 = (int16_t)read_minus_soft_clips->seq_length();
    pBases read_seq = read_minus_soft_clips->decode_to_str(pool);
    memcpy(sw->seq2, (char*)read_seq->data, sw->len2);
    sw->avx_function(sw);

    if (sw->alignment_offset == -1) {
        // RovacaLogger::warn("Smith-Waterman failed, returning original read");
        // pRealignedResult realigned_result = new ALLOC_TYPE_IN_POOL(pool, RealignedResult) RealignedResult{};
        //
        // realigned_result->new_start = original_read->get_start();
        // realigned_result->new_stop = original_read->get_stop();
        //
        // pCigarBuilder cigar_rebuiler = CigarBuilder::create(pool);
        //
        // for (uint32_t i = 0; i < sw->cigar_count; i++) {
        //     cigar_rebuiler->add(sw->cigar[i]);
        // }
        //
        // pCigar rebuilt_cigar = cigar_rebuiler->make();
        // realigned_result->realigned_cigar = rebuilt_cigar;
        //
        // original_read->_realigned_result = realigned_result;
        return;
    }

    pCigarBuilder sw_cigar_builder = CigarBuilder::create(pool);

    for (uint32_t i = 0; i < sw->cigar_count; i++) {
        sw_cigar_builder->add(sw->cigar[i]);
    }

    pCigar sw_cigar = sw_cigar_builder->make();

    pCigarBuilder hap_cigar_builder = CigarBuilder::create(pool);

    // compute here the read starts w.r.t. the reference from the SW result and the hap -> ref cigar
    hap_cigar_builder->add_all(haplotype->cigar()->data, haplotype->cigar()->num);
    hap_cigar_builder->add((PADDED_SIZE << BAM_CIGAR_SHIFT | BAM_CMATCH));
    pCigar right_padded_cigar = hap_cigar_builder->make();

    // this computes the number of reference bases before the read starts, based on the haplotype vs ref cigar
    // This cigar corresponds exactly to the readToRefCigarRaw, below.  One might wonder whether readToRefCigarRaw and
    // readToRefCigarClean ever imply different starts, which could occur if if the former has a leading deletion.  However,
    // according to the logic of applyCigarToCigar, this can only happen if the read has a leading deletion wrt its best haplotype,
    // which our SW aligner won't do, or if the read starts on a haplotype base that is in a deletion wrt to reference, which is nonsensical
    // since a base that exists is not a deletion.  Thus, there is nothing to worry about, in contrast to below where we do check
    // whether left-alignment shifted the start position.
    int64_t read_start_on_ref_hap = read_start_on_reference_haplotype(right_padded_cigar, sw->alignment_offset);

    int64_t read_start_on_reference = reference_start + haplotype->alignment_start_hap_wrt_ref() + read_start_on_ref_hap;

    // compute the read -> ref alignment by mapping read -> hap -> ref from the
    // SW of read -> hap mapped through the given by hap -> ref

    // This is the sub-cigar of the haplotype-to-ref alignment, with cigar elements before the read start removed.
    // Elements after the read end are kept.
    int64_t cigar_consume_read = bam_cigar2qlen((int32_t)right_padded_cigar->num, right_padded_cigar->data);
    pCigar haplotype_to_ref_cigar = trim_cigar_by_bases(right_padded_cigar, (int64_t)sw->alignment_offset, cigar_consume_read - 1, pool);

    pCigar read_to_ref_cigar = apply_cigar_to_cigar(sw_cigar, haplotype_to_ref_cigar, pool);

    // CigarBuilder::Result left_align_indels(pCigar cigar, pBases ref_bases, pBases read_bases, int64_t read_start, pMemoryPool pool);
    pBases ref_bases = ref_haplotype->get_display_string();
    pBases read_bases = read_minus_soft_clips->decode_to_str(pool);
    auto left_aligned_read_to_ref_cigar_result = left_align_indels(read_to_ref_cigar, ref_bases, read_bases, read_start_on_ref_hap, pool);
    pCigar lartrc = left_aligned_read_to_ref_cigar_result.cigar;  // left_aligned_read_to_ref_cigar

    pRealignedResult realigned_result = new ALLOC_TYPE_IN_POOL(pool, RealignedResult) RealignedResult{};

    uint32_t* original_cigar = original_read->cigar();
    uint32_t cigar_num = original_read->cigar_length();
    realigned_result->realigned_cigar = append_clipped_elements_from_cigar_to_cigar(lartrc, original_cigar, cigar_num, pool);

    realigned_result->new_start = (int64_t)read_start_on_reference + left_aligned_read_to_ref_cigar_result.leading_deletion_bases_removed;
    realigned_result->new_stop =
        realigned_result->new_start + bam_cigar2rlen(realigned_result->realigned_cigar->num, realigned_result->realigned_cigar->data) - 1;

    original_read->_realigned_result = realigned_result;
    CHECK_CONDITION_EXIT(bam_cigar2qlen((int32_t)lartrc->num, lartrc->data) + soft_clipped_bases != read_bases->num,
                         "create_read_aligned_to_ref boom!");
}

pBases AlignmentUtils::get_bases_covering_ref_interval(int64_t ref_start, int64_t ref_end, pBases bases, int64_t bases_start_on_ref,
                                                       pCigar bases_to_ref_cigar, pMemoryPool pool)
{
    CHECK_CONDITION_EXIT(ref_start < 0 || ref_end < ref_start, "bad start {} and/or stop {}", ref_start, ref_end);
    CHECK_CONDITION_EXIT(bases_start_on_ref < 0, "bases_start_on_ref must be >= 0");
    CHECK_CONDITION_EXIT(nullptr == bases || nullptr == bases_to_ref_cigar, "nullptr == bases || nullptr == bases_to_ref_cigar");
    CHECK_CONDITION_EXIT(bases->num != bam_cigar2qlen((int32_t)bases_to_ref_cigar->num, bases_to_ref_cigar->data),
                         "mismatch in length between reference and cigar");

    int64_t ref_pos = bases_start_on_ref;
    int64_t bases_pos = 0;
    int64_t bases_start = INVALID_INT;
    int64_t bases_stop = INVALID_INT;
    bool done = false;

    uint32_t ce, op, op_len;
    for (uint32_t iii = 0; !done && iii < bases_to_ref_cigar->num; iii++) {
        ce = bases_to_ref_cigar->data[iii];
        op = bam_cigar_op(ce);
        op_len = bam_cigar_oplen(ce);
        switch (op) {
            case BAM_CINS: {
                bases_pos += op_len;
                break;
            }
            case BAM_CDEL: {
                for (uint32_t i = 0; i < op_len; i++) {
                    if (ref_pos == ref_end || ref_pos == ref_start) {
                        return nullptr;
                    }
                    ref_pos++;
                }
                break;
            }
            case BAM_CMATCH:
            case BAM_CEQUAL:
            case BAM_CDIFF: {
                for (uint32_t i = 0; i < op_len; i++) {
                    if (ref_pos == ref_start) bases_start = bases_pos;
                    if (ref_pos == ref_end) {
                        bases_stop = bases_pos;
                        done = true;
                        break;
                    }
                    ref_pos++;
                    bases_pos++;
                }
                break;
            }
            case BAM_CREF_SKIP:
            case BAM_CSOFT_CLIP:
            case BAM_CHARD_CLIP:
            case BAM_CPAD:
            case BAM_CBACK:
            default: {
                RovacaLogger::error("unsupported operator: {}", op);
                exit(EXIT_FAILURE);
            }
        }
    }

    CHECK_CONDITION_EXIT(bases_start == INVALID_INT || bases_stop == INVALID_INT, "never found start or stop");

    char* start_ptr = (char*)bases->data + bases_start;
    uint32_t len = bases_stop - bases_start + 1;
    return BaseUtils::bases_create(start_ptr, len, pool);
}

std::pair<pBases, pBases> AlignmentUtils::get_bases_and_base_qualities_aligned_one_to_one(pReadRecord read, pMemoryPool pool)
{
    return get_bases_and_base_qualities_aligned_one_to_one(read, s_gap_base_character, s_gap_qual_character, pool);
}

std::pair<pBases, pBases> AlignmentUtils::get_bases_and_base_qualities_aligned_one_to_one(pReadRecord read, uint8_t base_pad,
                                                                                          uint8_t qual_pad, pMemoryPool pool)
{
    uint32_t* cigar = read->cigar();
    uint32_t cigar_num = read->cigar_length();
    bool saw_indel = false;
    for (uint32_t i = 0; i < cigar_num; ++i) {
        if (cigar_op_is_ins(bam_cigar_op(cigar[i])) || cigar_op_is_del(bam_cigar_op(cigar[i]))) {
            saw_indel = true;
            break;
        }
    }

    pBases bases, quals;
    if (!saw_indel) {
        uint32_t seq_len = uint32_t(read->seq_length());
        bases = new ALLOC_FLEXIBLE_IN_POOL(pool, Bases, seq_len, uint8_t) Bases{seq_len};
        quals = new ALLOC_FLEXIBLE_IN_POOL(pool, Bases, seq_len, uint8_t) Bases{seq_len};
        for (uint32_t i = 0; i < seq_len; ++i) {
            bases->data[i] = read->seq_i(i);
        }
        memcpy(quals->data, read->qual(), seq_len);
    }
    else {
        uint32_t seq_len = CigarUtils::count_ref_bases_and_soft_clips(read->cigar(), cigar_num, 0, cigar_num);
        bases = new ALLOC_FLEXIBLE_IN_POOL(pool, Bases, seq_len, uint8_t) Bases{seq_len};
        quals = new ALLOC_FLEXIBLE_IN_POOL(pool, Bases, seq_len, uint8_t) Bases{seq_len};
        uint32_t literal_pos = 0, padded_pos = 0;
        uint32_t ce, op, oplen;
        for (uint32_t i = 0; i < cigar_num; ++i) {
            ce = cigar[i];
            op = bam_cigar_op(ce);
            oplen = bam_cigar_oplen(ce);
            if (consumes_read_bases(op)) {
                if (!consumes_ref_bases(op)) {
                    literal_pos += oplen;  // skip inserted bases
                }
                else {
                    // note: 此处可能是瓶颈？
                    for (uint32_t j = 0; j < oplen; ++j) {
                        bases->data[padded_pos + j] = read->seq_i(literal_pos + j);
                        quals->data[padded_pos + j] = read->qual_i(literal_pos + j);
                    }
                    literal_pos += oplen;
                    padded_pos += oplen;
                }
            }
            else if (consumes_ref_bases(op)) {
                for (uint32_t j = 0; j < oplen; ++j) {
                    bases->data[padded_pos] = base_pad;
                    quals->data[padded_pos] = qual_pad;
                    padded_pos++;
                }
            }
        }
    }
    return {bases, quals};
}

}  // namespace rovaca
