#include "genotype_allele_counts.h"

#include "rovaca_logger.h"
#include "genotype_macors.h"
#include "index_range.hpp"
#include "math_utils.h"

namespace rovaca
{

pGenotypeAlleleCounts GenotypeAlleleCounts::first(int32_t ploidy, pMemoryPool pool)
{
    CHECK_CONDITION_EXIT(ploidy < 0, "the ploidy must be 0 or greater");
    if (ROVACA_UNLIKELY(0 == ploidy)) {
        Int32Vector sorted_allele_counts(pool);
        return new ALLOC_TYPE_IN_POOL(pool, GenotypeAlleleCounts) GenotypeAlleleCounts{0, 0, sorted_allele_counts, 0, pool};
    }

    Int32Vector sorted_allele_counts{{0, ploidy}, pool};
    return new ALLOC_TYPE_IN_POOL(pool, GenotypeAlleleCounts) GenotypeAlleleCounts{0, ploidy, std::move(sorted_allele_counts)};
}

pGenotypeAlleleCounts GenotypeAlleleCounts::next(pMemoryPool pool)
{
    if (0 == _distinct_allele_count) {
        return this;
    }
    else if (1 == _distinct_allele_count && 1 == _ploidy) {
        Int32Vector sorted_allele_counts{{_sorted_allele_counts.at(0) + 1, 1}, pool};
        return new ALLOC_TYPE_IN_POOL(pool, GenotypeAlleleCounts) GenotypeAlleleCounts{_index + 1, 1, std::move(sorted_allele_counts)};
    }
    else if (1 == _distinct_allele_count) {
        Int32Vector sorted_allele_counts{{0, _ploidy - 1, _sorted_allele_counts.at(0) + 1, 1}, pool};
        return new ALLOC_TYPE_IN_POOL(pool, GenotypeAlleleCounts)
            GenotypeAlleleCounts{_index + 1, _ploidy, std::move(sorted_allele_counts)};
    }

    // the following logic avoids dynamically sizing the new sorted-allele-counts array, which would be very slow at this point
    // distinct_allele_count >= 2 thus sorted_allele_counts.length >= 4.
    // we only need to look at the two lowest allele indices to decide what to do.

    int32_t freq0 = _sorted_allele_counts.at(1);
    int32_t allele0plus1 = _sorted_allele_counts.at(0) + 1;
    bool allele0and1are_consecutive = allele0plus1 == _sorted_allele_counts[2];

    Int32Vector new_sorted_allele_counts(pool);
    int32_t sorted_allele_counts_length = _distinct_allele_count << 1;
    if (freq0 == 1) {
        // in this case allele0 won't be present in the result and all its frequency should go to allele0 + 1.
        if (allele0and1are_consecutive) {  // need just to remove the first allele and 1 to the frequency of the second (freq1 += 1).
            std::copy(_sorted_allele_counts.begin() + 2, _sorted_allele_counts.end(), std::back_inserter(new_sorted_allele_counts));
            new_sorted_allele_counts.at(1)++;
        }
        else {
            // just need to mutate allele0 to allele0 + 1.
            std::copy(_sorted_allele_counts.begin(), _sorted_allele_counts.end(), std::back_inserter(new_sorted_allele_counts));
            new_sorted_allele_counts[0] = allele0plus1;
            // new_sorted_allele_counts[1] = 1; // :) no need to do it because it is already the case (freq0 == 1).
        }
    }
    else {
        // && freq0 > 1 as per sorted_allele_counts format restrictions. in this case allele0 will muttated to '0' with frequency decreased
        // by 1.
        if (allele0and1are_consecutive) {  // we don't need to add a component for allele0 + 1 since it already exists.
            std::copy(_sorted_allele_counts.begin(), _sorted_allele_counts.end(), std::back_inserter(new_sorted_allele_counts));
            new_sorted_allele_counts[0] = 0;
            new_sorted_allele_counts[1] = freq0 - 1;
            new_sorted_allele_counts[3]++;
        }
        else {  // we need to insert allele0 + 1 in the sorted-allele-counts array.
            new_sorted_allele_counts.resize(sorted_allele_counts_length + 2);
            new_sorted_allele_counts[0] = 0;
            new_sorted_allele_counts[1] = freq0 - 1;
            new_sorted_allele_counts[2] = allele0plus1;
            new_sorted_allele_counts[3]++;  // = 1 as the array was freshly created with 0s.
            std::copy(_sorted_allele_counts.begin() + 2, _sorted_allele_counts.end(), std::back_inserter(new_sorted_allele_counts));
        }
    }

    return new ALLOC_TYPE_IN_POOL(pool, GenotypeAlleleCounts)
        GenotypeAlleleCounts{_index + 1, _ploidy, std::move(new_sorted_allele_counts)};
}

double GenotypeAlleleCounts::log10combination_count()
{
    if (_log10combination_count == NEGATIVE_INFINITY) {
        _log10combination_count = MathUtils::log10factorial(_ploidy) - IndexRange(0, _distinct_allele_count).sum([&](int32_t n) -> double {
            return MathUtils::log10factorial(_sorted_allele_counts.at(2 * n + 1));
        });
    }
    return _log10combination_count;
}

int32_t GenotypeAlleleCounts::allele_count_for(int32_t index) const
{
    int32_t rank = allele_rank_for(index);
    return rank < 0 ? 0 : allele_count_at(rank);
}

void GenotypeAlleleCounts::copy_allele_counts(int32_t offset, Int32Vector& dest) const
{
    int32_t sorted_allele_counts_length = _distinct_allele_count << 1;
    CHECK_CONDITION_EXIT(sorted_allele_counts_length + offset > (int32_t)dest.size(), "the input array does not have enough capacity");
    std::copy_n(_sorted_allele_counts.begin(), sorted_allele_counts_length, dest.begin() + (long)offset);
}

void GenotypeAlleleCounts::increase()
{
    if (0 == _distinct_allele_count) {
        return;
    }

    if (1 == _distinct_allele_count) {
        if (_ploidy == 1) {
            ++_sorted_allele_counts.at(0);
        }
        else {
            if (_sorted_allele_counts.size() < 4) {
                _sorted_allele_counts.resize(4);
            }
            _sorted_allele_counts[2] = _sorted_allele_counts[0] + 1;
            _sorted_allele_counts[3] = 1;
            _sorted_allele_counts[0] = 0;
            _sorted_allele_counts[1] = _ploidy - 1;
            _distinct_allele_count = 2;
        }
    }
    else {
        int32_t allele0 = _sorted_allele_counts[0];
        int32_t freq0 = _sorted_allele_counts[1];
        int32_t allele1 = _sorted_allele_counts[2];
        int32_t allele0plus1 = allele0 + 1;
        bool allele0and1are_consecutive = allele0plus1 == allele1;
        int32_t sorted_allele_counts_length = _distinct_allele_count << 1;

        if (freq0 == 1) {
            // in this case allele0 wont be present in the result and all is frequency should go to allele0 + 1
            if (allele0and1are_consecutive) {
                // need just to remove the first allele and add 1 to the frequency of the second (freq1 += 1)

                // shift left the first component away
                std::copy_n(_sorted_allele_counts.begin() + 2, sorted_allele_counts_length - 2, _sorted_allele_counts.begin());

                // freq1 has become freq0.
                _sorted_allele_counts[1]++;
                _distinct_allele_count--;
            }
            else {
                // just need to mutate allele0 to allele0 + 1
                _sorted_allele_counts[0] = allele0plus1;
            }
        }
        else {
            // && freq0 > 1 as per sorted_allele_counts format restrictions. in this case allele0 will mutated to '0' with frequency
            // decreased by 1.
            if (allele0and1are_consecutive) {
                // we don't need to add a component for allele0 + 1 since it already exists.
                _sorted_allele_counts[0] = 0;
                _sorted_allele_counts[1] = freq0 - 1;
                _sorted_allele_counts[3]++;
            }
            else {
                // we need to insert allele0 + 1 in the sorted-allele-counts array and give it frequency 1.
                if ((int32_t)_sorted_allele_counts.size() < sorted_allele_counts_length + 2) {
                    // make room for the new component.
                    _sorted_allele_counts.resize(sorted_allele_counts_length + 2);
                }
                std::copy_n(_sorted_allele_counts.begin() + 2, sorted_allele_counts_length - 2, _sorted_allele_counts.begin() + 4);
                _sorted_allele_counts[0] = 0;
                _sorted_allele_counts[1] = freq0 - 1;
                _sorted_allele_counts[2] = allele0plus1;
                _sorted_allele_counts[3] = 1;
                _distinct_allele_count++;
            }
        }
    }
    ++_index;
    _log10combination_count = NEGATIVE_INFINITY;
}

void GenotypeAlleleCounts::increase(int32_t times)
{
    CHECK_CONDITION_EXIT(times < 0, "times");
    for (int32_t i = 0; i < times; ++i) {
        increase();
    }
}

pGenotypeAlleleCounts GenotypeAlleleCounts::copy(pMemoryPool pool) const
{
    Int32Vector new_sorted_allele_counts(_sorted_allele_counts, pool);
    return new ALLOC_TYPE_IN_POOL(pool, GenotypeAlleleCounts)
        GenotypeAlleleCounts{_index, _ploidy, _distinct_allele_count, std::move(new_sorted_allele_counts)};
}

double GenotypeAlleleCounts::sum_over_allele_indices_and_counts(const std::function<double(int32_t, int32_t)>& func)
{
    return IndexRange(0, _distinct_allele_count).sum([&](int32_t n) -> double {
        return func(_sorted_allele_counts.at(2 * n), _sorted_allele_counts.at(2 * n + 1));
    });
}

void GenotypeAlleleCounts::for_each_allele_index_and_count(const std::function<void(int32_t, int32_t)>& func)
{
    IndexRange(0, _distinct_allele_count).for_each([&](int32_t n) {
        func(_sorted_allele_counts.at(2 * n), _sorted_allele_counts.at(2 * n + 1));
    });
}

void GenotypeAlleleCounts::for_each_absent_allele_index(const std::function<void(int32_t)>& func, int32_t allele_count)
{
    int32_t present_allele_index = 0;
    int32_t present_allele = _sorted_allele_counts.front();
    for (int32_t n = 0; n < allele_count; n++) {
        // if we find n in sorted_allele_counts, it is present, so we move present_allele to the next index in sorted_allele_counts and skip
        // the allele; otherwise the allele is absent and we perform the action on it.
        if (n == present_allele) {
            // if we haven't exhausted all the present alleles, move to the next one.
            // note that distinct_allele_count == sorted_allele_counts.length/2
            if (++present_allele_index < _distinct_allele_count) {
                // every other entry in sorted_allele_counts is an allele index; hence we multiply by 2
                present_allele = _sorted_allele_counts.at(2 * present_allele_index);
            }
            continue;
        }
        func(n);
    }
}

AlleleVector GenotypeAlleleCounts::as_allele_list(const AlleleVector& alleles_to_use)
{
    CHECK_CONDITION_EXIT(alleles_to_use.empty(), "the input allele list cannot be empty");
    CHECK_CONDITION_EXIT((int32_t)alleles_to_use.size() < maximum_allele_index(),
                         "the provided alleles to use does not contain an element for the maximum allele");
    if (1 == _distinct_allele_count) {
        return {(size_t)_ploidy, alleles_to_use.at(_sorted_allele_counts.at(0)), alleles_to_use.get_allocator()};
    }
    pAllele a;
    int32_t i, repeats;
    AlleleVector result(alleles_to_use.get_allocator());
    result.reserve(_distinct_allele_count * 2);
    for (i = 0; i < _distinct_allele_count; ++i) {
        a = alleles_to_use.at(_sorted_allele_counts.at(2 * i));
        repeats = _sorted_allele_counts.at(2 * i + 1);
        std::fill_n(std::back_inserter(result), repeats, a);
    }
    return result;
}

int32_t GenotypeAlleleCounts::allele_index_to_rank(int32_t index, int32_t from, int32_t to) const
{
    if (to <= from) {
        return -from - 1;
    }
    if (from == to - 1) {
        int32_t only_index = _sorted_allele_counts.at(from << 1);
        return only_index == index ? from : (only_index > index) ? -from - 1 : -to - 1;
    }

    int32_t mid = (to + from) >> 1;
    int32_t mid_index = _sorted_allele_counts.at(mid << 1);
    if (mid_index == index) {
        return mid;
    }
    else if (mid_index < index) {
        return allele_index_to_rank(index, mid + 1, to);
    }
    else {
        return allele_index_to_rank(index, 0, mid);
    }
}

Int32Vector GenotypeAlleleCounts::allele_counts_by_index(int32_t maximum_allele_index, pMemoryPool pool)
{
    CHECK_CONDITION_EXIT(maximum_allele_index < 0, "the requested allele count cannot be less than 0");
    Int32Vector result{size_t(maximum_allele_index + 1), pool};
    copy_allele_counts_by_index(result, 0, 0, maximum_allele_index);
    return result;
}

void GenotypeAlleleCounts::copy_allele_counts_by_index(Int32Vector& dest, int32_t offset, int32_t minimum_allele_index,
                                                       int32_t maximum_allele_index)
{
    // first we determine what section of the sorted_allele_counts array contains the counts of interest, by the present allele rank range
    // of interest.
    int32_t minimum_allele_rank = allele_rank_for(minimum_allele_index);
    int32_t maximum_allele_rank = allele_rank_for(maximum_allele_index);

    // if the min or max allele index are absent (returned rank < 0) we note where the would be inserted; that way we avoid going through
    // the rest of positions in the sorted_allele_counts array. the range of interest is then [start_rank,end_rank].
    int32_t start_rank = minimum_allele_rank < 0 ? -minimum_allele_rank - 1 : minimum_allele_rank;
    int32_t end_rank = maximum_allele_rank < 0 ? -maximum_allele_rank - 2 : maximum_allele_rank;

    // iteration variables:
    int32_t next_index = minimum_allele_index;                  // next index that we want to output the count for.
    int32_t next_rank = start_rank;                             // next rank to query in sorted_allele_counts.
    int32_t next_sorted_allele_counts_offset = next_rank << 1;  // offset in sorted_allele_counts where the info is present for the nextrank
    int32_t next_dest_offset = offset;  // next offset in destination array where to set the count for the next_index.

    while (next_rank++ <= end_rank) {
        int32_t allele_index = _sorted_allele_counts[next_sorted_allele_counts_offset++];
        // fill non-present allele counts with 0s.
        while (allele_index > next_index) {
            dest[next_dest_offset++] = 0;
            next_index++;
        }
        // it is guaranteed that at this point allele_index == next_index
        // thanks to the condition of the enclosing while: there must be at least one index of interest that
        // is present in the remaining (next_rank,end_rank] interval as otherwise end_rank would be less than next_rank.
        dest[next_dest_offset++] = _sorted_allele_counts[next_sorted_allele_counts_offset++];
        next_index++;
    }
    // finally we take care of trailing requested allele indices.
    while (next_index++ <= maximum_allele_index) {
        dest[next_dest_offset++] = 0;
    }
}

}  // namespace rovaca