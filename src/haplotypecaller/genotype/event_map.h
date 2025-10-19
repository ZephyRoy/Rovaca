#ifndef ROVACA_HC_EVENT_MAP_H_
#define ROVACA_HC_EVENT_MAP_H_
#include <cstdint>

#include "forward.h"
#include "genotype_macors.h"

namespace rovaca
{

class EventMap
{
    /* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
private:
    pMemoryPool _pool;
    Int64ToVariantMap _events;
    int32_t _source_index;

    /* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
public:
    EventMap(int32_t source, pMemoryPool pool)
        : _pool(pool)
        , _events(pool)
        , _source_index(source)
    {}
    ~EventMap() = default;
    DISALLOW_COPY_AND_ASSIGN(EventMap);

    /*!
     * @brief Build event maps for each haplotype, returning the sorted set of all of the starting positions of all
     * events across all haplotypes
     */
    static Int64Set build_event_maps_for_haplotypes(HaplotypeVector& haplotypes, pRefFragment ref, pSimpleInterval ref_loc,
                                                    uint32_t max_mnp_distance, pMemoryPool pool);

    /*! @brief returns any events in the map that overlap loc, including spanning deletions and events that start at loc */
    VariantVector get_overlapping_events(int64_t loc) const;

    const Int64ToVariantMap& get_events() const { return _events; }

    /**
     * Create a block substitution out of two variant contexts that start at the same position
     *
     * vc1 can be SNP, and vc2 can then be either a insertion or deletion.
     * If vc1 is an indel, then vc2 must be the opposite type (vc1 deletion => vc2 must be an insertion)
     *
     * @param vc1 the first variant context we want to merge
     * @param vc2 the second
     * @return a block substitution that represents the composite substitution implied by vc1 and vc2
     */
    pVariant make_block(pVariant vc1, pVariant vc2);

    /* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
private:
    /*!
     * @param max_mnp_distance Phased substitutions separated by this distance or less are merged into MNPs.  More than two substitutions
     * occurring in the same alignment block (ie the same M/X/EQ CIGAR element) are merged until a substitution is separated from the
     * previous one by a greater distance. That is, if maxMnpDistance = 1, substitutions at 10,11,12,14,15,17 are partitioned into a MNP at
     * 10-12, a MNP at 14-15, and a SNP at 17.  May not be negative.
     */
    void process_cigar_for_initial_events(uint32_t max_mnp_distance, pRefFragment ref, pSimpleInterval ref_loc, pHaplotype h);

    /*! @brief add variant to this map, merging events with the same start sites if necessary */
    void add_variant(pVariant vc);
};

}  // namespace rovaca

#endif  // ROVACA_HC_EVENT_MAP_H_
