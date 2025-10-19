#include <gtest/gtest.h>

#include <algorithm>
#include <string>
#include <tuple>

#include "allele.h"
#include "event_map.h"
#include "forward.h"
#include "haplotype.h"
#include "simple_interval.h"
#include "utils/base_utils.h"
#include "utils/cigar_utils.h"
#include "variant.h"

using namespace std;
using namespace rovaca;

static constexpr size_t s_buffer_size = 1024 * 1024 * 10;

class EventMapUnitTest : public ::testing::Test
{
protected:
    void SetUp() override
    {
        _buffer = new uint8_t[s_buffer_size]{};
        _pool = new std::pmr::monotonic_buffer_resource(_buffer, s_buffer_size, std::pmr::null_memory_resource());
    }

    void TearDown() override
    {
        delete _pool;
        delete[] _buffer;
    }

    pVariant make_from_alleles(int32_t tid, int64_t start, const pair<string, string>& allele_str);

    vector<tuple<string, string, int64_t, pAllele, pAllele>> test_get_overlapping_events_data();

    uint8_t* _buffer{};
    pMemoryPool _pool{};
};

TEST_F(EventMapUnitTest, testGetNeighborhood)
{
    vector<pair<string, string>> first{{"A", "G"}, {"A", "G"}, {"AC", "A"}, {"ACGTA", "A"}, {"AC", "A"}, {"A", "ACGTA"}, {"A", "AC"}};
    vector<pair<string, string>> second{{"AGT", "A"}, {"A", "AGT"}, {"A", "AGT"}, {"A", "AG"}, {"A", "AGCGT"}, {"AG", "A"}, {"AGCGT", "A"}};
    vector<pair<string, string>> expected{{"AGT", "G"},    {"A", "GGT"},    {"AC", "AGT"},  {"ACGTA", "AG"},
                                          {"AC", "AGCGT"}, {"AG", "ACGTA"}, {"AGCGT", "AC"}};
    pVariant vc_1, vc_2, vc_expected, block;
    pEventMap event;
    for (size_t i = 0, len = first.size(); i < len; ++i) {
        const pair<string, string>& f1 = first.at(i);
        const pair<string, string>& s2 = second.at(i);
        const pair<string, string>& e3 = expected.at(i);
        vc_1 = make_from_alleles(20, 10, f1);
        vc_2 = make_from_alleles(20, 10, s2);
        vc_expected = make_from_alleles(20, 10, e3);
        event = new ALLOC_TYPE_IN_POOL(_pool, EventMap) EventMap{0, _pool};
        block = event->make_block(vc_1, vc_2);

        ASSERT_EQ(block->get_start(), vc_expected->get_start());

        const AlleleVector& block_alleles = block->alleles();
        const AlleleVector& expected_alleles = vc_expected->alleles();
        for (size_t j = 0; j < 2; ++j) {
            pAllele block_a = block_alleles.at(j);
            pAllele expected_a = expected_alleles.at(j);
            ASSERT_TRUE(block_a->equals(*expected_a)) << "i:" << i << " j:" << j << " block_a:" << block_a->get_display_string()->data
                                                      << " expected_a:" << expected_a->get_display_string()->data;
        }
    }
}

TEST_F(EventMapUnitTest, testGetOverlappingEvents)
{
    string ref_bases = "AAAAAAAAAACGGTCA";
    int64_t hap_start_wrt_ref = 7;
    pSimpleInterval ref_loc = SimpleInterval::create(0, 1, (int64_t)ref_bases.size(), _pool);
    vector<tuple<string, string, int64_t, pAllele, pAllele>> data = test_get_overlapping_events_data();
    RefFragment ref{(uint32_t)ref_bases.size(), (uint8_t*)ref_bases.c_str()};

    for (size_t i = 0, len = data.size(); i < len; ++i) {
        const tuple<string, string, int64_t, pAllele, pAllele>& tup = data.at(i);
        const string& haplotype_bases = get<0>(tup);
        const string& cigar_str = get<1>(tup);
        int64_t query_loc = get<2>(tup);
        pAllele expected_ref = get<3>(tup);
        pAllele expected_alt = get<4>(tup);

        uint32_t cigar_num = std::count(cigar_str.begin(), cigar_str.end(), 'M') + std::count(cigar_str.begin(), cigar_str.end(), 'D') +
                             std::count(cigar_str.begin(), cigar_str.end(), 'I');
        pCigar cigar = CigarUtils::str2uint(cigar_str.c_str(), cigar_num, _pool);

        pHaplotype h = Haplotype::create(_pool);
        h->init_haplotype(haplotype_bases.c_str(), haplotype_bases.size(), 0, _pool);
        h->set_alignment_start_hap_wrt_ref(hap_start_wrt_ref);
        h->set_cigar(cigar);

        HaplotypeVector haplotypes{{h}, _pool};

        EventMap::build_event_maps_for_haplotypes(haplotypes, &ref, ref_loc, 1, _pool);

        pEventMap events = h->event_map();
        VariantVector overlapping_events = events->get_overlapping_events(query_loc);

        bool events_expected = expected_alt != nullptr || expected_ref != nullptr;
        ASSERT_EQ(overlapping_events.size(), events_expected ? 1 : 0) << "i=" << i;

        if (events_expected) {
            ASSERT_TRUE(overlapping_events.at(0)->ref_allele()->equals(*expected_ref)) << "i=" << i;
            ASSERT_TRUE(overlapping_events.at(0)->alternate_allele_at(0)->equals(*expected_alt)) << "i=" << i;
        }
    }
}

pVariant EventMapUnitTest::make_from_alleles(int32_t tid, int64_t start, const pair<string, string>& allele_str)
{
    pBases ref_bases = BaseUtils::bases_create(allele_str.first.c_str(), _pool);
    pAllele ref_allele = Allele::create_allele(ref_bases, 1, _pool);
    pBases alt_bases = BaseUtils::bases_create(allele_str.second.c_str(), _pool);
    pAllele alt_allele = Allele::create_allele(alt_bases, 0, _pool);
    AlleleVector alleles{{ref_allele, alt_allele}, _pool};

    int64_t stop = start + ref_allele->length() - 1;

    pVariant result = Variant::create(_pool);
    result->set_tid(tid);
    result->set_start(start);
    result->set_stop(stop);
    result->set_alleles(alleles);
    return result;
}

vector<tuple<string, string, int64_t, pAllele, pAllele>> EventMapUnitTest::test_get_overlapping_events_data()
{
    pAllele deletion_ref_allele = Allele::create_allele("ACGG", true, _pool);
    pAllele deletion_alt_allele = Allele::create_allele("A", false, _pool);
    pAllele insertion_ref_allele = Allele::create_allele("G", true, _pool);
    pAllele insertion_alt_allele = Allele::create_allele("GTT", false, _pool);
    pAllele snp_ref_allele = Allele::create_allele("G", true, _pool);
    pAllele snp_alt_allele = Allele::create_allele("A", false, _pool);

    vector<tuple<string, string, int64_t, pAllele, pAllele>> data{
        // hap1
        {"AAATTTCA", "3M3D2I3M", 10, deletion_ref_allele, deletion_alt_allele},
        {"AAATTTCA", "3M3D2I3M", 11, deletion_ref_allele, deletion_alt_allele},
        {"AAATTTCA", "3M3D2I3M", 12, deletion_ref_allele, deletion_alt_allele},
        {"AAATTTCA", "3M3D2I3M", 13, insertion_ref_allele, insertion_alt_allele},
        // hap2
        {"AAATCA", "3M3D3M", 10, deletion_ref_allele, deletion_alt_allele},
        {"AAATCA", "3M3D3M", 11, deletion_ref_allele, deletion_alt_allele},
        {"AAATCA", "3M3D3M", 12, deletion_ref_allele, deletion_alt_allele},
        {"AAATCA", "3M3D3M", 13, deletion_ref_allele, deletion_alt_allele},
        // hap3
        {"AAACGATCA", "9M", 10, nullptr, nullptr},
        {"AAACGATCA", "9M", 11, nullptr, nullptr},
        {"AAACGATCA", "9M", 12, nullptr, nullptr},
        {"AAACGATCA", "9M", 13, snp_ref_allele, snp_alt_allele},
        // hap4
        {"AAACGGTTTCA", "6M2I3M", 10, nullptr, nullptr},
        {"AAACGGTTTCA", "6M2I3M", 11, nullptr, nullptr},
        {"AAACGGTTTCA", "6M2I3M", 12, nullptr, nullptr},
        {"AAACGGTTTCA", "6M2I3M", 13, insertion_ref_allele, insertion_alt_allele}};

    return data;
}
