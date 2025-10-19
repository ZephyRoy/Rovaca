#include "fake_assembler.h"

#include <algorithm>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

#include "../assemble_argument.h"
#include "assemble_engine.h"
#include "assemble_testcase_iterator.h"
#include "haplotype.h"
#include "hc_assemble_main.h"
#include "htslib/sam.h"
#include "read_record.h"
#include "simple_interval.h"
#include "testcase_iterator.h"
#include "testcase_loader.h"
#include "utils/base_utils.h"

static constexpr size_t s_buffer_size = 1024 * 1024 * 10;
static constexpr uint32_t k_default_region_padding = 100;

FakeAssembler::FakeAssembler(const char* casePath, const char* refPath, const char* resultPath, const AssembleArgument& argument,
                             pMemoryPool targetMem)
    : arguments_(argument)
    , assemble_debug_streamer_(resultPath, std::ofstream::out)
    , test_case_(casePath)
    , target_mem_(targetMem)
{
    arguments_.read_threading_argument.kmer = {10, 25};
    arguments_.debugAssembly = true;
    AssembleEngine::init_assemble_argument(&arguments_);
    std::ifstream ref_streamr(refPath);
    std::getline(ref_streamr, chr_ref);
    ref_streamr.close();
}

FakeAssembler::~FakeAssembler()
{
    AssembleEngine::finit_assemble_argument();
    assemble_debug_streamer_.close();
}

void FakeAssembler::run()
{
    for (auto& one_case : test_case_) {
        // Data prepair
        hc_region_active_storage region{
            tid : 0,
            active : 1,
            start_index : one_case.region_beg,
            end_index : one_case.region_end,
            activeSpan : {one_case.region_beg, one_case.region_end},
            paddedSpan : {one_case.region_beg - k_default_region_padding, one_case.region_end + k_default_region_padding}
        };
        std::vector<uint8_t> ref_data(one_case.reference.begin(), one_case.reference.end());
        const uint8_t* ref = reinterpret_cast<const uint8_t*>(chr_ref.data());
        chr_ref_len = chr_ref.size();
        std::vector<bam1_t*> region_reads = std::move(one_case.region_reads);

        // Call local asemble.
        // TODO:暂时new出，应指向RegionSource的缓存资源.
        p_hc_apply assembler = hc_apply_init();
        AssembleReadsBuffer reads_mem{buffer_ : new uint8_t[k_reads_buffer_capacity]{}, used_ : 0, capacity_ : k_reads_buffer_capacity};
        AssembleResult* untrimed_result =
            AssembleEngine::local_assemble(assembler, &region, ref, chr_ref_len, region_reads, target_mem_, &reads_mem);
        hc_apply_finit(assembler);
        // Get Assemble results.
        auto reads = untrimed_result->get_reads();
        auto haplotypes = untrimed_result->get_haplotypes();

        // Sort in Lexicographic order for convenient compair.
        std::sort(haplotypes.begin(), haplotypes.end(), [](rovaca::pHaplotype a, rovaca::pHaplotype b) {
            return strcmp((const char*)(a->get_display_string()->data), (const char*)(b->get_display_string()->data)) < 0;
        });
        if (assemble_debug_streamer_.good()) {
            assemble_debug_streamer_ << "chr1:" << region.activeSpan.start << "-" << region.activeSpan.end << "\t";
            assemble_debug_streamer_ << "reads: " << reads.size() << "\n";
            assemble_debug_streamer_ << "haplotypes: " << haplotypes.size() << std::endl;
            for (auto hap : haplotypes) {
                assemble_debug_streamer_ << hap->get_display_string()->data << "\n";
            }
            assemble_debug_streamer_ << std::endl;
            assemble_debug_streamer_.flush();
        }

        // 暂时释放要传下去的内存。
        untrimed_result->~AssembleResult();
        delete[] reads_mem.buffer_;
    }
}

int main(int argc, char* argv[])
{
    if (argc < 4) {
        std::cerr << "\nUsage:\t<assembler> [read_path] [chr_ref_path] [result_path]\n " << std::endl;
        return -1;
    }
    const char* dataPath = argv[1];
    const char* refPath = argv[2];
    const char* resultPath = argv[3];
    uint8_t* buffer = new uint8_t[s_buffer_size]{};
    pMemoryPool target_mem = new std::pmr::monotonic_buffer_resource(buffer, s_buffer_size, std::pmr::null_memory_resource());
    AssembleArgument argument = ASSEMBLE_DEFAULT_ARGUMENTS;
    FakeAssembler{dataPath, refPath, resultPath, argument, target_mem}.run();
    delete[] buffer;
    delete target_mem;
    return 0;
}