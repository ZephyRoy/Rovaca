#include <getopt.h>

#include <boost/asio.hpp>
#include <boost/asio/thread_pool.hpp>
#include <boost/lockfree/queue.hpp>
#include <iostream>
#include <memory>
#include <string>
#include <thread>

#include "ActiveMainThread.h"
#include "BlockingQueue.h"
#include "downsampler_hc.h"
#include "fasta_loader.h"
#include "reads_filter_hc.h"
#include "reads_stream_builder.h"
#include "reference_manager.h"

int usage()
{
    std::cerr << "ActiveRegionTest region" << std::endl;
    std::cerr << "--bam <FILE> bam file path" << std::endl;
    std::cerr << "--bed <FILE> bed file path" << std::endl;
    std::cerr << "--thread <INT>  number of threads " << std::endl;
    std::cerr << "--reference <FILE> referecen file path" << std::endl;
    std::cerr << "--output <FILE> region output file path" << std::endl;
    std::cerr << "--force  output nonActive Region" << std::endl;
    return -1;
}

const int WES_BUFFER = 2 * 1024 * 1024;
const int WES_BAM_POOL = 4 * 1024 * 1024;
const int WES_ACTIVE_BASE = 16 * 1024 * 1024;

const int WGS_BUFFER = 32 * 1024 * 1024;
const int WGS_BAM_POOL = 64 * 1024 * 1024;
const int WGS_ACTIVE_BASE = 1 * 1024 * 1024;
const int WGS_TARGET_MEM = 1;

int main(int argc, char *argv[])
{
    struct option long_options[]{{"bam", required_argument, 0, 'i'},
                                 {"bed", required_argument, 0, 'b'},
                                 {"reference", required_argument, 0, 'R'},
                                 {"thread", required_argument, 0, 't'},
                                 {"help", no_argument, 0, 'h'},
                                 {"output", required_argument, 0, 'o'},
                                 {"force", no_argument, 0, 'f'},
                                 {"padding", required_argument, 0, 'p'},
                                 {"max-reads-per-alignment-start", required_argument, 0, 'm'}};
    const char *options = "fo:t:i:b:R:m:p:h";
    int c = 0;
    int option_index = 0;
    int n_threads = 4;
    std::string bedfilename;
    std::string bam_filename;
    std::string reference;
    std::string outfilename;
    BedLoader *bed_loader = nullptr;
    bool force = false;
    int max_reads_per_alignment_start = 0;
    int padding = 0;

    while ((c = getopt_long(argc, argv, options, long_options, &option_index)) >= 0) {
        switch (c) {
            case 'h' /* constant-expression */: return usage(); break;
            case 'R': reference = std::string(optarg); break;
            case 'i': bam_filename = std::string(optarg); break;
            case 'b': bedfilename = std::string(optarg); break;
            case 't': n_threads = atoi(optarg); break;
            case 'o': outfilename = std::string(optarg); break;
            case 'f': force = true; break;
            case 'p': padding = atoi(optarg); break;
            case 'm': max_reads_per_alignment_start = atoi(optarg); break;
            default: break;
        }
    }
    int n_resource = n_threads * 2;
    if (reference.empty() || bam_filename.empty()) return usage();
    contig_info_t fasta_info;
    FastaLoader::get_fasta_dict(reference, &fasta_info);

    if (!bedfilename.empty()) {
        bed_loader = new BedLoader(bedfilename, padding, fasta_info);
        std::cerr << "padding=" << padding << std::endl;
        n_resource = n_threads * 4;
    }
    ReferenceManager fasta_manager_inst(reference, bed_loader);
    htsThreadPool tpool{nullptr, 0};
    tpool.pool = hts_tpool_init(bedfilename.empty() ? 4 : 10);
    BamLoader bam_info{bam_filename.c_str(), &tpool, !bedfilename.empty()};
    ReadFilter *reads_filter = new HCReadFilter(bam_info.get_sam_hdr()[0], false);
    ReadStreamBuilder builder(&bam_info);
    ActiveRegionBamBlockListSource block_resource(10);
    Downsampler *reads_downsampler = nullptr;
    if (max_reads_per_alignment_start != 0) {
        reads_downsampler = new HCDownsampler(max_reads_per_alignment_start);
    }
    ReadStream reads = builder.set_filter(reads_filter);

    if (bedfilename.empty()) {
        // bam_info.set_target("chr1:535879-535989");
    }
    else {
        // const std::map<std::string, p_bed_intervals> &all_interval = bed_loader->get_bed_intervals();
        // auto it = all_interval.find(std::string("chr1"));
        // p_bed_intervals current_interval = it->second;
        bam_info.set_target(bed_loader->get_all_intervals());
    }

    reads.init_reads_streamr();

    BlockingQueue<std::shared_ptr<ActiveBaseResource>> *base_resource = new BlockingQueue<std::shared_ptr<ActiveBaseResource>>(n_resource);

    BlockingQueue<BaseSource> *base_source = new BlockingQueue<BaseSource>(n_resource);
    BlockingQueue<std::shared_ptr<RegionResource>> *region_source = new BlockingQueue<std::shared_ptr<RegionResource>>(n_resource);
    BlockingQueue<std::shared_ptr<RegionSource>> *region_queue = new BlockingQueue<std::shared_ptr<RegionSource>>(n_resource);
    BlockingQueue<std::shared_ptr<BamSource>> *m_bam_resource = new BlockingQueue<std::shared_ptr<BamSource>>(n_resource + 20);

    for (int i = 0; i < n_resource; i++) {
        std::shared_ptr<ActiveBaseResource> active_base_resource =
            std::make_shared<ActiveBaseResource>(1024, bedfilename.empty() ? WGS_ACTIVE_BASE : WES_ACTIVE_BASE);
        std::shared_ptr<BamSource> bam_source = std::make_shared<BamSource>();
        std::shared_ptr<RegionResource> region_resource = std::make_shared<RegionResource>(
            i, bedfilename.empty() ? WGS_BAM_POOL : WES_BAM_POOL, bedfilename.empty() ? WGS_BUFFER : WES_BUFFER, WGS_TARGET_MEM);
        m_bam_resource->push(bam_source);
        base_resource->push(active_base_resource);
        region_source->push(region_resource);
    }
    for (int i = 0; i < 10; i++) {
        std::shared_ptr<BamSource> bam_source = std::make_shared<BamSource>();
        m_bam_resource->push(bam_source);
    }
    boost::asio::thread_pool thread_pool(n_threads);

    std::thread fasta_process(&ReferenceManager::run, &fasta_manager_inst);

    ActiveMainThreadDispatchTasks main_dispatch_process(&reads, bed_loader, &fasta_manager_inst, &fasta_manager_inst.get_contig(),
                                                        &block_resource, &thread_pool, base_resource, base_source, m_bam_resource);

    FakeConsumer fake_process(region_source, region_queue, outfilename, bam_info.get_sam_hdr()[0], &fasta_manager_inst);

    ActiveMainThreadReduce main_reduce_process(&fasta_manager_inst, &fasta_manager_inst.get_contig(), bed_loader, region_source,
                                               base_resource, &block_resource, base_source, region_queue, force);

    std::thread thread_dispatch(&ActiveMainThreadDispatchTasks::run, &main_dispatch_process);

    std::thread thread_reduce(&ActiveMainThreadReduce::run, &main_reduce_process);

    std::thread thread_fake(&FakeConsumer::run, &fake_process);

    pthread_setname_np(fasta_process.native_handle(), "Reference");
    pthread_setname_np(thread_dispatch.native_handle(), "RegionDispatch");
    pthread_setname_np(thread_reduce.native_handle(), "RegionReduce");
    pthread_setname_np(thread_fake.native_handle(), "FakeConsumer");

    if (max_reads_per_alignment_start != 0) {
        delete reads_downsampler;
    }
    fasta_process.join();
    thread_dispatch.join();
    thread_reduce.join();
    thread_fake.join();
    thread_pool.join();

    delete base_resource;
    delete base_source;
    delete region_source;
    delete region_queue;

    delete reads_filter;
    delete reads_downsampler;

    delete bed_loader;
    delete m_bam_resource;
    // delete main_reduce_process;
    return 0;
}