#ifndef ROVACA_PLUGIN_H
#define ROVACA_PLUGIN_H

#include <functional>
#include <memory>
#include <sstream>
#include <string>
#include <unordered_map>
#include <vector>

#include "bam_loader.h"
#include "bed_loader.h"
#include "downsampler_hc.h"
#include "fasta_loader.h"
#include "rovaca_tool_args.h"
#include "reads_filter_hc.h"
#include "reads_stream.h"
#include "transformer.h"

using UniqueStream = std::unique_ptr<ReadStream>;
using RGHashMap = std::unordered_map<std::string, std::string>;
static htsFile *k_outfile = nullptr;
class RovacaTool
{
public:
    RovacaTool();
    virtual ~RovacaTool();
    virtual void do_work() = 0;
    virtual void clear_and_exit() = 0;
    bool initialize_args(int argc, char *argv[]);
    void run();

private:
    bool start_up();
    bool initialize_stream_pool();
    bool initialize_bam();
    bool initialize_regions();
    bool initialize_reference();
    bool matched_reference();
    virtual const UniqueStream &custom_streamer() = 0;
    const UniqueStream &make_streamer();

public:
    std::unique_ptr<RovacaToolArgs> rovaca_args_;
    std::unique_ptr<htsThreadPool> stream_pool_;
    std::unique_ptr<BamLoader> bam_loader_;
    std::unique_ptr<BedLoader> bed_loader_;
    UniqueStream streamer_;
    std::unique_ptr<ReadFilter> filter_;
    std::unique_ptr<Transformer> transformer_;
    std::unique_ptr<Downsampler> downsampler_;
    std::map<std::string, p_bed_intervals> bed_intervals_;
    p_bed_intervals all_intervals_;
    p_bed_intervals all_padded_intervals_;
    std::vector<std::pair<std::string, int64_t>> bam_contigs_;
    std::vector<RGHashMap> read_group_map_;
    std::vector<std::string> chr_list_;
    std::vector<bam_hdr_t *> bam_headers_;
    std::vector<std::string> samplelist_;
    std::vector<std::string> read_group_id_;
    contig_info_t fast_dict_;
    bam_hdr_t *merged_header_;
    bool wes_;
    const char *user_interval_;
    std::string read_path_;
    std::string interval_path_;
    std::string reference_path_;
    bool apply_filter_;
    bool apply_transformer_;
    bool apply_downsampler_;
    uint32_t downsample_threshold_;
    uint32_t min_base_quality_;
    uint32_t interval_padding_;
};
RovacaTool::RovacaTool()
    : rovaca_args_(nullptr)
    , stream_pool_(nullptr)
    , bam_loader_(nullptr)
    , bed_loader_(nullptr)
    , streamer_(nullptr)
    , filter_(nullptr)
    , transformer_(nullptr)
    , downsampler_(nullptr)
    , bed_intervals_({})
    , all_intervals_(nullptr)
    , all_padded_intervals_(nullptr)
    , bam_contigs_({})
    , chr_list_({})
    , bam_headers_({})
    , samplelist_({})
    , fast_dict_({})
    , wes_(false)
    , user_interval_(nullptr)
    , read_path_("")
    , interval_path_("")
    , reference_path_("")
    , apply_filter_(true)
    , apply_transformer_(true)
    , apply_downsampler_(true)
    , downsample_threshold_(0)
    , min_base_quality_(0)
    , interval_padding_(0)
{}

RovacaTool::~RovacaTool()
{
    bam_loader_.reset();
    filter_.reset();
    hts_tpool_destroy(stream_pool_->pool);
    RovacaLogger::info("Rovaca-{} finished", rovaca_args_->tool());
}

bool RovacaTool::initialize_args(int argc, char *argv[])
{
    rovaca_args_ = std::make_unique<RovacaToolArgs>(argc, argv);
    return rovaca_args_->valid_check();
}

bool RovacaTool::start_up()
{
    if (!initialize_stream_pool()) {
        RovacaLogger::error("failed to initialize stream pool");
        return false;
    }
    if (!initialize_bam()) {
        RovacaLogger::error("failed to initialize bam file");
        return false;
    }
    if (!initialize_reference()) {
        RovacaLogger::error("failed to initialize fasta file");
        return false;
    }
    if (!matched_reference()) {
        RovacaLogger::error("unmatched bam and reference");
        return false;
    }
    // Reference must be initialized before Bed.
    if (!initialize_regions()) {
        RovacaLogger::error("failed to initialize bed file");
        return false;
    }
    const auto &streamer = make_streamer();
    return streamer != nullptr;
}

void RovacaTool::run()
{
    if (start_up()) {
        do_work();
    }
}

bool RovacaTool::initialize_stream_pool()
{
    const uint32_t wes_stream_pool_size = 4;
    wes_ = !interval_path_.empty();
    stream_pool_ = std::make_unique<htsThreadPool>();
    stream_pool_->pool = hts_tpool_init(wes_ ? wes_stream_pool_size : rovaca_args_->stream_pool_size());
    stream_pool_->qsize = 0;
    return stream_pool_->pool != nullptr;
}

bool RovacaTool::initialize_bam()
{
    read_path_ = rovaca_args_->bam_path()[0];
    bam_loader_ = std::make_unique<BamLoader>(read_path_.c_str(), stream_pool_.get(), wes_);
    // Get bamhdr from bam.
    bam_headers_ = bam_loader_->get_sam_hdr();
    // TODO: For multiple bam files, header need to be merged.
    merged_header_ = bam_headers_[0];
    auto hdr = merged_header_;
    for (int i = 0; i < hdr->n_targets; ++i) {
        std::string chr_name = *(hdr->target_name + i);
        uint32_t chr_len = *(hdr->target_len + i);
        bam_contigs_.push_back({chr_name, chr_len});
    }

    // Extract all RG info from bamhdr.
    const char *p_rg = strstr(hdr->text, "@RG");
    while (p_rg) {
        RGHashMap read_group;
        std::string line = std::string(p_rg, strstr(p_rg + 3, "\n") - p_rg);
        std::istringstream iss(line);
        std::string token;
        while (iss >> token) {
            size_t pos = token.find(':');
            if (pos != std::string::npos) {
                std::string key = token.substr(0, pos);
                std::string value = token.substr(pos + 1);
                read_group.insert({key, value});
            }
        }
        read_group_map_.emplace_back(read_group);
        p_rg = strstr(p_rg + 1, "@RG");
    }

    // Get sample list.
    for (auto &it : read_group_map_) {
        if (it.count("SM")) samplelist_.emplace_back(it["SM"]);
        if (it.count("ID")) read_group_id_.emplace_back(it["ID"]);
    }
    return !(bam_loader_ == nullptr || merged_header_ == nullptr);
}

bool RovacaTool::initialize_regions()
{
    interval_path_ = rovaca_args_->bed_path();
    interval_padding_ = rovaca_args_->interval_padding();
    wes_ = !interval_path_.empty();
    bed_loader_ = wes_ ? std::make_unique<BedLoader>(interval_path_, interval_padding_, fast_dict_) : nullptr;
    if (bed_loader_) {
        bed_intervals_ = bed_loader_->get_bed_intervals();
        chr_list_ = bed_loader_->get_bed_keys();
        all_intervals_ = bed_loader_->get_all_intervals();
        all_padded_intervals_ = bed_loader_->get_all_padded_intervals();
    }
    else {
        chr_list_ = fast_dict_.key;
    }
    bool ret = true;
    // check whether -L interval beyond chrome length.
    for (const auto &pair : bed_intervals_) {
        const std::string &chr = pair.first;
        p_bed_intervals intervals = pair.second;
        uint32_t tid = sam_hdr_name2tid(merged_header_, chr.c_str());
        uint64_t ref_len = fast_dict_.idict[tid];
        for (int i = 0; i < intervals->n; ++i) {
            long int interval_len = intervals->end[i] - intervals->start[i];

            if (interval_len < 0) {
                RovacaLogger::error("interval {} get a negative length", intervals->regarry[i]);
                ret = false;
            }
            if (interval_len > static_cast<long int>(ref_len)) {
                RovacaLogger::error("interval {} is beyond ref length {}", intervals->regarry[i], ref_len);
                ret = false;
            }
        }
    }
    return !chr_list_.empty() && ret;
}

bool RovacaTool::initialize_reference()
{
    reference_path_ = rovaca_args_->reference_path();
    FastaLoader::get_fasta_dict(reference_path_, &fast_dict_);
    return !(fast_dict_.dict.empty() || fast_dict_.key.empty());
}

bool RovacaTool::matched_reference()
{
    if (bam_contigs_.size() != fast_dict_.dict.size()) {
        return false;
    }
    for (size_t i = 0; i < bam_contigs_.size(); ++i) {
        if (bam_contigs_[i].first != fast_dict_.key[i] ||
            bam_contigs_[i].second != static_cast<int64_t>(fast_dict_.dict[fast_dict_.key[i]])) {
            return false;
        }
    }
    return true;
}

// set_target() supports two modes: complete WES reading and specified chromosome/region reading. WGS does not need to be set.
const UniqueStream &RovacaTool::make_streamer()
{
    downsample_threshold_ = rovaca_args_->max_reads_depth();
    apply_downsampler_ = downsample_threshold_ != 0;
    apply_transformer_ = !rovaca_args_->recal_table().empty();
    streamer_ = std::make_unique<ReadStream>(bam_loader_.get());

    user_interval_ = rovaca_args_->target_span();
    const auto &streamer = custom_streamer();
    if (user_interval_ && !wes_) streamer->set_target(user_interval_);
    if (wes_) streamer->set_target(all_padded_intervals_);
    return streamer;
}

class RovacaToolRegister
{
public:
    using Creator = std::function<std::unique_ptr<RovacaTool>()>;
    static RovacaToolRegister &instance()
    {
        static RovacaToolRegister instance;
        return instance;
    }

    bool register_tool(const std::string &tool_name, Creator creator)
    {
        creators[tool_name] = std::move(creator);
        return !creators.empty();
    }

    std::unique_ptr<RovacaTool> creat_tool(const std::string &tool_name)
    {
        const auto &it = creators.find(tool_name);
        if (it != creators.end()) {
            RovacaLogger::info("Rovaca-{} launched", tool_name);
            return it->second();
        }
        return nullptr;
    }

    int supported_tools()
    {
        std::string all_tools;
        for (const auto &it : creators) {
            all_tools.append(it.first).append(",");
        }
        all_tools.back() = '\0';
        RovacaLogger::info("supported tools: {}", all_tools);

        return 0;
    }

private:
    RovacaToolRegister() = default;
    std::map<std::string, Creator> creators;
};

#endif  // Rovaca_PLUGIN_H