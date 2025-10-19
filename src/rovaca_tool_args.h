#ifndef ROVACA_TOOL_ARGS_H
#define ROVACA_TOOL_ARGS_H

#include <boost/program_options.hpp>
#include <cstdint>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

#include "common/valid_file.h"
#include "haplotypecaller/common/enum.h"
#include "rovaca_logger/rovaca_logger.h"
#include "version.h"

namespace po = boost::program_options;
class RovacaToolArgs
{
public:
    struct AssembleArgument
    {};
    struct DbsnpArgument
    {};
    struct LikelihoodArgument
    {};
    struct PileupDetectionArgument
    {};
    struct SWParameters
    {
        int matchValue;
        int mismatchPenalty;
        int gapOpenPenalty;
        int gapExtendPenalty;
    };

    enum WriterType { ALL_POSSIBLE_HAPLOTYPES, CALLED_HAPLOTYPES, NO_HAPLOTYPES, CALLED_HAPLOTYPES_NO_READS };
    enum SmithWatermanAlignerType { NORMAL, INTEL_AVX };
    enum LogLevel { ERROR, WARNING, INFO, DEBUG };

    SWParameters DEFAULT_DANGLING_END_SMITH_WATERMAN_PARAMETERS{25, -50, -110, -6};
    SWParameters DEFAULT_HAPLOTYPE_TO_REFERENCE_SMITH_WATERMAN_PARAMETERS{200, -150, -260, -11};
    SWParameters DEFAULT_READ_TO_HAPLOTYPE_SMITH_WATERMAN_PARAMETERS{10, -15, -30, -5};

    float PREFILTER_QUAL_THRESHOLD{30};
    float PREFILTER_SOR_THRESHOLD{3};

    std::string bamOutputPath;
    std::string alleleLikelihoodMatrixPath;
    std::string alleleLikelihoodMatrixInterval;
    bool dontUseSoftClippedBases{false};
    bool overrideSoftclipFragmentCheck{false};
    uint8_t minBaseQualityScore{10};

    WriterType bamWriterType{CALLED_HAPLOTYPES};
    SmithWatermanAlignerType smithWatermanImplementation{INTEL_AVX};
    SWParameters getDanglingEndSWParameters() { return DEFAULT_DANGLING_END_SMITH_WATERMAN_PARAMETERS; }
    SWParameters getHaplotypeToReferenceSWParameters() { return DEFAULT_HAPLOTYPE_TO_REFERENCE_SMITH_WATERMAN_PARAMETERS; }
    SWParameters getReadToHaplotypeSWParameters() { return DEFAULT_READ_TO_HAPLOTYPE_SMITH_WATERMAN_PARAMETERS; }

    // FeatureInput<VariantContext> alleles;

    uint8_t refModelDelQual{30};  // ReferenceConfidenceModel.REF_MODEL_DELETION_QUAL;  // 30
    int informativeReadOverlapMargin{2};
    int flowAssemblyCollapseHKerSize{0};
    int maxMnpDistance;
    bool softClipLowQualityEnds{false};
    bool flowAssemblyCollapsePartialMode{false};
    bool filterAlleles{false};
    bool filterLoneAlleles{false};
    bool writeFilteringGraphs{false};
    bool forceCallFiltered{false};
    float prefilterQualThreshold{PREFILTER_QUAL_THRESHOLD};
    float prefilterSorThreshold{PREFILTER_SOR_THRESHOLD};
    struct argument_range
    {
        int32_t min;
        int32_t max;
    };

public:
    RovacaToolArgs(int argc, char* argv[]);
    bool valid_check();
    static void valid_range(const std::string& name, int value, const argument_range& range);
    static void usage();
    const std::string& tool() const { return tool_name_; }
    const std::string& version() const { return main_version_; }
    const std::vector<std::string>& bam_path() const { return input_file_; }
    const std::string& vcf_path() const { return out_file_; }
    const std::string& reference_path() const { return reference_file_; }
    const std::string& bed_path() const { return bed_file_; }
    int32_t interval_padding() const { return interval_padding_; }
    int32_t max_reads_depth() const { return max_reads_depth_; }
    int32_t base_quality_score_threshold() const { return base_quality_score_threshold_; }
    int32_t run_pool_size() const { return run_pool_size_; }
    int32_t stream_pool_size() const { return stream_pool_size_; }
    bool create_output_index() const { return create_output_index_; }
    bool create_output_md5() const { return create_output_md5_; }
    const char* target_span() const { return target_span_.empty() ? nullptr : target_span_.data(); }
    const std::string& command_line() const { return command_line_; }
    const std::vector<int32_t>& gq_bands() const { return gq_bands_; }
    const char* assemble_output_path() const { return assemble_output_path_.empty() ? nullptr : assemble_output_path_.data(); }
    const std::string& pcr_model() const { return pcr_indel_model_; }
    const std::string& erc_model() const { return reference_confidence_mode_; }
    const std::string& recal_table() const { return recal_table_; }
    bool inspect_reads() const { return inspect_reads_; };
    bool old_pairhmm_engine() const { return old_pairhmm_engine_; }
    int32_t compression_level() const { return compression_level_; };
    const std::string& dbsnp_path() const { return dnsnp_file_; }

private:
    po::variables_map vm_;
    std::string tool_name_;
    std::vector<std::string> input_file_;
    std::string out_file_;
    std::string reference_file_;
    std::string bed_file_;
    std::string target_span_;
    std::string command_line_;
    std::string main_version_;
    std::string assemble_output_path_;
    std::string pcr_indel_model_;
    std::string reference_confidence_mode_;
    std::string recal_table_;
    int32_t interval_padding_;
    int32_t max_reads_depth_;
    int32_t base_quality_score_threshold_;
    int32_t run_pool_size_;
    int32_t stream_pool_size_;
    std::vector<int32_t> gq_bands_;
    bool create_output_index_;
    bool create_output_md5_;
    bool inspect_reads_;
    bool old_pairhmm_engine_;
    int32_t compression_level_;
    std::string dnsnp_file_;

    static constexpr const int32_t DEFAULT_MAX_READS_DEPTH = 50;
    static constexpr const int32_t DEFAULT_RUNPOOL_SIZE = 30;
    static constexpr const int32_t DEFAULT_INTERVAL_PADDING = 0;
    static constexpr const int32_t DEFAULT_BASE_QUALITY_SCORE_THRESHOLD = 18;
    static constexpr const int32_t DEFAULT_IOSTREAM_POOL_SIZE = 10;
    static constexpr const argument_range MAX_READS_DEPTH_RANGE = {0, INT32_MAX};
    static constexpr const argument_range RUNPOOL_SIZE_RANGE = {1, 128};
    static constexpr const argument_range INTERVAL_PADDING_RANGE = {0, INT32_MAX};
    static constexpr const argument_range BASE_QUALITY_SCORE_THRESHOLD_RANGE = {6, 127};
    static constexpr const argument_range IOSTREAM_POOL_SIZE_RANGE = {1, 20};
    static constexpr const argument_range COMPRESSION_LEVEL_RANGE = {0, 9};

    static constexpr const char* INPUT_PATH_NAME = "input,I";
    static constexpr const char* OUTUT_PATH_NAME = "output,O";
    static constexpr const char* REFERENCE_PATH_NAME = "reference,R";
    static constexpr const char* INTERVAL_PATH_NAME = "interval,L";
    static constexpr const char* INTERVAL_PADDING_NAME = "interval-padding,P";
    static constexpr const char* DOWNSAMPLE_THRESHOLD_NAME = "max-reads-depth,D";
    static constexpr const char* BASE_QUALITY_SCORE_THRESHOLD_NAME = "base-quality-score-threshold,Q";
    static constexpr const char* TARGET_SPAN_NAME = "target-span,T";
    static constexpr const char* GVCF_GQ_BANDS_NAME = "gvcf-gq-bands,G";
    static constexpr const char* GVCF_GQ_BANDS_LONGNAME = "gvcf-gq-bands";
    static constexpr const char* THREADPOOL_LONGNAME = "nthreads";
    static constexpr const char* IOSTREAM_POOL_SIZE = "nstreampool";
    static constexpr const char* ASSEMBLE_OUTPUT_PATH_LONGNAME = "assemble_output";
    static constexpr const char* HELP_NAME = "help,H";
    static constexpr const char* HELP_LONGNAME = "help";
    static constexpr const char* VERSION_NAME = "version,V";
    static constexpr const char* VERSION_LONGNAME = "version";
    static constexpr const char* PCR_INDEL_MODEL = "pcr-indel-model";
    static constexpr const char* EMIT_REF_CONFIDENCE = "emit-ref-confidence";
    static constexpr const char* CREATE_OUTPUT_INDEX = "index";
    static constexpr const char* CREATE_OUTPUT_MD5 = "md5";
    static constexpr const char* INSPECT_READS = "inspect-reads";
    static constexpr const char* PAIRHMM_ENGINE = "old-pairhmm-engine";
    static constexpr const char* COMPRESSION_LEVEL = "compression-level";
    static constexpr const char* BQSR_RECAL_TABLE = "bqsr-recal-table";
    static constexpr const char* DBSNP = "dbsnp";
};

// clang-format off
RovacaToolArgs::RovacaToolArgs(int argc, char* argv[])
{
    tool_name_ = argv[1];
    main_version_ = MAIN_VERSION;
    po::options_description haplotypecaller("HaplotypeCaller options");
    po::options_description bqsr("BQSR options");
    po::options_description all("All options");

    haplotypecaller.add_options()(HELP_NAME, "produce help message")(
        VERSION_NAME, "display version information")(
        INPUT_PATH_NAME, po::value<std::vector<std::string>>(&input_file_)->required()->multitoken(), "input file")(
        OUTUT_PATH_NAME, po::value<std::string>(&out_file_)->required(), "output file")(
        REFERENCE_PATH_NAME, po::value<std::string>(&reference_file_)->required(), "reference file")(
        INTERVAL_PATH_NAME, po::value<std::string>(&bed_file_), "interval file")(
        INTERVAL_PADDING_NAME, po::value<int32_t>(&interval_padding_)->default_value(DEFAULT_INTERVAL_PADDING)->notifier([](const int32_t& value){valid_range(INTERVAL_PADDING_NAME, value, INTERVAL_PADDING_RANGE);}), "padding around intervals")(
        DOWNSAMPLE_THRESHOLD_NAME, po::value<int32_t>(&max_reads_depth_)->default_value(DEFAULT_MAX_READS_DEPTH)->notifier([](const int32_t& value){valid_range(DOWNSAMPLE_THRESHOLD_NAME, value, MAX_READS_DEPTH_RANGE);}), "maximum reads per alignment start")(
        BASE_QUALITY_SCORE_THRESHOLD_NAME, po::value<int32_t>(&base_quality_score_threshold_)->default_value(DEFAULT_BASE_QUALITY_SCORE_THRESHOLD)->notifier([](const int32_t& value){valid_range(BASE_QUALITY_SCORE_THRESHOLD_NAME,value, BASE_QUALITY_SCORE_THRESHOLD_RANGE);}), "base qualities")(
        TARGET_SPAN_NAME, po::value<std::string>(&target_span_), "target span to call")(
        GVCF_GQ_BANDS_NAME,po::value<std::vector<int32_t>>()->multitoken(), "gvcf GQ bands")(
        THREADPOOL_LONGNAME,po::value<int32_t>(&run_pool_size_)->default_value(DEFAULT_RUNPOOL_SIZE)->notifier([](const int32_t& value){valid_range(THREADPOOL_LONGNAME,value, RUNPOOL_SIZE_RANGE);}), "n threads to launch")(
        ASSEMBLE_OUTPUT_PATH_LONGNAME, po::value<std::string>(&assemble_output_path_), "path to write assemble result")(
        PCR_INDEL_MODEL, po::value<std::string>(&pcr_indel_model_)->default_value("CONSERVATIVE"), "exclusive upper bounds for reference confidence gq bands")(
        EMIT_REF_CONFIDENCE, po::value<std::string>(&reference_confidence_mode_)->default_value("NONE"), "mode for emitting reference confidence scores")(
        IOSTREAM_POOL_SIZE, po::value<int32_t>(&stream_pool_size_)->default_value(DEFAULT_IOSTREAM_POOL_SIZE)->notifier([](const int32_t& value){valid_range(IOSTREAM_POOL_SIZE,value, IOSTREAM_POOL_SIZE_RANGE);}), "mode for emitting reference confidence scores")(
        CREATE_OUTPUT_INDEX, po::bool_switch(&create_output_index_)->default_value(true)->implicit_value(false), "create an index for output file")(
        CREATE_OUTPUT_MD5,po::bool_switch(&create_output_md5_)->default_value(false)->implicit_value(true), "create a md5 for output file")(
        INSPECT_READS,po::bool_switch(&inspect_reads_)->default_value(false)->implicit_value(true), "inspect reads")(
        PAIRHMM_ENGINE,po::bool_switch(&old_pairhmm_engine_)->default_value(false)->implicit_value(true), "old pairhmm engine(intel)")(
        COMPRESSION_LEVEL,po::value<int32_t>(&compression_level_)->default_value(6)->notifier([](int32_t value){valid_range(COMPRESSION_LEVEL, value, COMPRESSION_LEVEL_RANGE);}), "compression level")(
        DBSNP, po::value<std::string>(&dnsnp_file_), "dbSNP file");
    
    bqsr.add_options()(BQSR_RECAL_TABLE, po::value<std::string>(&recal_table_),"bqsr recal table file.");
    all.add(haplotypecaller).add(bqsr);
    // clang-format on
    try {
        po::store(po::parse_command_line(argc, argv, all), vm_);

        if (vm_.count(VERSION_LONGNAME)) {
            fprintf(stdout,
                    "version                        = %s\n"
                    "dev version                    = %s\n"
                    "dev branch                     = %s\n"
                    "dev id                         = %s\n\n",
                    MAIN_VERSION, GIT_DES, GIT_BRANCH, GIT_ID_LONG);
            exit(0);
        }

        if (vm_.count(HELP_LONGNAME)) {
            usage();
            exit(0);
        }

        po::notify(vm_);
    }
    catch (const po::error& e) {
        RovacaLogger::error("{}", e.what());
        usage();
        exit(EXIT_FAILURE);
    }

    if (vm_.count(GVCF_GQ_BANDS_LONGNAME)) {
        gq_bands_ = vm_[GVCF_GQ_BANDS_LONGNAME].as<std::vector<int32_t>>();
    }

    std::string path{argv[0]};
    auto pos = path.find_last_of('/');
    command_line_ = pos == std::string::npos ? path : path.substr(pos + 1);
    for (int32_t i = 1; i < argc; ++i) {
        command_line_.append(1, ' ').append(argv[i]);
    }
}

// clang-format off
void RovacaToolArgs::usage()
{
    std::cout << std::endl;
    std::cout << "Usage: rovaca <tool> [-option]" << std::endl;
    std::cout << "Required:" << std::endl;
    std::cout << "  -I, --input <file1> <file2> ...             input file paths" << std::endl;
    std::cout << "  -O, --output <file>                         output file path" << std::endl;
    std::cout << "  -R, --reference <file>                      reference file path" << std::endl;
    std::cout << "Options:" << std::endl;
    std::cout << "  -H, --help                                  display help message" << std::endl;
    std::cout << "  -V, --version                               display version message" << std::endl;
    std::cout << "  -L, --interval <file>                       interval file path" << std::endl;
    std::cout << "  -P, --interval-padding <int>                interval padding size, must be a non-negative integer (default: 0)" << std::endl;
    std::cout << "  -Q, --base-quality-score-threshold <int>    base qualities threshold, must be in [6, 127](default: 18)" << std::endl;
    // std::cout << "  -T, --target-span <str>                     specify target chromosome interval (e.g., chrM, chr1:1-10000)" << std::endl;
    std::cout << "  -D, --max-reads-depth <int>                 maximum reads depth per alignment start position (default: 50)" << std::endl;
    std::cout << "                                              must be a non-negative integer,set to 0 to disable" << std::endl;
    std::cout << "  -G, --gvcf-gq-bands <int1> <int2>...        specify GQ bands for merging non-variant sites in GVCF mode" << std::endl;
    std::cout << "                                              must be in [1, 100], (default: 1, 2, 3,... 60, 70, 80, 90, 99)" << std::endl;
    std::cout << "      --bqsr-recal-table <file>               specify the recalibration.table file for base quality scores recalibration(BQSR)." << std::endl;
    std::cout << "      --nthreads <int>                        number of threads to use, must be in [1, 128] (default: 30)" << std::endl;
    std::cout << "      --pcr-indel-model <str>                 PCR indel model (default: CONSERVATIVE)" << std::endl;
    std::cout << "                                              available options: {NONE, HOSTILE, CONSERVATIVE, AGGRESSIVE}" << std::endl;
    std::cout << "      --emit-ref-confidence <str>             emit reference confidence score mode (default: NONE)" << std::endl;
    std::cout << "                                              available options: {NONE, GVCF}" << std::endl;
    std::cout << "      --nstreampool <int>                     iostream pool size, must be in [1, 20] (default: 10)" << std::endl;
    // std::cout << "      --index                                 create index for output file (default: true), ture if specified" << std::endl;
    std::cout << "      --inspect-reads                         strictly inspect input reads (default: false)" << std::endl;
    std::cout << "      --old-pairhmm-engine                    use intel pairhmm engine (default: rovaca pairhmm engine)" << std::endl;
    std::cout << "      --compression-level                     compression level, must be in [0, 9] (default: 6)" << std::endl;
    // std::cout << "      --md5                                   create MD5 for output file (default: fasle), ture if specified" << std::endl;
    std::cout << "      --dbsnp <file>                          dbSNP file  Default value: null." << std::endl;
    std::cout << std::endl;
}
// clang-format on

bool RovacaToolArgs::valid_check()
{
    // required option has checked.
    for (const auto& read_path : input_file_) {
        CHECK_FILE_EXIST(read_path);
        CHECK_FILE_READALBE(read_path);
    }
    CHECK_FILE_EXIST(reference_file_);
    CHECK_FILE_READALBE(reference_file_);
    CHECK_FILE_WRITEABLE(out_file_.c_str());
    return true;
}

void RovacaToolArgs::valid_range(const std::string& name, int value, const argument_range& range)
{
    if (value < range.min || value > range.max) {
        RovacaLogger::error("the argument for option: {}, out of range", name);
        usage();
        exit(EXIT_FAILURE);
    }
}

#endif  // Rovaca_TOOL_ARGS_H