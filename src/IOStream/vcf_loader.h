#ifndef VCF_LOADER_H
#define VCF_LOADER_H
#include <mutex>
#include <set>
#include <string>
#include <vector>

#include "header.h"
#include "htslib/tbx.h"
#include "htslib/vcf.h"

class VcfLoader
{
private:
    /// @brief When a chromosome is finished, return all its bcf1_t to avoid repeated malloc
    class VcfPool
    {
    private:
        std::vector<bcf1_t*> data_{};
        std::mutex push_mutex_{};

    public:
        VcfPool() = default;
        ~VcfPool() = default;

    /// @brief Resource initialization
        void initialization();

    /// @brief Resource release
        void terminalization();

        bool push(bcf1_t* b);
        bool pop(bcf1_t** b);
    };

private:
    htsFile* fd_;
    bcf_hdr_t* hdr_;
    tbx_t* tbx_;
    VcfPool pool_;
    kstring_t cache_line_;
    const char** seqnames_;
    int32_t n_seqname_;
    std::set<std::string> seqnames_set_;

public:
    /// @brief Create a VcfLoader
    /// @param vcf_db dbsnp file
    VcfLoader()
        : fd_(nullptr)
        , hdr_(nullptr)
        , tbx_(nullptr)
        , pool_()
        , cache_line_{0, 0, nullptr}
        , seqnames_{nullptr}
        , n_seqname_{0}
        , seqnames_set_{}
    {}

    ~VcfLoader() = default;

    /// @brief Resource initialization
    /// @param vcf_db vcf file
    /// @param tp public thread pool, used for reading vcf
    /// @param info contig info from fasta, used to check for match
    void initialization(const std::string& vcf_db, htsThreadPool* tp, contig_info_t* info);

    /// @brief Resource release
    void terminalization();

    /// @brief Read vcf according to given parameters, two modes: if bed is nullptr, read the whole chromosome by chr; if bed is not nullptr, read all intervals corresponding to tid according to bed
    /// @param tid
    /// @param chr
    /// @param bed
    /// @return
    std::vector<bcf1_t*> load_vcf(const char* chr, p_bed_intervals bed);

    /// @brief When a chromosome is finished, db_manager returns the corresponding data
    /// @param b
    /// @return
    bool recover_bcf(bcf1_t* b) { return pool_.push(b); }
};

#endif  // VCF_LOADER_H