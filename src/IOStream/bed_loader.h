#ifndef BED_READER_H
#define BED_READER_H

#include <map>
#include <stdexcept>
#include <string>
#include <vector>

#include "fasta_loader.h"
#include "header.h"
#include "htslib/hts.h"

static const int32_t WES_PREPADDING = 100;
class BedLoader
{
public:
    explicit BedLoader(const std::string& path, const int padding, const contig_info_t& fasta_dict)
        : bed_lists({})
        , bed_key({})
        , all_intervals(nullptr)
        , all_padded_intervals_(nullptr)
        , reference_dict_(fasta_dict)
    {
        if (fasta_dict.dict.empty() || fasta_dict.key.empty()) {
            throw(std::runtime_error("Bed must be initialized after fasta."));
        }
        load_all_bed_intervals(path, padding, fasta_dict);
    }
    ~BedLoader() { release_bed_intervals(); }

    const std::map<std::string, p_bed_intervals>& get_bed_intervals() { return bed_lists; }
    const std::vector<std::string>& get_bed_keys() { return bed_key; }
    p_bed_intervals get_all_intervals() const { return all_intervals; }
    p_bed_intervals get_all_padded_intervals() const { return all_padded_intervals_; }
    static p_bed_intervals init_bed_intervals();

private:
    std::map<std::string, p_bed_intervals> bed_lists;
    std::vector<std::string> bed_key;
    p_bed_intervals all_intervals;
    p_bed_intervals all_padded_intervals_;
    const contig_info_t& reference_dict_;

private:
    void load_all_bed_intervals(const std::string& path, const int padding, const contig_info_t& fasta_dict);
    void release_bed_intervals();
};

#endif  // BED_READER_H