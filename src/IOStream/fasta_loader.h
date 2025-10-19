#ifndef FASTA_LOADER_H
#define FASTA_LOADER_H

#include <map>
#include <string>
#include <vector>

#include "header.h"
#include "htslib/faidx.h"
#include "htslib/sam.h"

const char seq_nt16_str_refine[] = "NACNGNNNTNNNNNNN";

class FastaLoader
{
private:
    faidx_t* fai;

public:
    explicit FastaLoader(const std::string& path);
    ~FastaLoader();

    std::vector<char*> fetch_target_seq(const char* chr_name, const hts_pos_t* start, const hts_pos_t* end, int n);

    char* fetch_target_seq(const char* chr_name, int p_beg_i, int p_end_i);

    char* fetch_target_seq(int tid, int p_beg_i, int p_end_i);

    char* fetch_chr_seq(const char* chr_name);

    char* fetch_chr_seq(int tid);

    static void get_fasta_dict(const std::string& fasa_path, contig_info_t* contig);
};
#endif  // FASTA_LOADER_H
