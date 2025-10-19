#include "fasta_loader.h"

#include <unistd.h>

#include <fstream>
#include <iostream>
#include <ostream>

FastaLoader::FastaLoader(const std::string& path)
    : fai(nullptr)
{
    fai = fai_load(path.c_str());
    if (!fai) {
        // statics_print_error("failed to open fasta file %s!", path.c_str());
        std::cerr << "" << std::endl;
        exit(-1);
    }
}

FastaLoader::~FastaLoader()
{
    if (fai) fai_destroy(fai);
}

std::vector<char*> FastaLoader::fetch_target_seq(const char* chr_name, const hts_pos_t* start, const hts_pos_t* end, int n)
{
    int len = 0;
    std::vector<char*> ret;
    for (int i = 0; i < n; ++i) {
        char* ref_bases = faidx_fetch_seq(fai, chr_name, *(start + i), *(end + i), &len);
        if (!ref_bases) {
            // statics_print_log("fetch ref seq failed. %s:%d-%d", chr_name, p_beg_i, p_end_i);
        }
        for (int i = 0; i < len; i++) {
            ref_bases[i] = seq_nt16_str_refine[seq_nt16_table[ref_bases[i] & 0xFF]];
        }
        ret.emplace_back(ref_bases);
    }
    return ret;
}

char* FastaLoader::fetch_target_seq(const char* chr_name, int p_beg_i, int p_end_i)
{
    int len = 0;
    char* ref_bases = faidx_fetch_seq(fai, chr_name, p_beg_i, p_end_i, &len);
    if (!ref_bases) {
        // statics_print_log("fetch ref seq failed. %s:%d-%d", chr_name, p_beg_i, p_end_i);
    }
    for (int i = 0; i < len; i++) {
        ref_bases[i] = seq_nt16_str_refine[seq_nt16_table[ref_bases[i] & 0xFF]];
    }
    return ref_bases;
}

char* FastaLoader::fetch_target_seq(int tid, int p_beg_i, int p_end_i)
{
    const char* chrom_name = faidx_iseq(fai, tid);
    return fetch_target_seq(chrom_name, p_beg_i, p_end_i);
}

char* FastaLoader::fetch_chr_seq(int tid)
{
    const char* chrom_name = faidx_iseq(fai, tid);
    return fetch_chr_seq(chrom_name);
}

char* FastaLoader::fetch_chr_seq(const char* chr_name)
{
    int len = 0;
    int target_len = faidx_seq_len(fai, chr_name);
    char* ref_bases = faidx_fetch_seq(fai, chr_name, 0, target_len, &len);
    if (!ref_bases) {
        // statics_print_log("fetch ref seq failed. %s", chr_name);
    }
    for (int i = 0; i < len; i++) {
        ref_bases[i] = seq_nt16_str_refine[seq_nt16_table[ref_bases[i] & 0xFF]];
    }
    return ref_bases;
}

void FastaLoader::get_fasta_dict(const std::string& fasa_path, contig_info_t* contig)
{
    std::string fai_path = fasa_path + ".fai";

    std::ifstream ifs(fai_path);
    std::string line;
    char contig_cache[64];
    uint64_t chr_len;

    int i = 0;

    while (std::getline(ifs, line)) {
        sscanf(line.c_str(), "%s\t%lu", contig_cache, &chr_len);
        contig->dict.insert({std::string(contig_cache), chr_len});
        contig->key.emplace_back(std::string(contig_cache));
        contig->idict[i] = chr_len;
        i++;
    }
    ifs.close();
}
