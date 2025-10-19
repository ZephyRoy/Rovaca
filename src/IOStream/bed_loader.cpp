#include "bed_loader.h"

#include <cmath>
#include <fstream>
#include <iostream>
#include <sstream>

void BedLoader::load_all_bed_intervals(const std::string &path, const int padding, const contig_info_t& fasta_dict)
{
    std::ifstream ifs(path);
    if (!ifs.is_open()) {
        throw(std::runtime_error("Init bed file error."));
    }
    std::string line;
    char contig_cache[64];
    hts_pos_t bed_start = 0;
    hts_pos_t bed_end = 0;
    all_intervals = init_bed_intervals();
    all_padded_intervals_ = init_bed_intervals();

    while (std::getline(ifs, line)) {
        if (line.front() != 'c') {
            continue;
        }
        sscanf(line.c_str(), "%s\t%ld\t%ld", contig_cache, &bed_start, &bed_end);
        std::string contig = std::string(contig_cache);

        if (!fasta_dict.dict.count(contig)) {
            throw(std::runtime_error("contig exists in bed, but cannot be found in fa"));
        }

        // User padding
        bed_start = std::max(bed_start - padding, (int64_t)0);
        bed_end = std::min(bed_end + padding, (int64_t)reference_dict_.dict.at(contig));

        if (!bed_lists.count(contig)) {
            p_bed_intervals one_chr_interval = init_bed_intervals();
            // first insert.
            one_chr_interval->n++;
            one_chr_interval->start[0] = bed_start;
            one_chr_interval->end[0] = bed_end;
            // one_chr_interval->regarry[0] = region;
            bed_lists.insert({contig, one_chr_interval});
            bed_key.emplace_back(contig);
        }
        else {
            p_bed_intervals bi = bed_lists[contig];
            // combine overlap.
            if (bi->n > 0 && bi->end[bi->n - 1] >= bed_end && bi->start[bi->n - 1] <= bed_start) {
                continue;
            }
            else if (bi->n > 0 && bi->end[bi->n - 1] > bed_start) {
                bi->end[bi->n - 1] = bed_end;
                continue;
            }
            else {
                bi->n++;
            }
            if (bi->n - 1 == bi->m) {
                bi->m = bi->m << 1;
                bi->start = (hts_pos_t *)realloc(bi->start, sizeof(hts_pos_t) * bi->m);
                bi->end = (hts_pos_t *)realloc(bi->end, sizeof(hts_pos_t) * bi->m);
                bi->regarry = (char **)realloc(bi->regarry, sizeof(char *) * bi->m);
            }
            bi->start[bi->n - 1] = bed_start;
            bi->end[bi->n - 1] = bed_end;
        }
    }

    for (auto &chrom : bed_key) {
        p_bed_intervals bi = bed_lists[chrom];
        if (bi->n + all_intervals->n > all_intervals->m) {
            all_intervals->m = bi->n + all_intervals->n;
            all_intervals->m = kroundup64(all_intervals->m);
            all_intervals->regarry = (char **)realloc(all_intervals->regarry, sizeof(char *) * all_intervals->m);
            if (all_intervals->regarry == nullptr) {
                throw(std::runtime_error("realloc error"));
            }
        }

        if (bi->n + all_padded_intervals_->n > all_padded_intervals_->m) {
            all_padded_intervals_->m = bi->n + all_padded_intervals_->n;
            all_padded_intervals_->m = kroundup64(all_padded_intervals_->m);
            all_padded_intervals_->regarry = (char **)realloc(all_padded_intervals_->regarry, sizeof(char *) * all_padded_intervals_->m);
            if (all_padded_intervals_->regarry == nullptr) {
                throw(std::runtime_error("realloc error"));
            }
        }

        for (int i = 0; i < bi->n; i++) {
            bed_start = bi->start[i];
            bed_end = bi->end[i];
            int reglength = strlen(chrom.c_str()) + trunc(log(bed_start + 1)) + 2 + trunc(log(bed_end)) + 3;
            char *reg = new char[reglength];
            sprintf(reg,
                    "%s:"
                    "%" PRId64
                    "-"
                    "%" PRId64 "",
                    chrom.c_str(), bed_start, bed_end);
            bi->regarry[i] = reg;
            all_intervals->regarry[all_intervals->n] = reg;
            all_intervals->n++;

            // Pre-padding with 100bps around.
            char *rg = new char[reglength + 8];
            sprintf(rg,
                    "%s:"
                    "%" PRId64
                    "-"
                    "%" PRId64 "",
                    chrom.c_str(), std::max(bed_start - WES_PREPADDING, int64_t(1)),
                    std::min(bed_end + WES_PREPADDING, int64_t(reference_dict_.dict.at(chrom))));
            all_padded_intervals_->regarry[all_padded_intervals_->n] = rg;
            all_padded_intervals_->n++;
        }
    }
    ifs.close();
}

p_bed_intervals BedLoader::init_bed_intervals()
{
    p_bed_intervals intervals = (p_bed_intervals)calloc(1, sizeof(bed_intervals));
    intervals->n = 0;
    intervals->m = 1;
    intervals->start = (hts_pos_t *)calloc(intervals->m, sizeof(hts_pos_t));
    intervals->end = (hts_pos_t *)calloc(intervals->m, sizeof(hts_pos_t));
    intervals->regarry = (char **)calloc(intervals->m, sizeof(char *));
    return intervals;
}

void BedLoader::release_bed_intervals()
{
    for (auto &chr_intervals : bed_lists) {
        p_bed_intervals intervals = chr_intervals.second;
        if (intervals == nullptr) return;
        if (intervals->start != nullptr) free(intervals->start);
        if (intervals->end != nullptr) free(intervals->end);
        if (intervals->regarry != nullptr) free(intervals->regarry);
        free(intervals);
    }
    if (all_intervals != nullptr) {
        if (all_intervals->start != nullptr) free(all_intervals->start);
        if (all_intervals->end != nullptr) free(all_intervals->end);
        for (int i = 0; i < all_intervals->n; ++i) {
            delete[] all_intervals->regarry[i];
        }
        if (all_intervals->regarry != nullptr) free(all_intervals->regarry);
        free(all_intervals);
    }

    if (all_padded_intervals_ != nullptr) {
        if (all_padded_intervals_->start != nullptr) free(all_padded_intervals_->start);
        if (all_padded_intervals_->end != nullptr) free(all_padded_intervals_->end);
        for (int i = 0; i < all_padded_intervals_->n; ++i) {
            delete[] all_padded_intervals_->regarry[i];
        }
        if (all_padded_intervals_->regarry != nullptr) free(all_padded_intervals_->regarry);
        free(all_padded_intervals_);
    }
}