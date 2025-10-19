#ifndef HEADER_H
#define HEADER_H
#include <map>
#include <string>
#include <vector>

#include "htslib/hts.h"

typedef struct bed_interval_t
{
    char** regarry;  // for htslib
    hts_pos_t* start;
    hts_pos_t* end;
    int m;  // capability of number
    int n;  // elements of number
} bed_intervals, *p_bed_intervals;

typedef struct contig_info_t
{
    std::map<int, uint64_t> idict;
    std::map<std::string, uint64_t> dict;
    std::vector<std::string> key;
}* p_contig_info_t;

#endif  // HEADER_H