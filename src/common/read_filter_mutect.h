#ifndef READ_FILTER_MUTECT_H
#define READ_FILTER_MUTECT_H

#include <iostream>

#include "htslib/sam.h"
#include "reads_filter.h"

class MutectReadFilter : public ReadFilter
{
public:
    MutectReadFilter(bam_hdr_t* hdr);
    MutectReadFilter(bam_hdr_t* hdr, int flags, int map_qual);
    ~MutectReadFilter() {}
    bool test(bam1_t* read) override;
};

MutectReadFilter::MutectReadFilter(bam_hdr_t* hdr)
    : ReadFilter(hdr)
{
    flag = (BAM_FUNMAP | BAM_FSECONDARY | BAM_FDUP | BAM_FQCFAIL);
    min_mqual = 20;
}

MutectReadFilter::MutectReadFilter(bam_hdr_t* hdr, int flags, int map_qual)
    : ReadFilter(hdr)
{
    flag = flags;
    min_mqual = map_qual;
}

bool MutectReadFilter::test(bam1_t* read)
{
    if (!read) {
        return false;
    }
    return !(read->core.flag & flag) && (read->core.qual >= min_mqual);
}

#endif  // READ_FILTER_MUTECT_H