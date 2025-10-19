#ifndef READS_FILTER_HC_H
#define READS_FILTER_HC_H

#include <iostream>

#include "htslib/sam.h"
#include "../rovaca_logger/rovaca_logger.h"
#include "reads_filter.h"
#include "reads_filter_lib.h"

#define HC_READ_FILTER_MINMQ (20)
#define HC_READ_FILTER_MAXMQ (255)

/************************************************
   HC中，对read主要做如下方面的过滤：
   1. mapping quality. (min-255)
   2. not unmapped.
   3. not secondary alignment.
   4. not duplicate.
   5. read fails platform/vendor quality checks.
   6. non zero reference length.
   7. good cigar operation().
   8. well formed.
*************************************************
*/
class HCReadFilter : public ReadFilter
{
public:
    HCReadFilter(bam_hdr_t* hdr, bool inspect);
    HCReadFilter(bam_hdr_t* hdr, int flags, int map_qual);
    ~HCReadFilter() { RovacaLogger::info("ReadsFilter done, {} reads were filtered", filterd_reads_count); }
    bool test(bam1_t* read) override;

private:
    bool read_mapq_check(bam1_t* read);

private:
    bool inspect_{false};
    uint32_t filterd_reads_count{0};
};

HCReadFilter::HCReadFilter(bam_hdr_t* hdr, bool inspect)
    : ReadFilter(hdr)
    , inspect_(inspect)
{
    flag = (BAM_FUNMAP | BAM_FSECONDARY | BAM_FDUP | BAM_FQCFAIL);
    min_mqual = HC_READ_FILTER_MINMQ;
    max_mqual = HC_READ_FILTER_MAXMQ;
}

HCReadFilter::HCReadFilter(bam_hdr_t* hdr, int flags, int map_qual)
    : ReadFilter(hdr)
{
    flag = flags;
    min_mqual = map_qual;
    max_mqual = HC_READ_FILTER_MAXMQ;
}

bool HCReadFilter::read_mapq_check(bam1_t* read)
{
    bool ret = (read->core.qual >= min_mqual) && (read->core.qual < max_mqual);
    return ret;
}

bool HCReadFilter::test(bam1_t* read)
{
    if (!read) {
        return false;
    }
    bool ret = !(read->core.flag & flag) && read_mapq_check(read) &&
               (!inspect_ || (ReadFilterLib::is_well_formed(read, bam_hdr) && ReadFilterLib::is_good_cigar(read)));
    if (!ret) filterd_reads_count++;
    return ret;
}

#endif  // READS_FILTER_HC_H