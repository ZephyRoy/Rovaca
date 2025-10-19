#ifndef READS_FILTER_H
#define READS_FILTER_H

#include <iostream>

// #include "../rovaca_logger/rovaca_logger.h"
#include "htslib/sam.h"

class ReadFilter
{
public:
    bam_hdr_t* bam_hdr;

protected:
    int flag;
    int min_mqual;
    int max_mqual;

public:
    ReadFilter(bam_hdr_t* hdr)
        : bam_hdr(hdr)
    {
        if (hdr == nullptr) {
            // RovacaLogger::error("nullptr header catched in read filter");
            exit(EXIT_FAILURE);
        }
    }
    virtual ~ReadFilter() {}
    virtual bool test(bam1_t* read) = 0;
};

#endif  // READS_FILTER_H