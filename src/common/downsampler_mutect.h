#ifndef DOWNSAMPLER_MUTECT_H
#define DOWNSAMPLER_MUTECT_H

#include "downsampler.h"

class MutectDownsampler : public Downsampler
{
public:
    MutectDownsampler();
    ~MutectDownsampler();
};

MutectDownsampler::MutectDownsampler() {}

MutectDownsampler::~MutectDownsampler() {}

#endif  // DOWNSAMPLER_MUTECT_H