#ifndef TRANSFORMER_H
#define TRANSFORMER_H
#include <htslib/sam.h>

class Transformer
{
public:
    Transformer() {}
    virtual ~Transformer() {}
    virtual void apply(bam1_t* read) = 0;
    virtual void load_report() = 0;

private:
};

#endif  // TRANSFORMER_H