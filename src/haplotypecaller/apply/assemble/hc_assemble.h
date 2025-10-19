#ifndef __HC_ASSEMBLE_H__
#define __HC_ASSEMBLE_H__

#include <stdint.h>
// #include "hc_main.h"
#include "utlist.h"

#define HC_ASSEMBLE_MIN_TAIL_QUALITY_TO_USE        (9)
#define HC_ASSEMBLE_GRAPH_MIN_QUAL                 (10)
#define HC_ASSEMBLE_GRAPH_HASH_MUL                 (31)
#define HC_ASSEMBLE_KMER_10                        (10)
#define HC_ASSEMBLE_KMER_25                        (25)
#define HC_ASSEMBLE_KMER_ITERATION_INCREASE        (10)
#define HC_ASSEMBLE_MAX_KMER_ITERATIONS_TO_ATTEMPT (6)
#define HC_ASSEMBLE_MAX_KMER                       (HC_ASSEMBLE_KMER_25 + HC_ASSEMBLE_KMER_ITERATION_INCREASE * HC_ASSEMBLE_MAX_KMER_ITERATIONS_TO_ATTEMPT)
#define HC_ASSEMBLE_GRAPH_CHAIN_PRUNE_FACTOR       (2)

#endif  // !__HC_ASSEMBLE_H__