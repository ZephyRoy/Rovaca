#ifndef HC_ASSEMBLE_CIGAR_UTILS_H
#define HC_ASSEMBLE_CIGAR_UTILS_H

#include "assemble_interface.h"

hts_pos_t hc_assemble_utils_read_get_end(p_hc_apply_one_read read);
void merge_consecutive_identical_cigar(p_hc_apply_one_read read);

#endif  // HC_ASSEMBLE_CIGAR_UTILS_H