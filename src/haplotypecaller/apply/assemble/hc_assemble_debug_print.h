#ifndef HC_ASSEMBLE_DEBUG_PRINT_H
#define HC_ASSEMBLE_DEBUG_PRINT_H
#include "hc_apply.h"

void hc_assemble_utils_print_reads_after_finalize_region(p_hc_apply apply);
void hc_assemble_utils_print_path(p_hc_apply apply);
void hc_assemble_graph_layout(p_assemble_graph graph, const char* path);
void hc_assemble_graph_layout_vertex(p_assemble_graph graph, const char* path);

#endif  // HC_ASSEMBLE_DEBUG_PRINT_H