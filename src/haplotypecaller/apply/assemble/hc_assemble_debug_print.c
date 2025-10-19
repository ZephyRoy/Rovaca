#include "hc_assemble_debug_print.h"

#include "assemble_interface.h"
#include "rbtree.h"
/**
 * @brief
 *
 * @param apply
 */
void hc_assemble_utils_print_reads_after_finalize_region(p_hc_apply apply)
{
    p_hc_apply_one_read read;
    struct rb_node* node;
    printf("\n=================================================\n");
    printf("Region after finalized:\t%d-%d\treads count: %d\n", apply->span_start, apply->span_end, apply->reads_count);

    for (node = rb_first(&apply->read_sorted); node; node = rb_next(node)) {
        read = rb_entry(node, hc_apply_one_read, node);
        printf("%s\t%ld-%ld\t%s\n", read->read_data, read->pos_start, read->pos_end, read->cigar_str);
    }
}

void hc_assemble_utils_print_path(p_hc_apply apply)
{
    p_assemble_dijkstra_path all_st = &apply->read_graph.dijkstra_path_finder;
    p_assemble_dijkstra_one_path_node kBestHaplotype;
    p_hc_shared_lib_list_item item;
    uint32_t hsum = 0, i;
    int cached_length;
    const char* cigar_string = "MIDNSHP=XBsf";

    puts("\n");
    printf(
        "=================================================\nassemblyRegion: "
        "%d:%" PRId64 "-%" PRId64 "\n=================================================\n\n",
        apply->region->tid, apply->region->start_index, apply->region->end_index);
    CDL_COUNT(all_st->result->nodes, item, hsum);
    printf("Unclipped Haplotypes(%u):\n", hsum);

    CDL_FOREACH(all_st->result->nodes, item)
    {
        kBestHaplotype = item->data;

        printf("[%u-%u] k=%u len: %lu ", apply->span_start, apply->span_end, kBestHaplotype->is_reference ? 0 : kBestHaplotype->kmer_size,
               kBestHaplotype->st->seq_len - strlen(HC_ASSEMBLE_SW_PAD_STRING) * 2);
        for (i = 0; i < kBestHaplotype->st->cigar_len; i++) {
            printf("%u%c", bam_cigar_oplen(kBestHaplotype->st->cigar[i]), cigar_string[bam_cigar_op(kBestHaplotype->st->cigar[i])]);
        }
        if (kBestHaplotype->is_reference) {
            puts("ref");
        }
        else {
            puts("");
        }
        cached_length = kBestHaplotype->st->seq_len - strlen(HC_ASSEMBLE_SW_PAD_STRING) * 2;
        printf("%.*s\n", cached_length, kBestHaplotype->st->seq + strlen(HC_ASSEMBLE_SW_PAD_STRING));
    }
}

void hc_assemble_graph_layout(p_assemble_graph graph, const char* path)
{
    FILE* fp = fopen(path, "w");
    p_assemble_edge edge;
    p_assemble_graph_hash_node vertex_st;
    p_assemble_vertex vertex;
    unsigned int weight;

    if (!fp) {
        return;
    }

    // head
    fwrite("digraph assemblyGraphs {\n", strlen("digraph assemblyGraphs {\n"), sizeof(char), fp);
    // edge
    CDL_FOREACH2(graph->edge_list, edge, link.all_st.next)
    {
        weight = edge->weight;
        fprintf(fp, "vertex%p -> vertex%p [label=\"%u\", color=%s];\n", edge->link.from.vertex, edge->link.to.vertex, weight,
                edge->is_ref ? "red" : "black");
    }
    // vertex
    CDL_FOREACH(graph->vertics_hash->node_list, vertex_st)
    {
        vertex = vertex_st->value;
        fprintf(fp, "vertex%p [label=\"%.*s\", shape=box]\n", vertex, vertex->data_len, (char*)vertex->data);
    }
    // end
    fputs("}", fp);
    fclose(fp);
}

void hc_assemble_graph_layout_vertex(p_assemble_graph graph, const char* path)
{
    FILE* fp = fopen(path, "w");
    p_assemble_edge edge;
    p_assemble_graph_hash_node vertex_st;
    p_assemble_vertex vertex;

    if (!fp) {
        return;
    }

    // vertex
    CDL_FOREACH(graph->vertics_hash->node_list, vertex_st)
    {
        vertex = vertex_st->value;
        fprintf(fp, "vertex%p [label=\"%.*s\", shape=box] ", vertex, vertex->data_len, (char*)vertex->data);
        CDL_FOREACH2(vertex->from_edges, edge, link.to.next)
        {
            if (edge->link.from.vertex) {
                fprintf(fp, " from vertex %p [label=\"%.*s\", shape=box] ", edge->link.from.vertex, edge->link.from.vertex->data_len,
                        (char*)edge->link.from.vertex->data);
            }
        }
        CDL_FOREACH2(vertex->out_edges, edge, link.from.next)
        {
            if (edge->link.to.vertex) {
                fprintf(fp, " to vertex %p [label=\"%.*s\", shape=box] ", edge->link.to.vertex, edge->link.to.vertex->data_len,
                        (char*)edge->link.to.vertex->data);
            }
        }
        fprintf(fp, " \n");
    }
    fclose(fp);
}