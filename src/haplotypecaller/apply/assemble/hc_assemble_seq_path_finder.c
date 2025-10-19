/**
 * @file hc_assemble_seq_path_finder.c
 * @author Cao Ce (caoce@genomics.cn)
 * @brief Seq Graph 生成Path及Cigar
 * @version 0.1
 * @date 2021-12-14
 *
 * @copyright Copyright (c) 2021 Rovaca SDK
 *
 */
#include "hc_assemble_seq_path_finder.h"

#include "hc_apply.h"
#include "hc_assemble.h"
#include "hc_assemble_dijkstra_shortest_path.h"
#include "hc_func.h"
#include "hc_marco.h"

#define MIN_HAPLOTYPE_REFERENCE_LENGTH (30)

static int hc_assemble_seq_finder_is_cigar_failure(p_hc_assemble_dijkstra_out_path_storage st, p_assemble_dijkstra_path all_st);

/**
 * Find discover paths by using KBestHaplotypeFinder over each graph.
 *
 * This method has the side effect that it will annotate all of the AssemblyResults objects with the derived haplotypes
 * which can be used for basing kmer graph pruning on the discovered haplotypes.
 *
 * @param assemblyResultByGraph assembly results objects keyed by graphs used to construct them
 * @param refHaplotype          reference haplotype
 * @param refLoc                location of reference haplotype
 * @param activeRegionWindow    window of the active region (without padding)
 * @param resultSet             (can be null) the results set into which to deposit discovered haplotypes
 * @param aligner               SmithWaterman aligner to use for aligning the discovered haplotype to the reference haplotype
 * @return A list of discovered haplotyes (note that this is not currently used for anything)
 */
void hc_assemble_seq_finder_find_best_paths(p_hc_apply apply)
{
    // add the reference haplotype separately from all the others to ensure that it is present in the list of haplotypes
    p_assemble_dijkstra_path all_st = &apply->read_graph.dijkstra_path_finder;
    p_assemble_graph_hash_node source_el, sink_el;
    p_assemble_vertex source, sink;
    p_hc_shared_lib_list_item item, tmp1, tmp2;
    p_assemble_dijkstra_one_path_node kBestHaplotype;

    // Validate that the graph is valid with extant source and sink before operating
    source_el = hc_assemble_utils_get_reference_source_vertex_with_head(
        apply->read_thread_graph->vertics_hash->node_list, apply->read_thread_graph->vertics_hash->node_list, apply->read_thread_graph);
    sink_el = hc_assemble_utils_get_reference_sink_vertex_with_head(
        apply->read_thread_graph->vertics_hash->node_list, apply->read_thread_graph->vertics_hash->node_list, apply->read_thread_graph);
    if (!source_el || !source_el->value) {
        return;
    }
    else {
        source = source_el->value;
    }
    if (!sink_el || !sink_el->value) {
        sink = NULL;
    }
    else {
        sink = sink_el->value;
    }

    hc_assemble_dijkstra_find_best_haplotypes(source, sink, apply);

    CDL_FOREACH_SAFE(all_st->result->nodes, item, tmp1, tmp2)
    {
        kBestHaplotype = item->data;
        if (kBestHaplotype->done) {
            continue;
        }

        kBestHaplotype->kmer_size = apply->read_graph.kmer;

        // TODO this score seems to be irrelevant at this point...
        if (kBestHaplotype->is_reference) {
            all_st->ref_haplotype = kBestHaplotype;
        }
        if (!kBestHaplotype->st) {
            CDL_DELETE(all_st->result->nodes, item);
            continue;
        }

        hc_assemle_cigar_cacl_calculate_cigar(apply, kBestHaplotype);

        // N cigar elements means that a bubble was too divergent from the reference so skip over this path
        if (!kBestHaplotype->st->cigar_len || hc_assemble_seq_finder_is_cigar_failure(kBestHaplotype->st, all_st)) {
            CDL_DELETE(all_st->result->nodes, item);
            continue;
        }

        kBestHaplotype->st->alignmentStartHapwrtRef = apply->read_thread_graph->origin_ref_start + apply->sw_avx->alignment_offset;
        kBestHaplotype->st->genomeLocation.from = apply->span_start;
        kBestHaplotype->st->genomeLocation.to = apply->span_end;
        kBestHaplotype->done = 1;
    }
}

/**
 * @return The number of reference bases that the read covers, excluding padding.
 */
static int hc_assemble_seq_finder_is_cigar_failure(p_hc_assemble_dijkstra_out_path_storage st, p_assemble_dijkstra_path all_st)
{
    uint32_t length = 0;
    uint32_t op, opr, opl, i = 0;

    for (; i < st->cigar_len; i++) {
        op = st->cigar[i];
        opr = bam_cigar_op(op);
        opl = bam_cigar_oplen(op);

        switch (opr) {
            case WORKER_SIGAR_STATUS_NSTATUS: return 1;
            case WORKER_SIGAR_STATUS_MATCH:
            case WORKER_SIGAR_STATUS_DELETION:
            case WORKER_SIGAR_S_MATCH:
            case WORKER_SIGAR_S_MISMATCH: length += opl; break;
            default: break;
        }
    }

    return length < MIN_HAPLOTYPE_REFERENCE_LENGTH || length != (all_st->ref_len - strlen(HC_ASSEMBLE_SW_PAD_STRING) * 2);
}

/**
 * Release Grapgh内Cycle的检测修复
 *
 * 实际极少数据会走这里，有特例再进行实现
 */
#if 0
/**
 * Removes edges that produces cycles and also dead vertices that do not lead to any sink vertex.
 * @return never {@code null}.
 */
private BaseGraph<V, E> removeCyclesAndVerticesThatDontLeadToSinks(final BaseGraph<V, E> original, final Collection<V> sources, final Set<V> sinks) {
    final Set<E> edgesToRemove = new HashSet<>(original.edgeSet().size());
    final Set<V> vertexToRemove = new HashSet<>(original.vertexSet().size());

    boolean foundSomePath = false;
    for (final V source : sources) {
        final Set<V> parentVertices = new HashSet<>(original.vertexSet().size());
        foundSomePath = findGuiltyVerticesAndEdgesToRemoveCycles(original, source, sinks, edgesToRemove, vertexToRemove, parentVertices) || foundSomePath;
    }

    Utils.validate(foundSomePath, () -> "could not find any path from the source vertex to the sink vertex after removing cycles: "
            + Arrays.toString(sources.toArray()) + " => " + Arrays.toString(sinks.toArray()));

    Utils.validate(!(edgesToRemove.isEmpty() && vertexToRemove.isEmpty()), "cannot find a way to remove the cycles");

    final BaseGraph<V, E> result = original.clone();
    result.removeAllEdges(edgesToRemove);
    result.removeAllVertices(vertexToRemove);
    return result;
}

/**
 * Recursive call that looks for edges and vertices that need to be removed to get rid of cycles.
 *
 * @param graph the original graph.
 * @param currentVertex current search vertex.
 * @param sinks considered sink vertices.
 * @param edgesToRemove collection  of edges that need to be removed in order to get rid of cycles.
 * @param verticesToRemove collection of vertices that can be removed.
 * @param parentVertices collection of vertices that preceded the {@code currentVertex}; i.e. the it can be
 *                       reached from those vertices using edges existing in {@code graph}.
 *
 * @return {@code true} to indicate that the some sink vertex is reachable by {@code currentVertex},
 *  {@code false} otherwise.
 */
private boolean findGuiltyVerticesAndEdgesToRemoveCycles(p_hc_apply hc_apply,
                                                            p_assemble_vertex currentVertex,
                                                            final Set<V> sinks,
                                                            final Set<E> edgesToRemove,
                                                            final Set<V> verticesToRemove,
                                                            final Set<V> parentVertices) 
{
    if (sinks.contains(currentVertex)) {
        return true;
    }

    final Set<E> outgoingEdges = graph.outgoingEdgesOf(currentVertex);
    parentVertices.add(currentVertex);

    boolean reachesSink = false;
    for (final E edge : outgoingEdges) {
        final V child = graph.getEdgeTarget(edge);
        if (parentVertices.contains(child)) {
            edgesToRemove.add(edge);
        } else {
            final boolean childReachSink = findGuiltyVerticesAndEdgesToRemoveCycles(graph, child, sinks, edgesToRemove, verticesToRemove, parentVertices);
            reachesSink = reachesSink || childReachSink;
        }
    }
    if (!reachesSink) {
        verticesToRemove.add(currentVertex);
    }
    return reachesSink;
}
#endif