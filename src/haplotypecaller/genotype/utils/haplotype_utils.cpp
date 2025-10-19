#include "haplotype_utils.h"

#include "genotype_macors.h"
#include "haplotype.h"
#include "simple_interval.h"
#include "utils/alignment_utils.h"
#include "utils/cigar_utils.h"

namespace rovaca
{

pHaplotype HaplotypeUtils::trim(pHaplotype h, pSimpleInterval interval, bool ignore_ref_state, pMemoryPool pool)
{
    CHECK_CONDITION_EXIT(!h->interval()->contains(*interval), "can only trim a Haplotype to a containing span");

    int64_t new_start = interval->get_start() - h->interval()->get_start();
    int64_t new_stop = new_start + interval->get_stop() - interval->get_start();

    pBases old_bases = h->get_bases();
    pCigar old_cigar = h->cigar();
    pBases new_bases = AlignmentUtils::get_bases_covering_ref_interval(new_start, new_stop, old_bases, 0, old_cigar, pool);
    if (nullptr == new_bases || 0 == new_bases->num) {
        return nullptr;
    }

    pCigar new_cigar = AlignmentUtils::trim_cigar_by_reference(old_cigar, new_start, new_stop, pool);
    bool leading_insertion = !consumes_ref_bases(bam_cigar_op(new_cigar->data[0]));
    bool trailing_insertion = !consumes_ref_bases(bam_cigar_op(new_cigar->data[new_cigar->num - 1]));
    uint32_t first_index_to_keep_inclusive = leading_insertion ? 1 : 0;
    uint32_t last_index_to_keep_exclusive = new_cigar->num - (trailing_insertion ? 1 : 0);
    if (last_index_to_keep_exclusive <= first_index_to_keep_inclusive) {
        // edge case of entire cigar is insertion
        return nullptr;
    }

    pCigar leading_indel_trimmed_new_cigar;
    if (!(leading_insertion || trailing_insertion)) {
        leading_indel_trimmed_new_cigar = new_cigar;
    }
    else {
        leading_indel_trimmed_new_cigar =
            CigarBuilder::create(false, pool)
                ->add_all(new_cigar->data + first_index_to_keep_inclusive, last_index_to_keep_exclusive - first_index_to_keep_inclusive)
                ->make();
    }

    pHaplotype new_h = Haplotype::create(pool);
    new_h->init_haplotype(new_bases, !ignore_ref_state && h->is_reference());
    new_h->set_cigar(leading_indel_trimmed_new_cigar);
    new_h->set_interval(interval);
    new_h->set_score(h->score());
    new_h->set_kmer_size(h->kmer_size());
    new_h->set_alignment_start_hap_wrt_ref(new_start + h->alignment_start_hap_wrt_ref());
    return new_h;
}

pHaplotype HaplotypeUtils::trim(pHaplotype h, pSimpleInterval interval, pMemoryPool pool) { return trim(h, interval, false, pool); }

}  // namespace rovaca