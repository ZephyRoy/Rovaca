#include "read_record_utils.h"

#include <algorithm>

#include "bam_data_pool.hpp"
#include "cigar_utils.h"
#include "genotype_struct.h"
#include "htslib/sam.h"
#include "rovaca_logger.h"
#include "read_record.h"
#include "simple_interval.h"

namespace rovaca
{

pReadRecord ReadRecordUtils::copy_read(pReadRecord const &original, uint8_t *seq, uint8_t *qual, uint32_t seq_len, pCigar cigar,
                                       pBamDataPool bampool, pMemoryPool pool)
{
    sam_hdr_t *header = original->header();
    bam1_t *new_bam = bampool->alloc_bam_struct_pool();

    uint16_t qname_len = original->qname_len();
    const char *qname = original->qname();
    uint16_t flag = original->flag();
    int32_t tid = original->get_tid(), mtid = original->mate_tid();
    int64_t pos = std::max(original->get_start() - 1, (int64_t)0), mpos = original->mate_pos();
    uint8_t mapping_quality = original->mapping_quality();
    int64_t isize = original->insert_size();

    int32_t r = bam_set1(new_bam, qname_len, qname, flag, tid, pos + 1, mapping_quality, cigar->num, cigar->data, mtid, mpos, isize,
                         seq_len, (char *)seq, (char *)qual, 0);
    if (ROVACA_UNLIKELY(r < 0)) {
        RovacaLogger::error("bam_set1 error");
        exit(EXIT_FAILURE);
    }

    // bampool->assign_bam_mem_pool(new_bam);
    bampool->peek_bam(new_bam);

    return ReadRecord::create(pool, header, new_bam);
}

pReadRecord ReadRecordUtils::copy_read(pReadRecord const &original, pBamDataPool bampool, pMemoryPool gpool)
{
    pBases seq = original->decode_to_str(gpool);
    uint8_t *qual = original->qual();

    sam_hdr_t *header = original->header();
    bam1_t *new_bam = bampool->alloc_bam_struct_pool();

    uint16_t qname_len = original->qname_len();
    const char *qname = original->qname();
    uint16_t flag = original->flag();
    int32_t tid = original->get_tid(), mtid = original->mate_tid();
    int64_t pos = std::max(original->get_start() - 1, (int64_t)0), mpos = original->mate_pos();
    uint8_t mapping_quality = original->mapping_quality();
    int64_t isize = original->insert_size();

    int32_t r = bam_set1(new_bam, qname_len, qname, flag, tid, pos + 1, mapping_quality, original->cigar_length(), original->cigar(), mtid,
                         mpos, isize, seq->num, (char *)seq->data, (char *)qual, 0);
    if (ROVACA_UNLIKELY(r < 0)) {
        RovacaLogger::error("bam_set1 error");
        exit(EXIT_FAILURE);
    }

    // bampool->assign_bam_mem_pool(new_bam);
    bampool->peek_bam(new_bam);

    return ReadRecord::create(gpool, header, new_bam);
}

Int64Uint32Pair ReadRecordUtils::get_read_index_for_reference_coordinate(int64_t alignment_start, uint32_t cigar_num, const uint32_t *cigar,
                                                                         int64_t ref_coord)
{
    if (ref_coord < alignment_start) {
        return {s_read_index_not_found, BAM_CUNINITIALIZE};
    }

    int64_t first_read_pos_of_element, first_ref_pos_of_element;                      // inclusive
    int64_t last_read_pos_of_element = 0, last_ref_pos_of_element = alignment_start;  // exclusive

    // 遍历cigar，直到某个CIgarElement将ref_coord括起来
    uint32_t i, ce, op, op_len;
    for (i = 0; i < cigar_num; ++i) {
        ce = cigar[i];
        op = bam_cigar_op(ce);
        op_len = bam_cigar_oplen(ce);
        first_read_pos_of_element = last_read_pos_of_element;
        first_ref_pos_of_element = last_ref_pos_of_element;
        last_read_pos_of_element += consumes_read_bases(op) ? op_len : 0;
        last_ref_pos_of_element += (consumes_ref_bases(op) || cigar_op_is_soft_clip(op)) ? op_len : 0;

        // ref_coord falls within this cigar element
        if (first_ref_pos_of_element <= ref_coord && ref_coord < last_ref_pos_of_element) {
            int64_t pos = first_read_pos_of_element + (consumes_read_bases(op) ? (ref_coord - first_ref_pos_of_element) : 0);
            return {pos, op};
        }
    }
    return {s_read_index_not_found, BAM_CUNINITIALIZE};
}

Int64Uint32Pair ReadRecordUtils::get_read_index_for_reference_coordinate(pReadRecord const &read, int64_t ref_coord)
{
    return get_read_index_for_reference_coordinate(read->get_soft_start(), read->cigar_length(), read->cigar(), ref_coord);
}

std::vector<std::string> ReadRecordUtils::get_samples_from_header(const std::vector<bam_hdr_t *> &headers)
{
    std::vector<std::string> samples;
    const char *q, *p;
    for (const auto &h : headers) {
        p = strstr(strstr(h->text, "@RG"), "SM");
        for (p += 3, q = p; *q && *q != '\t' && *q != '\n'; ++q) {}
        std::string sample_name(p, q - p);
        if (std::find(samples.begin(), samples.end(), sample_name) == samples.end()) {
            samples.push_back(std::move(sample_name));
        }
    }
    return samples;
}

OptionalUint8 ReadRecordUtils::get_read_base_quality_at_reference_coordinate(pReadRecord read, int64_t ref_coord)
{
    if (ref_coord < read->get_start() || read->get_stop() < ref_coord) {
        return {false, 0};
    }

    int64_t soft_start = read->get_soft_start();
    uint32_t cigar_num = read->cigar_length();
    uint32_t *cigar = read->cigar();
    Int64Uint32Pair offset_and_operator = get_read_index_for_reference_coordinate(soft_start, cigar_num, cigar, ref_coord);

    if (offset_and_operator.second != BAM_CUNINITIALIZE && consumes_read_bases(offset_and_operator.second)) {
        return {true, read->qual_i((uint32_t)offset_and_operator.first)};
    }
    return {false, 0};
}

}  // namespace rovaca