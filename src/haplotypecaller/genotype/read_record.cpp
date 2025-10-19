#include "read_record.h"

#include "assemble_interface.h"
#include "genotype_struct.h"
#include "htslib/sam.h"
#include "utils/cigar_utils.h"

namespace rovaca
{

#define PAIRED(b)                     (((b)->core.flag & BAM_FPAIRED) != 0)
#define PROPERLY_PAIRED(b)            (((b)->core.flag & BAM_FPROPER_PAIR) != 0)
#define UNMAPPED(b)                   (((b)->core.flag & BAM_FUNMAP) != 0)
#define MATE_UNMAPPED(b)              (((b)->core.flag & BAM_FMUNMAP) != 0)
#define REVERSE_STRAND(b)             (((b)->core.flag & BAM_FREVERSE) != 0)
#define MATE_REVERSE_STRAND(b)        (((b)->core.flag & BAM_FMREVERSE) != 0)
#define FIRST_OF_PAIR(b)              (((b)->core.flag & BAM_FREAD1) != 0)
#define SECOND_OF_PAIR(b)             (((b)->core.flag & BAM_FREAD2) != 0)
#define SECONDARY_ALIGNMENT(b)        (((b)->core.flag & BAM_FSECONDARY) != 0)
#define FAILS_VENDOR_QUALITY_CHECK(b) (((b)->core.flag & BAM_FQCFAIL) != 0)
#define DUPLICATE(b)                  (((b)->core.flag & BAM_FDUP) != 0)
#define SUPPLEMENTARY_ALIGNMENT(b)    (((b)->core.flag & BAM_FSUPPLEMENTARY) != 0)

int32_t ReadRecord::get_tid() const { return _raw_bam->core.tid; }

int64_t ReadRecord::get_start() const
{
    if (_realigned_result) {
        return _realigned_result->new_start;
    }
    else if (INVALID_INT == _start) {
        const_cast<ReadRecord*>(this)->_start = _is_native ? _raw_bam->core.pos : _assemble_bam->pos_start;
    }
    return _start;
}

int64_t ReadRecord::get_stop() const
{
    if (_realigned_result) {
        return _realigned_result->new_stop;
    }
    else if (INVALID_INT == _stop) {
        if (INVALID_INT == _consumes_reference) {
            const_cast<ReadRecord*>(this)->init_consumes();
        }
        const_cast<ReadRecord*>(this)->_stop = get_start() + _consumes_reference - 1;
    }
    return _stop;
}

void ReadRecord::set_tid(int32_t tid) { _raw_bam->core.tid = tid; }

void ReadRecord::set_start(int64_t start) { _start = start; }

void ReadRecord::set_stop(int64_t stop) { _stop = stop; }

pReadRecord ReadRecord::create(pMemoryPool pool, bam_hdr_t* header, bam1_t* hts_read)
{
    return new ALLOC_TYPE_IN_POOL(pool, ReadRecord) ReadRecord{hts_read, header, pool};
}

pReadRecord ReadRecord::create(pMemoryPool pool, bam_hdr_t* header, p_hc_apply_one_read assemble_read)
{
    return new ALLOC_TYPE_IN_POOL(pool, ReadRecord) ReadRecord{assemble_read, header, pool};
}

uint16_t ReadRecord::flag() const { return _raw_bam->core.flag; }

uint32_t* ReadRecord::cigar() const
{
    return _realigned_result ? _realigned_result->realigned_cigar->data : (_is_native ? bam_get_cigar(_raw_bam) : _assemble_bam->cigar);
}

int32_t ReadRecord::seq_length() const { return _is_native ? _raw_bam->core.l_qseq : (int32_t)_assemble_bam->read_len; }

uint32_t ReadRecord::cigar_length() const
{
    return _realigned_result ? _realigned_result->realigned_cigar->num : (_is_native ? _raw_bam->core.n_cigar : _assemble_bam->cigar_len);
}

int32_t ReadRecord::mate_tid() const { return _raw_bam->core.mtid; }

int64_t ReadRecord::mate_pos() const { return _raw_bam->core.mpos; }

int64_t ReadRecord::insert_size() const { return _is_native ? _raw_bam->core.isize : _assemble_bam->insert_size; }

uint8_t ReadRecord::mapping_quality() const { return _raw_bam->core.qual; }

uint8_t* ReadRecord::qual() const { return _is_native ? bam_get_qual(_raw_bam) : _assemble_bam->qual; }

uint8_t ReadRecord::seq_i(uint32_t i) const
{
    return _is_native ? seq_nt16_str[bam_seqi(bam_get_seq(_raw_bam), i)] : _assemble_bam->seq[i];
}

uint8_t ReadRecord::qual_i(uint32_t i) const { return _is_native ? bam_get_qual(_raw_bam)[i] : _assemble_bam->qual[i]; }

uint32_t ReadRecord::cigar_i(uint32_t i) const { return cigar()[i]; }

const char* ReadRecord::qname() const { return bam_get_qname(_raw_bam); }

uint16_t ReadRecord::qname_len() const { return _raw_bam->core.l_qname; }

bool ReadRecord::mate_on_same_contig_or_no_mapped_mate_read() const
{
    return !is_paired() || is_unmapped() || (!is_unmapped() && get_tid() == mate_tid());
}

int64_t ReadRecord::get_soft_start()
{
    if (_soft_start == INVALID_INT) {
        _soft_start = get_start();
        uint32_t i = 0, ce, op, op_len, cigar_len = cigar_length();
        for (; i < cigar_len; ++i) {
            ce = cigar_i(i);
            op = bam_cigar_op(ce);
            op_len = bam_cigar_oplen(ce);
            if (cigar_op_is_soft_clip(op)) {
                _soft_start -= op_len;
            }
            else if (!cigar_op_is_hard_clip(op)) {
                break;
            }
        }
    }
    return _soft_start;
}

int64_t ReadRecord::get_soft_end()
{
    if (_soft_stop == INVALID_INT) {
        bool found_aligned_base = false;
        _soft_stop = get_stop();
        uint32_t ce, op, op_len, cigar_len = cigar_length();
        uint32_t i = cigar_len - 1;
        for (; i < cigar_len; --i) {  // uint32负数越界
            ce = cigar_i(i);
            op = bam_cigar_op(ce);
            op_len = bam_cigar_oplen(ce);
            if (cigar_op_is_soft_clip(op)) {
                _soft_stop += op_len;
            }
            else if (!cigar_op_is_hard_clip(op)) {
                found_aligned_base = true;
                break;
            }
        }
        if (!found_aligned_base) {
            _soft_stop = get_stop();
        }
    }
    return _soft_stop;
}

int32_t ReadRecord::unclipped_read_length() const
{
    int32_t soft_clipped_bases = 0;
    for (uint32_t i = 0, len = cigar_length(); i < len; ++i) {
        if (cigar_op_is_soft_clip(bam_cigar_op(cigar_i(i)))) {
            soft_clipped_bases += bam_cigar_oplen(cigar_i(i));
        }
    }
    return seq_length() - soft_clipped_bases;
}

pBases ReadRecord::decode_to_str(pMemoryPool pool) const
{
    uint32_t len = (uint32_t)seq_length();
    pBases result = new ALLOC_FLEXIBLE_IN_POOL(pool, Bases, len, uint8_t) Bases{len};
    for (uint32_t i = 0, idx = 0; i < len; i++) {
        result->data[idx++] = seq_i(i);
    }
    return result;
}

void ReadRecord::decode_base(pBases base) const
{
    uint32_t len = base->num = (uint32_t)seq_length();
    for (uint32_t i = 0, idx = 0; i < len; i++) {
        base->data[idx++] = seq_i(i);
    }
}

void ReadRecord::decode_qual(pBases base) const
{
    uint32_t len = base->num = (uint32_t)seq_length();
    base->data[len] = '\0';
    memcpy(base->data, this->qual(), len * sizeof(uint8_t));
}

void ReadRecord::ins_gops(pBases base) const
{
    uint32_t len = base->num = (uint32_t)seq_length();
    base->data[len] = '\0';
    memset(base->data, '-', len);
}

void ReadRecord::del_gops(pBases base) const
{
    uint32_t len = base->num = (uint32_t)seq_length();
    base->data[len] = '\0';
    memset(base->data, '-', len);
}

void ReadRecord::gap_conts(pBases base) const
{
    uint32_t len = base->num = (uint32_t)seq_length();
    base->data[len] = '\0';
    memset(base->data, '+' - 33, len);
}

int32_t ReadRecord::get_adaptor_bounder()
{
    if (_realigned_result && _realigned_result->_new_adaptor_bounder == INVALID_INT) {
        _realigned_result->_new_adaptor_bounder = get_adaptor_boundary(1);
    }
    else if (INVALID_INT == _adaptor_bounder) {
        _adaptor_bounder = get_adaptor_boundary(1);
    }
    return _realigned_result ? _realigned_result->_new_adaptor_bounder : _adaptor_bounder;
}

bool ReadRecord::is_paired() const { return PAIRED(_raw_bam); }

bool ReadRecord::is_properly_paired() const { return PAIRED(_raw_bam) && PROPERLY_PAIRED(_raw_bam); }

bool ReadRecord::is_unmapped() const { return UNMAPPED(_raw_bam); }

bool ReadRecord::mate_is_unmapped() const { return PAIRED(_raw_bam) && MATE_UNMAPPED(_raw_bam); }

bool ReadRecord::is_reverse_strand() const { return REVERSE_STRAND(_raw_bam); }

bool ReadRecord::mate_is_reverse_strand() const { return PAIRED(_raw_bam) && MATE_REVERSE_STRAND(_raw_bam); }

bool ReadRecord::is_second_of_pair() const { return PAIRED(_raw_bam) && SECOND_OF_PAIR(_raw_bam); }

bool ReadRecord::is_first_of_pair() const { return PAIRED(_raw_bam) && FIRST_OF_PAIR(_raw_bam); }

bool ReadRecord::is_secondary_alignment() const { return SECONDARY_ALIGNMENT(_raw_bam); }

bool ReadRecord::fails_vendor_quality_check() const { return FAILS_VENDOR_QUALITY_CHECK(_raw_bam); }

bool ReadRecord::is_duplicate() const { return DUPLICATE(_raw_bam); }

bool ReadRecord::is_supplementary_alignment() const { return SUPPLEMENTARY_ALIGNMENT(_raw_bam); }

int32_t ReadRecord::get_sam_flags_for_read() const
{
    int32_t sam_flags = 0;
    if (is_paired()) {
        sam_flags |= BAM_FPAIRED;
    }
    if (is_properly_paired()) {
        sam_flags |= BAM_FPROPER_PAIR;
    }
    if (is_unmapped()) {
        sam_flags |= BAM_FUNMAP;
    }
    if (is_paired() && mate_is_unmapped()) {
        sam_flags |= BAM_FMUNMAP;
    }
    if (!is_unmapped() && is_reverse_strand()) {
        sam_flags |= BAM_FREVERSE;
    }
    if (is_paired() && !mate_is_unmapped() && mate_is_reverse_strand()) {
        sam_flags |= BAM_FMREVERSE;
    }
    if (is_first_of_pair()) {
        sam_flags |= BAM_FREAD1;
    }
    if (is_second_of_pair()) {
        sam_flags |= BAM_FREAD2;
    }
    if (is_secondary_alignment()) {
        sam_flags |= BAM_FSECONDARY;
    }
    if (fails_vendor_quality_check()) {
        sam_flags |= BAM_FQCFAIL;
    }
    if (is_duplicate()) {
        sam_flags |= BAM_FDUP;
    }
    if (is_supplementary_alignment()) {
        sam_flags |= BAM_FSUPPLEMENTARY;
    }
    return sam_flags;
}

void ReadRecord::set_read_unmapped_flag(bool flag) { set_flag(flag, BAM_FUNMAP); }

void ReadRecord::init_consumes()
{
    _consumes_reference = _consumes_query = 0;
    uint32_t* cigar = this->cigar();
    auto cigar_num = (int32_t)this->cigar_length();
    _consumes_reference = (int32_t)bam_cigar2rlen(cigar_num, cigar);
    _consumes_query = (int32_t)bam_cigar2qlen(cigar_num, cigar);
}

ReadRecord::ReadRecord(bam1_t* b, bam_hdr_t* header, pMemoryPool pool)
    : _pool(pool)
    , _raw_bam(b)
    , _header(header)
    , _sample_index(INVALID_INT)
    , _is_native(true)
{}

ReadRecord::ReadRecord(p_hc_apply_one_read b, bam_hdr_t* header, pMemoryPool pool)
    : _pool(pool)
    , _raw_bam(&b->read)
    , _header(header)
    , _assemble_bam(b)
    , _sample_index(INVALID_INT)
    , _is_native(false)
{
    _consumes_query = static_cast<int32_t>(b->read_len);
    _consumes_reference = static_cast<int32_t>(b->ref_len);
}

bool ReadRecord::has_well_defined_fragment_size() const
{
    if (0 == insert_size() || !is_paired() || (is_unmapped() || mate_is_unmapped()) || (is_reverse_strand() == mate_is_reverse_strand())) {
        return false;
    }
    if (is_reverse_strand()) {
        return get_stop() > mate_pos();
    }
    else {
        return get_start() <= mate_pos() + insert_size();
    }
}

int32_t ReadRecord::get_adaptor_boundary(int)
{
    if (!has_well_defined_fragment_size()) {
        return s_cannot_compute_adaptor_boundary;
    }
    else if (is_reverse_strand()) {
        return int32_t(mate_pos() - 1);
    }
    else {
        auto ins_size = int32_t(std::abs(insert_size()));
        return int32_t(get_start()) + ins_size;
    }
}

}  // namespace rovaca