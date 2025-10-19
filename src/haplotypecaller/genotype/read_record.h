#ifndef ROVACA_HC_READ_RECORD_H_
#define ROVACA_HC_READ_RECORD_H_
#include <string>

#include "forward.h"
#include "genotype_macors.h"
#include "htslib/sam.h"
#include "interface/interface_locatable.hpp"

typedef struct hc_apply_one_read_t hc_apply_one_read, *p_hc_apply_one_read;

namespace rovaca
{

static constexpr int32_t s_cannot_compute_adaptor_boundary = INT32_MIN;

typedef struct RealignedResult
{
    pCigar realigned_cigar{nullptr};
    int64_t new_start{INVALID_INT};
    int64_t new_stop{INVALID_INT};
    int32_t _new_adaptor_bounder{INVALID_INT};
} RealignedResult, *pRealignedResult;

class ReadRecord : public InterfaceLocatable
{
    /* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
private:
    pMemoryPool _pool;
    bam1_t *_raw_bam{nullptr};
    bam_hdr_t *_header{nullptr};
    p_hc_apply_one_read _assemble_bam{nullptr};
    int64_t _start{INVALID_INT}, _stop{INVALID_INT}, _soft_start{INVALID_INT}, _soft_stop{INVALID_INT};
    int32_t _consumes_query{INVALID_INT}, _consumes_reference{INVALID_INT};
    int32_t _adaptor_bounder{INVALID_INT};
    int32_t _sample_index : 31;
    int32_t _is_native : 1;
    pInformativeSet _informative_set{nullptr};

    /* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
public:
    pRealignedResult _realigned_result{nullptr};

    /* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
public:
    ~ReadRecord() override = default;
    int32_t get_tid() const override;
    int64_t get_start() const override;
    int64_t get_stop() const override;
    void set_tid(int32_t tid) override;
    void set_start(int64_t start) override;
    void set_stop(int64_t stop) override;

    /* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
public:
    DISALLOW_COPY_AND_ASSIGN(ReadRecord);

    /*!
     * @brief 对象的唯一生成接口，确保所有对象均源于pMemoryPool
     * @note ReadRecord可以由两种源构建：htslib中的bam1_t，或是assemble产出的p_hc_apply_one_read
     */
    static pReadRecord create(pMemoryPool pool, bam_hdr_t *header, bam1_t *hts_read);
    static pReadRecord create(pMemoryPool pool, bam_hdr_t *header, p_hc_apply_one_read assemble_read);

    uint16_t flag() const;
    uint32_t *cigar() const;
    int32_t seq_length() const;
    uint32_t cigar_length() const;
    int32_t mate_tid() const;
    int64_t mate_pos() const;
    int64_t insert_size() const;
    uint8_t mapping_quality() const;
    uint8_t *qual() const;
    uint8_t seq_i(uint32_t i) const;
    uint8_t qual_i(uint32_t i) const;
    uint32_t cigar_i(uint32_t i) const;
    const char *qname() const;
    uint16_t qname_len() const;
    bool is_empty() const { return seq_length() == 0; }
    bam_hdr_t *header() const { return _header; }
    int32_t sample_index() const { return _sample_index; }
    p_hc_apply_one_read assemble_read() { return _assemble_bam; }
    void set_sample_index(int32_t si) { _sample_index = si; }
    bool mate_on_same_contig_or_no_mapped_mate_read() const;
    bam1_t *raw_data() const { return _raw_bam; }

    int64_t get_soft_start() const { return const_cast<ReadRecord *>(this)->get_soft_start(); }
    int64_t get_soft_end() const { return const_cast<ReadRecord *>(this)->get_soft_end(); }
    int64_t get_soft_start();
    int64_t get_soft_end();

    /*!
     * @brief 计算一个read序列的长度减去软剪切（soft clip）的碱基数
     */
    int32_t unclipped_read_length() const;

    /*! @brief 此处仅供PairHMM转换为普通字符串时使用 */
    pBases decode_to_str(pMemoryPool pool) const;
    void decode_base(pBases base) const;
    void decode_qual(pBases base) const;
    void ins_gops(pBases base) const;
    void del_gops(pBases base) const;
    void gap_conts(pBases base) const;

    /*! @brief 获取当前 reads 的 adaptor_bounder */
    int32_t get_adaptor_bounder();

    pInformativeSet get_informative_set() const { return _informative_set; }
    void set_informative_set(pInformativeSet is) { _informative_set = is; }
    void clear_informative_set()
    {
        if (_informative_set) {
            using namespace boost;
            _informative_set->~dynamic_bitset<>();
            _informative_set = nullptr;
        }
    }

    /*! @brief 当前reads是不是成对的 */
    bool is_paired() const;
    /*! @brief 当前reads和mate-reads是不是一组合理范围的reads */
    bool is_properly_paired() const;
    /*! @brief 当前reads是不是未比对上。（用不着，前面的筛掉了） */
    bool is_unmapped() const;
    /*! @brief 配对reads是不是未比对上。（用不着，前面的筛掉了） */
    bool mate_is_unmapped() const;
    /*! @brief 当前reads是否是正负链 */
    bool is_reverse_strand() const;
    /*! @brief 配对reads是否是正负链 */
    bool mate_is_reverse_strand() const;
    /*! @brief 当前reads 是不是Pair-reads的第一条reads */
    bool is_second_of_pair() const;
    /*! @brief 当前reads 是不是Pair-reads的第二条reads */
    bool is_first_of_pair() const;
    /*! @brief 当前比对的结果为reads的次优的结果 */
    bool is_secondary_alignment() const;
    /*! @brief 当前reads 没有过quality的检测 */
    bool fails_vendor_quality_check() const;
    /*! @brief 当前reads是不是PCR-Duplicate的reads（用不着，前面的筛掉了） */
    bool is_duplicate() const;
    /*! @brief 此reads的比对结果是补充的位置
     * 一条拆成两部分:一部分在primary上(离pair-reads比较近),另外一条就是supplementary_alignment */
    bool is_supplementary_alignment() const;

    /*! @brief 构造此Read的bit flag(每位代表一种状态)————此flag和Read中自身的flag有什么区别？ */
    int32_t get_sam_flags_for_read() const;

    /** the query sequence itself is unmapped */
    void set_read_unmapped_flag(bool flag);

    /* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
private:
    void init_consumes();

    /*! @brief 使用bam1_t构建 */
    ReadRecord(bam1_t *b, bam_hdr_t *header, pMemoryPool pool);
    /*! @brief 使用组装后的read构建 */
    ReadRecord(p_hc_apply_one_read b, bam_hdr_t *header, pMemoryPool pool);

    void set_flag(bool flag, uint16_t bit)
    {
        if (flag)
            _raw_bam->core.flag |= bit;
        else
            _raw_bam->core.flag &= ~bit;
    }

    /*! @brief 根据read和其mate的比对信息判断是否可以去除read中的适配器序列 */
    bool has_well_defined_fragment_size() const;
    int32_t get_adaptor_boundary(int);
};

}  // namespace rovaca

namespace std
{

template <>
struct less<rovaca::pReadRecord>
{
    bool operator()(rovaca::pReadRecord l, rovaca::pReadRecord r) const
    {
        if (l->get_tid() != r->get_tid()) {
            return l->get_tid() < r->get_tid();
        }
        if (l->get_start() != r->get_start()) {
            return l->get_start() < r->get_start();
        }
        if (l->is_reverse_strand() != r->is_reverse_strand()) {
            return !l->is_reverse_strand();
        }
        if (strcmp(l->qname(), r->qname()) != 0) {
            return strcmp(l->qname(), r->qname());
        }
        if (l->flag() != r->flag()) {
            return l->flag() < r->flag();
        }
        if (l->mapping_quality() != r->mapping_quality()) {
            return l->mapping_quality() < r->mapping_quality();
        }
        if (l->is_paired() && r->is_paired()) {
            if (l->mate_tid() != r->mate_tid()) {
                return l->mate_tid() < r->mate_tid();
            }
            if (l->mate_pos() != r->mate_pos()) {
                return l->mate_pos() < r->mate_pos();
            }
        }
        return l->insert_size() < r->insert_size();
    }
};

}  // namespace std

#endif  // ROVACA_HC_READ_RECORD_H_
