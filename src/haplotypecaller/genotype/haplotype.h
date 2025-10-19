#ifndef ROVACA_HC_HAPLOTYPE_H_
#define ROVACA_HC_HAPLOTYPE_H_
#include "allele.h"

namespace rovaca
{

class Haplotype : public Allele
{
    /* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
private:
    pCigar _cigar{nullptr};
    pSimpleInterval _interval{nullptr};
    pEventMap _events{nullptr};
    double _score{NEGATIVE_INFINITY};
    int64_t _alignment_start_hap_wrt_ref{INVALID_INT};
    int32_t _kmer_size{INVALID_INT};

    /* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
public:
    ~Haplotype() override = default;

    /*! @brief 对象的唯一生成接口，确保所有对象均源于pMemoryPool */
    static pHaplotype create(pMemoryPool pool);

    pCigar cigar() const { return _cigar; }
    double score() const { return _score; }
    int64_t alignment_start_hap_wrt_ref() const { return _alignment_start_hap_wrt_ref; }
    int32_t kmer_size() const { return _kmer_size; }
    pEventMap event_map() const { return _events; }
    pSimpleInterval interval() const { return _interval; }

    void set_cigar(pCigar c) { _cigar = c; }
    void set_score(double s) { _score = s; }
    void set_alignment_start_hap_wrt_ref(int64_t a) { _alignment_start_hap_wrt_ref = a; }
    void set_kmer_size(int32_t k) { _kmer_size = k; }
    void set_event_map(pEventMap e) { _events = e; }
    void set_interval(pSimpleInterval i) { _interval = i; }

    void init_haplotype(pBases bases, uint8_t is_ref) { init_allele(bases, is_ref); }

    /*! @brief 升级版初始化方式，不需要再手动构建 pBases */
    void init_haplotype(const char *bases_str, uint8_t is_ref, pMemoryPool pool);
    void init_haplotype(const char *bases_str, uint32_t num, uint8_t is_ref, pMemoryPool pool);

    /* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
protected:
    Haplotype() = default;
};

}  // namespace rovaca

#endif  // ROVACA_HC_HAPLOTYPE_H_
