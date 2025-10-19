#ifndef ROVACA_HC_GENOTYPE_H_
#define ROVACA_HC_GENOTYPE_H_
#include "forward.h"
#include "genotype_enum.h"
#include "genotype_macors.h"
#include "htslib/vcf.h"

namespace rovaca
{

typedef struct PhasedData
{
    pBases pid;
    const char* pgt;
    int32_t ps;
} PhasedData, *pPhasedData;

class Genotype
{
    /* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
private:
    AlleleVector _alleles;
    Int32Vector _ad;
    Int32Vector _pl;
    DoubleVector _gl;
    Int32Vector _sb;
    pPhasedData _phased_data{nullptr};  // phased 为 true 时才会有数据
    int32_t _dp{INVALID_INT}, _gq{INVALID_INT}, _min_dp{INVALID_INT};
    int32_t _phased : 1, _type : 8, _sample_id : 23;

    /* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
public:
    DISALLOW_COPY_AND_ASSIGN(Genotype);
    ~Genotype() = default;

    /*! @brief 对象的唯一生成接口，确保所有对象均源于pMemoryPool */
    static pGenotype create(pMemoryPool pool) { return new ALLOC_TYPE_IN_POOL(pool, Genotype) Genotype{pool}; }

    void set_pl(Int32Vector&& pls) { _pl = std::move(pls); }

    void set_alleles(AlleleVector&& allele) { _alleles = std::move(allele); }
    void set_ad(Int32Vector&& ad) { _ad = std::move(ad); }
    void set_sb(Int32Vector&& sb) { _sb = std::move(sb); }
    void set_dp(int32_t dp) { _dp = dp; }
    void set_gq(int32_t gq) { _gq = gq; }
    void set_id(int32_t id) { _sample_id = id; }
    void set_log10perror(double p_log10error) { _gq = (int32_t)std::round(p_log10error * -10); }
    void set_phased(pPhasedData data)
    {
        _phased = true;
        _phased_data = data;
    }

    const AlleleVector& alleles() const { return _alleles; }
    const Int32Vector& ad() const { return _ad; }
    const Int32Vector& pl() const { return _pl; }
    const Int32Vector& sb() const { return _sb; }
    int32_t sample_id() const { return _sample_id; }
    pAllele get_allele_at(uint32_t i) const { return _alleles.at(i); }
    pGenotypeLikelihoods get_likelihoods() const;

    bool has_likelihoods() const { return !_pl.empty(); }
    bool has_dp() const { return _dp != INVALID_INT; }
    bool has_gq() const { return _gq != INVALID_INT; }
    bool has_min_dp() const { return _min_dp != INVALID_INT; }
    bool has_sb() const { return !_sb.empty(); }
    bool has_ad() const { return !_ad.empty(); }

    /**
     * what is the ploidy of this sample?
     * @return the ploidy of this genotype.  0 if the site is no-called.
     */
    int32_t get_ploidy() const { return (int32_t)_alleles.size(); }
    int32_t get_dp() const { return _dp; }
    int32_t get_gq() const { return _gq; }
    int32_t get_min_dp() const { return _min_dp; }
    GenotypeType get_type();
    int32_t count_allele(pAllele a) const;

    bool is_phased() const { return _phased; }
    /** @return true if all observed alleles are the same (regardless of whether they are ref or alt); if any alleles are no-calls, this
     * method will return false. */
    bool is_hom() { return is_hom_ref() || is_hom_var(); }
    /** @return true if all observed alleles are ref; if any alleles are no-calls, this method will return false. */
    bool is_hom_ref() { return get_type() == GT_HOM_REF; }
    /** @return true if all observed alleles are alt; if any alleles are no-calls, this method will return false. */
    bool is_hom_var() { return get_type() == GT_HOM_VAR; }
    /** @return true if we're het (observed alleles differ); if the ploidy is less than 2 or if any alleles are no-calls, this method will
     * return false. */
    bool is_het() { return get_type() == GT_HET; }
    /** @return true if we're het (observed alleles differ) and neither allele is reference; if the ploidy is less than 2 or if any alleles
     * are no-calls, this method will return false. */
    bool is_het_non_ref();
    /** @return true if this genotype is not actually a genotype but a "no call" (e.g. './.' in vcf); if any alleles are not no-calls (even
     * if some are), this method will return false. */
    bool is_no_call() { return get_type() == GT_NO_CALL; }
    /** @return true if this genotype is comprised of any alleles that are not no-calls (even if some are). */
    bool is_called() { return !is_no_call() && get_type() != GT_UNAVAILABLE; }
    /** @return true if this genotype is comprised of both calls and no-calls. */
    bool is_mixed() { return get_type() == GT_MIXED; }
    /** @return true if the type of this genotype is set. */
    bool is_available() { return get_type() != GT_UNAVAILABLE; }

    std::string& append_genotype_data(std::string& s, const std::pmr::map<pAllele, char>& mm);

    void genotype2bcf(bcf_hdr_t* vcf_header, bcf1_t* rec, const std::pmr::map<pAllele, int32_t>& mm);

    /* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
private:
    explicit Genotype(pMemoryPool pool);

    GenotypeType determine_type();
};

}  // namespace rovaca

#endif  // ROVACA_HC_GENOTYPE_H_
