#ifndef ROVACA_HC_VARIANT_H_
#define ROVACA_HC_VARIANT_H_
#include "forward.h"
#include "genotype_enum.h"
#include "genotype_macors.h"
#include "htslib/vcf.h"
#include "interface/interface_locatable.hpp"

namespace rovaca
{

class Variant : public InterfaceLocatable
{
    /* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
private:
    AlleleVector _alleles;
    AlleleVector _alt_alleles;
    int64_t _start{INVALID_INT}, _stop{INVALID_INT};
    double _log10_error{NEGATIVE_INFINITY};
    pAllele _ref{nullptr}, _biallelic_alt{nullptr};
    pBases _db_id{nullptr};
    pInfoData _info{nullptr};
    pGenotypesContext _genotypes{nullptr};
    int32_t _tid{INVALID_INT}, _source_id{INVALID_INT}, _block_end{INVALID_INT};
    VariantType _type{VT_UNINITIALIZED}, _type_ignoring_non_ref{VT_UNINITIALIZED};

    /* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
public:
    DISALLOW_COPY_AND_ASSIGN(Variant);
    ~Variant() override = default;
    int32_t get_tid() const override { return _tid; }
    int64_t get_start() const override { return _start; }
    int64_t get_stop() const override { return _stop; }
    void set_tid(int32_t tid) override { _tid = tid; }
    void set_start(int64_t start) override;
    void set_stop(int64_t stop) override;

    /*! @brief 对象的唯一生成接口，确保所有对象均源于pMemoryPool */
    static pVariant create(pMemoryPool pool) { return new ALLOC_TYPE_IN_POOL(pool, Variant) Variant{pool}; }
    static pVariant create_dbsnp_source(bcf1_t* b, pMemoryPool pool);

    VariantType type(bool ignore_non_ref) const;
    VariantType type() const { return type(false); }
    const AlleleVector& alleles() const { return _alleles; }
    const AlleleVector& alt_alleles() const;
    pAllele ref_allele() const { return _ref; }
    pAllele biallelic_alt() const { return _biallelic_alt; }
    pAllele alternate_allele_at(uint32_t i) const { return _alleles.at(i + 1); }  // first is ref
    int32_t end_key() const { return _block_end; }
    pBases db_id() const { return _db_id; }
    pInfoData info() const { return _info; }
    pGenotypesContext genotype() { return _genotypes; }
    double log10_error() const { return _log10_error; }
    size_t allele_num() const { return _alleles.size(); }
    int32_t source_id() const { return _source_id; }
    size_t sample_count() const;
    int32_t get_max_ploidy(int32_t default_ploidy);
    int32_t get_called_chr_count() const;
    int32_t get_called_chr_count(pAllele a) const;
    int32_t get_allele_index(pAllele a) const;
    Int32Vector get_gl_indices_of_alternate_allele(pAllele a, pMemoryPool pool) const;
    double get_phred_scaled_qual() const { return (_log10_error * -10) + 0.0; }

    void set_alleles(const AlleleVector& alleles);
    void set_info(pInfoData info) { _info = info; }
    void set_id(pBases id) { _db_id = id; }
    void set_end_key(int32_t end) { _block_end = end; }
    void set_source_id(int32_t id) { _source_id = id; }
    void set_log10_error(double error) { _log10_error = error; }
    void set_genotype(pGenotypesContext genotypes) { _genotypes = genotypes; }

    /** @return true if the context is strictly bi-allelic */
    bool is_biallelic() const { return 2 == _alleles.size(); }
    bool is_snp() const { return type() == VT_SNP; }
    bool is_indel1() const { return type() == VT_INDEL; }
    bool is_simple_indel() const;
    bool is_simple_insertion() const;
    bool is_simple_deletion() const;
    bool is_variant() const { return type() != VT_NO_VARIATION; }

    void add_non_ref_symbolic_allele();

    static bool has_symbolic_alleles(const AlleleVector& alleles);
    bool has_symbolic_alleles() { return has_symbolic_alleles(_alleles); }
    bool has_allele(pAllele a) const;
    bool has_allele(pAllele a, bool ignore_ref_state, bool consider_ref_allele) const;
    bool has_genotypes() const { return nullptr != _genotypes; }
    bool has_log10_error() const { return _log10_error != NEGATIVE_INFINITY; }
    bool has_end_key() const { return _block_end != INVALID_INT; }
    bool equals(const Variant& v) const;

    /* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
private:
    explicit Variant(pMemoryPool pool)
        : _alleles(pool)
        , _alt_alleles(pool)
    {}

    VariantType determine_type(bool ignore_non_ref);
    VariantType determine_polymorphic_type(bool ignore_non_ref);

    static VariantType type_of_biallelic_variant(pAllele ref, pAllele a);

    /*! @brief 检测Variantstart、stop、alleles是否合格 */
    void validate_stop();
    void validate_alleles();
    void validate_variant()
    {
        validate_stop();
        validate_alleles();
    }
};

struct VariantHash
{
    size_t operator()(const pVariant& v) const;
};

struct VariantEqual
{
    bool operator()(const pVariant& l, const pVariant& r) const;
};

}  // namespace rovaca

#endif  // ROVACA_HC_VARIANT_H_
