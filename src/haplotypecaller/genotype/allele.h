#ifndef ROVACA_HC_ALLELE_H_
#define ROVACA_HC_ALLELE_H_
#include <memory>
#include <string>

#include "forward.h"
#include "genotype_macors.h"
#include "genotype_struct.h"

namespace rovaca
{

/*!
 * @brief StaticAllele 中为静态的 Allele 对象，程序运行期间持续存在
 */
struct StaticAllele
{
    static pConstantsStr s_constants;

    std::unique_ptr<Allele> _ref_a{nullptr};
    std::unique_ptr<Allele> _ref_t{nullptr};
    std::unique_ptr<Allele> _ref_g{nullptr};
    std::unique_ptr<Allele> _ref_c{nullptr};
    std::unique_ptr<Allele> _ref_n{nullptr};
    std::unique_ptr<Allele> _alt_a{nullptr};
    std::unique_ptr<Allele> _alt_t{nullptr};
    std::unique_ptr<Allele> _alt_g{nullptr};
    std::unique_ptr<Allele> _alt_c{nullptr};
    std::unique_ptr<Allele> _alt_n{nullptr};
    std::unique_ptr<Allele> _span_del{nullptr};
    std::unique_ptr<Allele> _no_call{nullptr};
    std::unique_ptr<Allele> _non_ref_allele{nullptr};
    std::unique_ptr<Allele> _unspecified_alternate_allele{nullptr};
    std::unique_ptr<Allele> _spanning_deletion_symbolic_allele_deprecated{nullptr};

    DISALLOW_COPY_AND_ASSIGN(StaticAllele);
    ~StaticAllele() = default;

    static StaticAllele* get_instance()
    {
        static std::unique_ptr<StaticAllele> ua(new StaticAllele());
        return ua.get();
    }

protected:
    StaticAllele();
};

typedef struct StaticAllele StaticAllele, *pStaticAllele;

/*!
 * @brief Allele 类表示一个等位基因，即一个基因座上的一种变异形式
 */
class Allele
{
    /* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
private:
    friend struct StaticAllele;

    static pConstantsStr s_constants;  // 初始化必须优先于s_static_obj
    static pStaticAllele s_static_obj;

    /* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
private:
    pBases _bases;
    uint8_t _is_ref : 1, _is_no_call : 1, _is_symbolic : 1, : 5;

    /* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
public:
    DISALLOW_COPY_AND_ASSIGN(Allele);
    virtual ~Allele() = default;

    /*! @brief 对象的唯一生成接口，确保所有对象均源于pMemoryPool */
    static pAllele create_allele(uint8_t bases, uint8_t is_ref);
    static pAllele create_allele(pBases bases, uint8_t is_ref, pMemoryPool pool);

    /*! @brief 外部手动构造 pBases 有点蠢，升级为可传入char* */
    static pAllele create_allele(const char* bases_str, uint8_t is_ref, pMemoryPool pool);
    static pAllele create_allele(const char* bases_str, uint32_t num, uint8_t is_ref, pMemoryPool pool);

    /*! @brief 判断该变异位点是否为参考等位基因 */
    bool is_reference() const { return _is_ref; }
    /*! @brief 判断该变异位点是否被召回（即是否被成功检测到） */
    bool is_called() const { return !_is_no_call; }
    /*! @brief 判断该变异位点是否为符号性变异（Insertion/Deletion/Inversion/Translocation/Repeat） */
    bool is_symbolic() const { return _is_symbolic; }

    /** @return true if this Allele is a breakpoint ( ex: G]17:198982] or ]13:123456]A ) */
    bool is_breakpoint() const;
    /** @return true if this Allele is a single breakend (ex: .A or A.) */
    bool is_single_breakend() const;
    /** @return true if Allele is either <NON_REF> or <*> */
    bool is_non_ref_allele() const;

    /*! @brief Returns true if this and other are equal. If ignore_ref_state is true, doesn't require both alleles has the same ref tag */
    bool equals(const Allele& other, bool ignore_ref_state) const;
    bool equals(const Allele& other) const { return equals(other, false); }

    /*!
     * @brief Return the DNA bases segregating in this allele.  Note this isn't reference polarized, so the Null allele is represented by a
     * vector of length 0
     */
    pBases get_bases() const;
    pBases get_display_string() const { return _bases; }

    /*! @return the length of this allele.  Null and NO_CALL alleles have 0 length */
    uint32_t length() const { return _is_symbolic ? 0 : _bases->num; }

    size_t hash() const;

    void init_allele(pBases bases, uint8_t is_ref);

    /* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
protected:
    /*! @brief 此函数供Haplotype调用 */
    Allele()
        : _bases(nullptr)
        , _is_ref(0)
        , _is_no_call(0)
        , _is_symbolic(0)
    {}

    /* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
private:
    /*!
     * @brief 此函数仅供StaticAllele中创建静态对象时调用
     * @param bases
     * @param is_ref
     */
    Allele(pBases bases, uint8_t is_ref)
        : _bases(nullptr)
        , _is_ref(0)
        , _is_no_call(0)
        , _is_symbolic(0)
    {
        init_allele(bases, is_ref);
    }

    static bool would_be_star_allele(pBases bases);
    static bool would_be_null_allele(pBases bases);
    static bool would_be_no_call_allele(pBases bases);
    static bool would_be_symbolic_allele(pBases bases);
    static bool would_be_breakpoint_allele(pBases bases);
    static bool would_be_single_breakend_allele(pBases bases);
    static bool acceptable_allele_bases(pBases bases, bool is_ref);
};

struct AlleleHash
{
    size_t operator()(const pAllele& a) const;
};

struct AlleleEqual
{
    bool operator()(const pAllele& l, const pAllele& r) const;
};

bool less_bases(pBases l_bases, pBases r_bases);

}  // namespace rovaca

namespace std
{

template <>
struct less<rovaca::pAllele>
{
    bool operator()(rovaca::pAllele l, rovaca::pAllele r) const
    {
        if (l == r) return false;
        if (l->is_reference() && !r->is_reference())
            return true;
        else if (!l->is_reference() && r->is_reference())
            return false;
        rovaca::pBases l_bases = l->get_display_string();
        rovaca::pBases r_bases = r->get_display_string();
        return rovaca::less_bases(l_bases, r_bases);
    }
};

}  // namespace std

#endif  // ROVACA_HC_ALLELE_H_
