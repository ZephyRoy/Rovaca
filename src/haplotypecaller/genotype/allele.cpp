#include "allele.h"

#include "constants_str.hpp"
#include "rovaca_logger.h"

constexpr char ALLELE_STR[] = "ACGTNacgtn";

namespace rovaca
{

pConstantsStr StaticAllele::s_constants = ConstantsStr::get_instance();

StaticAllele::StaticAllele()
{
    _ref_a.reset(new Allele(s_constants->k_bases_a, 1));
    _ref_t.reset(new Allele(s_constants->k_bases_t, 1));
    _ref_g.reset(new Allele(s_constants->k_bases_g, 1));
    _ref_c.reset(new Allele(s_constants->k_bases_c, 1));
    _ref_n.reset(new Allele(s_constants->k_bases_n, 1));
    _alt_a.reset(new Allele(s_constants->k_bases_a, 0));
    _alt_t.reset(new Allele(s_constants->k_bases_t, 0));
    _alt_g.reset(new Allele(s_constants->k_bases_g, 0));
    _alt_c.reset(new Allele(s_constants->k_bases_c, 0));
    _alt_n.reset(new Allele(s_constants->k_bases_n, 0));
    _span_del.reset(new Allele(s_constants->k_bases_span_del, 0));
    _no_call.reset(new Allele(s_constants->k_bases_no_call, 0));
    _non_ref_allele.reset(new Allele(s_constants->k_bases_non_ref_allele, 0));
    _unspecified_alternate_allele.reset(new Allele(s_constants->k_bases_unspecified_alternate_allele, 0));
    _spanning_deletion_symbolic_allele_deprecated.reset(new Allele(s_constants->k_spanning_deletion_symbolic_allele_deprecated, 0));
}

pConstantsStr Allele::s_constants = ConstantsStr::get_instance();
pStaticAllele Allele::s_static_obj = StaticAllele::get_instance();

pAllele Allele::create_allele(uint8_t bases, uint8_t is_ref)
{
    switch (bases) {
        case '.': {
            CHECK_CONDITION_EXIT(is_ref, "cannot tag a no_call allele as the reference allele");
            return s_static_obj->_no_call.get();
        }
        case '*': {
            CHECK_CONDITION_EXIT(is_ref, "cannot tag a spanning deletions allele as the reference allele");
            return s_static_obj->_span_del.get();
        }
        case 'A':
        case 'a': return is_ref ? s_static_obj->_ref_a.get() : s_static_obj->_alt_a.get();
        case 'C':
        case 'c': return is_ref ? s_static_obj->_ref_c.get() : s_static_obj->_alt_c.get();
        case 'G':
        case 'g': return is_ref ? s_static_obj->_ref_g.get() : s_static_obj->_alt_g.get();
        case 'T':
        case 't': return is_ref ? s_static_obj->_ref_t.get() : s_static_obj->_alt_t.get();
        case 'N':
        case 'n': return is_ref ? s_static_obj->_ref_n.get() : s_static_obj->_alt_n.get();
        default: break;
    }
    RovacaLogger::error("illegal base [{}] seen in the allele", bases);
    exit(EXIT_FAILURE);
}

pAllele Allele::create_allele(pBases bases, uint8_t is_ref, pMemoryPool pool)
{
    CHECK_CONDITION_EXIT(bases == nullptr, "the Allele base string cannot be nullptr");
    if (1 == bases->num) {
        return create_allele(bases->data[0], is_ref);
    }

    auto a = new ALLOC_TYPE_IN_POOL(pool, Allele) Allele();
    a->init_allele(bases, is_ref);
    return a;
}

pAllele Allele::create_allele(const char* bases_str, uint8_t is_ref, pMemoryPool pool)
{
    return create_allele(bases_str, strlen(bases_str), is_ref, pool);
}

pAllele Allele::create_allele(const char* bases_str, uint32_t num, uint8_t is_ref, pMemoryPool pool)
{
    if (1 == num) {
        return create_allele(bases_str[0], is_ref);
    }
    auto bases = new ALLOC_FLEXIBLE_IN_POOL(pool, Bases, num, uint8_t) Bases{num};
    memcpy(bases->data, bases_str, num * sizeof(uint8_t));
    bases->data[num] = 0;
    return create_allele(bases, is_ref, pool);
}

bool Allele::is_breakpoint() const { return would_be_breakpoint_allele(_bases); }

bool Allele::is_single_breakend() const { return would_be_single_breakend_allele(_bases); }

bool Allele::is_non_ref_allele() const
{
    return equals(*s_static_obj->_non_ref_allele) || equals(*s_static_obj->_unspecified_alternate_allele);
}

bool Allele::equals(const Allele& other, bool ignore_ref_state) const
{
    if (this == &other) {
        return true;
    }
    else if ((!ignore_ref_state && _is_ref != other._is_ref) || _bases->num != other._bases->num) {
        return false;
    }
    else {
        return _is_no_call == other._is_no_call && _is_symbolic == other._is_symbolic &&
               !strncmp((const char*)_bases->data, (const char*)other._bases->data, _bases->num);
    }
}

pBases Allele::get_bases() const { return _is_symbolic ? s_constants->k_empty_allele_bases : _bases; }

size_t Allele::hash() const
{
    size_t hash_value = _is_ref;
    for (uint32_t i = 0; i < _bases->num; ++i) {
        hash_value += (i + 1) * _bases->data[i];
    }
    return hash_value;
}

void Allele::init_allele(pBases bases, uint8_t is_ref)
{
    CHECK_CONDITION_EXIT(would_be_null_allele(bases), "null alleles are not supported");

    if (would_be_no_call_allele(bases)) {
        _bases = s_constants->k_empty_allele_bases;
        _is_no_call = true;
        CHECK_CONDITION_EXIT(is_ref, "cannot tag a no_call allele as the reference allele");
        return;
    }

    if (would_be_symbolic_allele(bases)) {
        _is_symbolic = true;
        CHECK_CONDITION_EXIT(is_ref, "cannot tag a symbolic allele as the reference allele");
    }
    else {
        for (uint32_t i = 0; i < bases->num; ++i) {
            bases->data[i] = toupper(bases->data[i]);
        }
    }

    _is_ref = is_ref;
    _bases = bases;
    _bases->data[_bases->num] = '\0';

    CHECK_CONDITION_EXIT(!acceptable_allele_bases(bases, is_ref), "unexpected base in allele bases {}", (char*)bases->data);
}

bool Allele::would_be_star_allele(pBases bases) { return 1 == bases->num && bases->data[0] == s_constants->k_spanning_deletion_allele; }

bool Allele::would_be_null_allele(pBases bases)
{
    return (1 == bases->num && bases->data[0] == s_constants->k_null_allele) || 0 == bases->num;
}

bool Allele::would_be_no_call_allele(pBases bases) { return 1 == bases->num && s_constants->k_no_call_allele == bases->data[0]; }

bool Allele::would_be_symbolic_allele(pBases bases)
{
    if (bases->num <= 1) {
        return false;
    }
    return s_constants->k_symbolic_allele_start == bases->data[0] || s_constants->k_symbolic_allele_end == bases->data[bases->num - 1] ||
           would_be_breakpoint_allele(bases) || would_be_single_breakend_allele(bases);
}

bool Allele::would_be_breakpoint_allele(pBases bases)
{
    if (bases->num <= 1) {
        return false;
    }
    for (uint32_t i = 0; i < bases->num; ++i) {
        if (s_constants->k_breakend_extending_left == bases->data[i] || s_constants->k_breakend_extending_right == bases->data[i]) {
            return true;
        }
    }
    return false;
}

bool Allele::would_be_single_breakend_allele(pBases bases)
{
    if (bases->num <= 1) {
        return false;
    }
    return s_constants->k_single_breakend_indicator == bases->data[0] ||
           s_constants->k_single_breakend_indicator == bases->data[bases->num - 1];
}

bool Allele::acceptable_allele_bases(pBases bases, bool is_ref)
{
    if (would_be_null_allele(bases)) {
        return false;
    }
    if (would_be_no_call_allele(bases) || would_be_symbolic_allele(bases)) {
        return true;
    }
    if (would_be_star_allele(bases)) {
        return !is_ref;
    }
    for (uint32_t i = 0; i < bases->num; ++i) {
        if (strchr(ALLELE_STR, bases->data[i]) == nullptr) {
            return false;
        }
    }
    return true;
}

size_t AlleleHash::operator()(pAllele const& a) const { return a->hash(); }

bool AlleleEqual::operator()(pAllele const& l, pAllele const& r) const
{
    if (l->is_reference() == r->is_reference() && l->is_called() == r->is_called() && l->is_symbolic() == r->is_symbolic() &&
        l->length() == r->length()) {
        pBases lb = l->get_display_string();
        pBases rb = r->get_display_string();
        return !strncmp((const char*)lb->data, (const char*)rb->data, lb->num);
    }
    return false;
}

bool less_bases(pBases l_bases, pBases r_bases)
{
    if (l_bases->num != r_bases->num)
        return l_bases->num < r_bases->num;
    else {
        for (uint32_t i = 0; i < l_bases->num; ++i) {
            if (l_bases->data[i] != r_bases->data[i]) {
                return l_bases->data[i] < r_bases->data[i];
            }
        }
    }
    return false;
}

}  // namespace rovaca