#include "variant.h"

#include <algorithm>

#include "allele.h"
#include "genotype_likelihoods.h"
#include "genotypes_context.hpp"
#include "rovaca_logger.h"

namespace rovaca
{

void Variant::set_start(int64_t start)
{
    _start = start;
    if (_start != INVALID_INT && _stop != INVALID_INT && !_alleles.empty()) {
        validate_variant();
    }
}

void Variant::set_stop(int64_t stop)
{
    _stop = stop;
    if (_start != INVALID_INT && _stop != INVALID_INT && !_alleles.empty()) {
        validate_variant();
    }
}

pVariant Variant::create_dbsnp_source(bcf1_t* b, pMemoryPool pool)
{
    AlleleVector alleles{pool};
    alleles.reserve(b->n_allele);
    for (uint32_t i{0}; i < b->n_allele; ++i) {
        alleles.push_back(Allele::create_allele(b->d.allele[i], 0 == i, pool));
    }

    uint32_t id_len = strlen(b->d.id);
    pBases id = new ALLOC_FLEXIBLE_IN_POOL(pool, Bases, id_len, uint8_t) Bases{id_len};
    id->data[id_len] = 0;
    id->num = id_len;
    memcpy(id->data, b->d.id, (id_len) * sizeof(uint8_t));

    pVariant v = Variant::create(pool);
    v->set_tid(b->rid);
    v->set_start(b->pos + 1);
    v->set_stop(b->pos + 1 + strlen(b->d.allele[0]) - 1);
    v->set_alleles(alleles);
    v->set_id(id);

    return v;
}

VariantType Variant::type(bool ignore_non_ref) const
{
    pVariant this_vc = const_cast<Variant*>(this);
    if (ignore_non_ref) {
        if (VT_UNINITIALIZED == _type_ignoring_non_ref) {
            this_vc->_type_ignoring_non_ref = this_vc->determine_type(ignore_non_ref);
        }
        return _type_ignoring_non_ref;
    }
    else {
        if (VT_UNINITIALIZED == _type) {
            this_vc->_type = this_vc->determine_type(ignore_non_ref);
        }
        return _type;
    }
}

const AlleleVector& Variant::alt_alleles() const
{
    if (_alt_alleles.empty()) {
        Variant* vc = const_cast<Variant*>(this);
        vc->_alt_alleles.reserve(_alleles.size() - 1);
        std::copy(_alleles.begin() + 1, _alleles.end(), std::back_inserter(vc->_alt_alleles));
    }
    return _alt_alleles;
}

size_t Variant::sample_count() const { return _genotypes->size(); }

int32_t Variant::get_max_ploidy(int32_t default_ploidy) { return _genotypes->get_max_ploidy(default_ploidy); }

int32_t Variant::get_called_chr_count() const
{
    int32_t count = 0;
    if (has_genotypes()) {
        pGenotype g;
        for (size_t i = 0, len = _genotypes->size(); i < len; ++i) {
            g = _genotypes->at(i);
            const AlleleVector& alleles_in_g = g->alleles();
            for (pAllele a : alleles_in_g) {
                count += a->is_called() ? 1 : 0;
            }
        }
    }
    return count;
}

int32_t Variant::get_called_chr_count(pAllele a) const
{
    int32_t count = 0;
    if (has_genotypes()) {
        for (size_t i = 0, len = _genotypes->size(); i < len; ++i) {
            count += _genotypes->at(i)->count_allele(a);
        }
    }
    return count;
}

int32_t Variant::get_allele_index(pAllele a) const
{
    int32_t index = INVALID_INT;
    for (int32_t i = 0, len = (int32_t)_alleles.size(); i < len; ++i) {
        if (_alleles.at(i)->equals(*a)) {
            index = i;
            break;
        }
    }
    return index;
}

Int32Vector Variant::get_gl_indices_of_alternate_allele(pAllele a, pMemoryPool pool) const
{
    int32_t index = get_allele_index(a);
    CHECK_CONDITION_EXIT(INVALID_INT == index, "Allele not in this Variant");
    return GenotypeLikelihoods::get_plindices_of_alleles(0, index, pool);
}

void Variant::set_alleles(const AlleleVector& alleles)
{
    _alleles.clear();
    _alleles.reserve(alleles.size());
    std::copy(alleles.begin(), alleles.end(), std::back_inserter(_alleles));
    _ref = _alleles.front();
    _biallelic_alt = 2 == _alleles.size() ? _alleles.at(1) : nullptr;
    if (_start != INVALID_INT && _stop != INVALID_INT && !_alleles.empty()) {
        validate_variant();
    }
}

bool Variant::is_simple_indel() const
{
    pAllele first_alt = alternate_allele_at(0);
    return VT_INDEL == type() && is_biallelic() && _ref->length() > 0 && first_alt->length() > 0 &&
           _ref->get_bases()->data[0] == first_alt->get_bases()->data[0] && (_ref->length() == 1 || first_alt->length() == 1);
}

bool Variant::is_simple_insertion() const { return is_simple_indel() && 1 == _ref->length(); }

bool Variant::is_simple_deletion() const { return is_simple_indel() && 1 == alternate_allele_at(0)->length(); }

void Variant::add_non_ref_symbolic_allele()
{
    _alt_alleles.clear();
    _type = _type_ignoring_non_ref = VT_UNINITIALIZED;
    _alleles.emplace_back(StaticAllele::get_instance()->_non_ref_allele.get());
    _biallelic_alt = 2 == _alleles.size() ? _alleles.at(1) : nullptr;
}

bool Variant::has_symbolic_alleles(const AlleleVector& alleles)
{
    return std::any_of(alleles.begin(), alleles.end(), [](pAllele a) { return a->is_symbolic(); });
}

bool Variant::has_allele(pAllele a) const { return has_allele(a, false, true); }

bool Variant::has_allele(pAllele a, bool ignore_ref_state, bool consider_ref_allele) const
{
    if ((consider_ref_allele && a->equals(*_ref)) || (_biallelic_alt && a->equals(*_biallelic_alt))) {
        return true;
    }
    const AlleleVector& as = consider_ref_allele ? alleles() : alt_alleles();
    return std::any_of(as.begin(), as.end(), [&](pAllele aa) { return aa->equals(*a, ignore_ref_state); });
}

bool Variant::equals(const Variant& v) const
{
    return this->get_start() == v.get_start() && this->ref_allele()->equals(*v.ref_allele()) &&
           this->biallelic_alt()->equals(*v.biallelic_alt());
}

VariantType Variant::determine_type(bool ignore_non_ref)
{
    size_t allele_num = _alleles.size();
    if (ROVACA_LIKELY(allele_num >= 1)) {
        return determine_polymorphic_type(ignore_non_ref);
    }
    else if (allele_num == 1) {
        return VT_UNINITIALIZED;
    }
    RovacaLogger::error("unexpected error: requested type of variant_context with no alleles");
    exit(EXIT_FAILURE);
}

VariantType Variant::determine_polymorphic_type(bool ignore_non_ref)
{
    VariantType type = VT_UNINITIALIZED;
    bool non_ref_allele_found = false;

    for (const auto& a : _alleles) {
        if (a == ref_allele()) {
            continue;
        }
        if (ignore_non_ref && a->is_non_ref_allele()) {
            non_ref_allele_found = true;
            continue;
        }
        VariantType biallelic_type = type_of_biallelic_variant(ref_allele(), a);
        if (VT_UNINITIALIZED == type) {
            type = biallelic_type;
        }
        else if (biallelic_type != type) {
            return VT_MIXED;
        }
    }

    if (VT_UNINITIALIZED == type && non_ref_allele_found) {
        return VT_NO_VARIATION;
    }

    return type;
}

VariantType Variant::type_of_biallelic_variant(pAllele ref, pAllele a)
{
    CHECK_CONDITION_EXIT(ref->is_symbolic(), "unexpected error: encountered a record with a symbolic reference allele");

    if (a->is_symbolic()) {
        return VT_SYMBOLIC;
    }

    if (ref->length() == a->length()) {
        if (a->length() == 1) {
            return VT_SNP;
        }
        else {
            return VT_MNP;
        }
    }
    return VT_INDEL;
}

void Variant::validate_stop()
{
    int64_t length = _stop - _start + 1;
    CHECK_CONDITION_EXIT(!has_symbolic_alleles() && length != ref_allele()->length(),
                         "bug: genome_loc {}:{}-{} has a size == {} but the variation reference allele has length {}", _tid, _start, _stop,
                         length, ref_allele()->length());
}

void Variant::validate_alleles()
{
    bool already_seen_ref = false;
    for (const auto& a : _alleles) {
        if (a->is_reference()) {
            CHECK_CONDITION_EXIT(already_seen_ref, "bug: received two reference");
            already_seen_ref = true;
        }
        CHECK_CONDITION_EXIT(!a->is_called(), "bug: cannot add a no call allele to a variant");
    }
    CHECK_CONDITION_EXIT(!already_seen_ref, "no reference allele found in variant");
    CHECK_CONDITION_EXIT(!_ref->is_reference(), "_ref is not reference");
}

size_t VariantHash::operator()(pVariant const& v) const
{
    size_t hash_value = 31 * v->get_start() + v->get_stop();
    const AlleleVector& alleles = v->alleles();
    for (pAllele a : alleles) {
        hash_value += a->hash();
    }
    return hash_value;
}

bool VariantEqual::operator()(pVariant const& l, pVariant const& r) const
{
    if (l == r) {
        return true;
    }
    if (l->get_start() != r->get_start()) {
        return false;
    }
    if (l->allele_num() != r->allele_num()) {
        return false;
    }
    return VariantHash()(l) == VariantHash()(r);
}

}  // namespace rovaca