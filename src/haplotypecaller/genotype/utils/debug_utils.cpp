#include "debug_utils.h"

#include <stdio.h>

#include <algorithm>
#include <iostream>
#include <sstream>

#include "bam_data_pool.hpp"
#include "constants_str.hpp"
#include "genotype_macors.h"
#include "genotype_struct.h"
#include "genotypes_context.hpp"
#include "haplotype.h"
#include "htslib/sam.h"
#include "indexed_sample_list.hpp"
#include "info_data.hpp"
#include "rovaca_logger.h"
#include "read_record.h"
#include "simple_interval.h"
#include "variant.h"

namespace rovaca
{

static pConstantsStr s_const = ConstantsStr::get_instance();

std::string DebugUtils::cigar2str(pCigar c)
{
    std::string ret;
    for (uint32_t i = 0; i < c->num; ++i) {
        uint32_t ele = c->data[i];
        ret.append(std::to_string(bam_cigar_oplen(ele))).append(1, BAM_CIGAR_STR[bam_cigar_op(ele)]);
    }
    return ret;
}

std::string DebugUtils::cigar2str(const uint32_t* cigar, uint32_t num)
{
    std::string ret;
    for (uint32_t i = 0; i < num; ++i) {
        uint32_t ele = cigar[i];
        ret.append(std::to_string(bam_cigar_oplen(ele))).append(1, BAM_CIGAR_STR[bam_cigar_op(ele)]);
    }
    return ret;
}

pCigar DebugUtils::str2cigar(const std::string& s, pMemoryPool pool)
{
    uint32_t cigar_num = std::count_if(s.begin(), s.end(), [](char c) { return bam_cigar_table[uint8_t(c)] != -1; });
    auto* new_ciagr = new ALLOC_FLEXIBLE_IN_POOL(pool, Cigar, cigar_num, uint32_t) Cigar{cigar_num};
    uint32_t* buffer = new_ciagr->data;
    auto num = (size_t)cigar_num;
    CHECK_CONDITION_EXIT(sam_parse_cigar(s.c_str(), nullptr, &buffer, &num) != cigar_num, "sam_parse_cigar error");
    return new_ciagr;
}

std::string DebugUtils::bases2str(pBases b)
{
    std::string ret(b->num, '0');
    for (uint32_t i = 0; i < b->num; ++i) {
        ret[i] = char(b->data[i]);
    }
    return ret;
}

std::string DebugUtils::bases2str(const uint8_t* b, int32_t len)
{
    std::string ret(len, '0');
    for (int32_t i = 0; i < len; ++i) {
        ret[i] = char(b[i]);
    }
    return ret;
}

pBases DebugUtils::str2bases(const std::string& s, pMemoryPool pool)
{
    uint32_t num = uint32_t(s.size());
    pBases b = new ALLOC_FLEXIBLE_IN_POOL(pool, Bases, num, uint8_t) Bases{num};
    for (uint32_t i = 0; i < num; ++i) {
        b->data[i] = uint8_t(s.at(i));
    }
    return b;
}

std::string DebugUtils::interval2str(pInterfaceLocatable i, bam_hdr_t* header)
{
    std::string ret(header->target_name[i->get_tid()]);
    ret.append(1, ':').append(std::to_string(i->get_start())).append(1, '-').append(std::to_string(i->get_stop()));
    return ret;
}

pSimpleInterval DebugUtils::str2interval(const std::string& str, bam_hdr_t* header, pMemoryPool pool)
{
    std::string chr = str.substr(0, str.find(':'));
    std::string pos = str.substr(str.find(':') + 1);
    int64_t start = std::stol(pos.substr(0, pos.find('-')));
    int64_t end = std::stol(pos.substr(pos.find('-') + 1));
    int32_t tid = chr_name2id(chr.c_str(), header);
    return SimpleInterval::create(tid, start, end, pool);
}

std::string DebugUtils::haplotype2str(pHaplotype h, bam_hdr_t* header)
{
    pConstantsStr s_con = ConstantsStr::get_instance();
    pCigar c = h->cigar();
    pSimpleInterval i = h->interval();
    pBases b = h->get_display_string();

    std::string ret;
    ret.append(cigar2str(c))
        .append(1, s_con->k_field_separator)
        .append(interval2str(i, header))
        .append(1, s_con->k_field_separator)
        .append(std::to_string(h->score()))
        .append(1, s_con->k_field_separator)
        .append(std::to_string(h->alignment_start_hap_wrt_ref()))
        .append(1, s_con->k_field_separator)
        .append(std::to_string(h->kmer_size()))
        .append(1, s_con->k_field_separator)
        .append(std::to_string(h->is_reference() ? 1 : 0))
        .append(1, s_con->k_field_separator)
        .append(bases2str(b))
        .append(1, '\n');
    return ret;
}

pHaplotype DebugUtils::str2haplotype(const std::string& s, bam_hdr_t* header, pMemoryPool pool)
{
    // 单倍体格式：cigar interval score ali kmer is_ref bases
    std::istringstream is(s);
    std::string cigar, interval, score, ali, kmer, is_ref, bases;
    is >> cigar >> interval >> score >> ali >> kmer >> is_ref >> bases;

    pCigar c = str2cigar(cigar, pool);
    pSimpleInterval i = str2interval(interval, header, pool);
    double core = (score == "nan") ? NAN : (score == "inf") ? INFINITY : (score == "-inf") ? double(-INFINITY) : std::stod(score);
    int64_t a = std::stol(ali);
    int32_t k = std::stoi(kmer);
    bool r = std::stoi(is_ref);
    pBases b = str2bases(bases, pool);

    pHaplotype h = Haplotype::create(pool);
    h->set_cigar(c);
    h->set_interval(i);
    h->set_score(core);
    h->set_alignment_start_hap_wrt_ref(a);
    h->set_kmer_size(k);
    h->init_haplotype(b, r);
    return h;
}

std::string DebugUtils::read2str(pReadRecord r, bam_hdr_t* header, pMemoryPool pool)
{
    pConstantsStr s_con = ConstantsStr::get_instance();

    const char* name = r->qname();
    uint16_t flag = r->flag();
    std::string interval = interval2str(r, header);
    uint8_t mq = r->mapping_quality();
    std::string cigar = cigar2str(r->cigar(), r->cigar_length());
    int32_t mtid = r->mate_tid();
    int64_t mpos = r->mate_pos();
    int64_t isize = r->insert_size();
    std::string seq = bases2str(r->decode_to_str(pool));
    std::string qual = bases2str(r->qual(), r->seq_length());
    std::for_each(qual.begin(), qual.end(), [](char& c) { c += 33; });

    std::string ret;
    ret.append(name)
        .append(1, s_con->k_field_separator)
        .append(std::to_string(flag))
        .append(1, s_con->k_field_separator)
        .append(interval)
        .append(1, s_con->k_field_separator)
        .append(std::to_string(mq))
        .append(1, s_con->k_field_separator)
        .append(cigar)
        .append(1, s_con->k_field_separator)
        .append(mtid >= 0 ? header->target_name[mtid] : std::to_string(mtid))
        .append(1, s_con->k_field_separator)
        .append(std::to_string(mpos))
        .append(1, s_con->k_field_separator)
        .append(std::to_string(isize))
        .append(1, s_con->k_field_separator)
        .append(seq)
        .append(1, s_con->k_field_separator)
        .append(qual)
        .append(1, '\n');
    return ret;
}

pReadRecord DebugUtils::str2read(const std::string& s, bam_hdr_t* header, pMemoryPool pool, pBamDataPool bpool)
{
    // Read格式：name flag interval mq cigar mtid mpos isize seq qual
    std::string name, flag, interval, mq, cigar, mtid, mpos, isize, seq, qual;
    std::istringstream is{s};
    is >> name >> flag >> interval >> mq >> cigar >> mtid >> mpos >> isize >> seq >> qual;

    std::for_each(qual.begin(), qual.end(), [](char& c) { c -= 33; });
    pSimpleInterval i = str2interval(interval, header, pool);
    pCigar c = str2cigar(cigar, pool);
    int32_t mtidn = chr_name2id(mtid.c_str(), header);

    bam1_t* b = bpool->alloc_bam_struct_pool();
    int ret = bam_set1(b, name.size(), name.c_str(), std::stoi(flag), i->get_tid(), i->get_start(), std::stoi(mq), c->num, c->data, mtidn,
                       std::stol(mpos), std::stol(isize), seq.size(), seq.c_str(), qual.c_str(), 0);

    CHECK_CONDITION_EXIT(ret < 0, "bam_set1 error");
    bpool->peek_bam(b);
    return ReadRecord::create(pool, header, b);
}

std::string DebugUtils::variant2str(pVariant v, bam_hdr_t* header, bool is_vcf_model)
{
    std::string ret;
    ret.reserve(300);

    // chr pos id ref
    ret.append(chr_id2name(v->get_tid(), header))
        .append(1, s_const->k_field_separator)
        .append(std::to_string(v->get_start()))
        .append(1, s_const->k_field_separator);
    if (ROVACA_LIKELY(v->db_id() == nullptr)) {
        ret.append(1, s_const->k_empty_id_field);
    }
    else {
        ret.append(bases2str(v->db_id()));
    }
    ret.append(1, s_const->k_field_separator);
    ret.append(bases2str(v->ref_allele()->get_display_string())).append(1, s_const->k_field_separator);

    // alt
    if (v->is_variant()) {
        const AlleleVector& alt_alleles = v->alt_alleles();
        ret.append(bases2str(alt_alleles.at(0)->get_display_string()));
        for (size_t i = 1, alt_alleles_num = alt_alleles.size(); i < alt_alleles_num; ++i) {
            ret.append(1, s_const->k_info_field_array_separator).append(bases2str(alt_alleles.at(i)->get_display_string()));
        }
    }
    else {
        ret.append(1, s_const->k_empty_alternate_allele_field);
    }
    ret.append(1, s_const->k_field_separator);

    // qual
    ret = v->has_log10_error() ? format_qual_value(v->get_phred_scaled_qual(), ret) : ret.append(1, s_const->k_missing_value_v4);
    ret.append(1, s_const->k_field_separator);

    // 当前永远不会有 Filter，暂未实现
    ret.append(1, s_const->k_unfiltered).append(1, s_const->k_field_separator);

    // info
    pInfoData info = v->info();
    if (!is_vcf_model && !info) {
        // gvcf模式，非变异位点，info信息只有 end_key,存储在 Variant 中
        ret.append(std::to_string(v->end_key())).append(1, s_const->k_field_separator);
    }
    else {
        info->to_string(is_vcf_model, ret);
        ret.back() = s_const->k_field_separator;  // 最后一个 ; -> \t
    }

    // format
    pGenotypesContext gc = v->genotype();
    if (gc) {
        pGenotype g = !gc->empty() ? gc->at(0) : nullptr;
        if (g) {
            std::pmr::map<pAllele, char> mm = build_allele_map(v->alleles());
            g->append_genotype_data(ret, mm);
        }
    }
    ret.back() = '\n';

    return ret;
}

pVariant DebugUtils::str2variant([[maybe_unused]] const std::string& s, [[maybe_unused]] pMemoryPool pool) { return nullptr; }

void DebugUtils::format_double(double d, const char* format, std::string& s)
{
    char str[32]{};
    sprintf(str, format, d);
    s.append(str, strlen(str));
}

std::string& DebugUtils::format_qual_value(double d, std::string& str)
{
    format_double(d, "%.2f", str);
    if (strcmp(str.c_str() + str.size() - 3, ".00") == 0) {
        str.resize(str.size() - 3);
    }
    return str;
}

void DebugUtils::format_vcf_double(double d, std::string& s)
{
    const char* format;
    if (d < 1) {
        if (d < 0.01) {
            if (std::abs(d) >= 1e-20) {
                format = "%.3e";
            }
            else {
                s.append("0.00");  // return a zero format
                return;
            }
        }
        else {
            format = "%.3f";
        }
    }
    else {
        format = "%.2f";
    }
    format_double(d, format, s);
}

std::string& DebugUtils::write_int_vector(const Int32Vector& data, std::string& s)
{
    if (data.empty()) {
        return s;
    }

    s.append(std::to_string(data.front()));
    for (size_t i = 1, len = data.size(); i < len; ++i) {
        s.append(1, s_const->k_info_field_array_separator).append(std::to_string(data.at(i)));
    }

    return s;
}

std::string& DebugUtils::write_long_vector(const LongVector& data, std::string& s)
{
    if (data.empty()) {
        return s;
    }

    s.append(std::to_string(data.front()));
    for (size_t i = 1, len = data.size(); i < len; ++i) {
        s.append(1, s_const->k_info_field_array_separator).append(std::to_string(data.at(i)));
    }

    return s;
}

std::string& DebugUtils::write_double_vector(const DoubleVector& data, std::string& s)
{
    if (data.empty()) {
        return s;
    }

    format_vcf_double(data.front(), s);
    for (size_t i = 1, len = data.size(); i < len; ++i) {
        s.append(1, s_const->k_info_field_array_separator);
        format_vcf_double(data.at(i), s);
    }

    return s;
}

std::pmr::map<pAllele, char> DebugUtils::build_allele_map(const AlleleVector& alleles)
{
    std::pmr::map<pAllele, char> mm{alleles.get_allocator()};
    mm.insert({StaticAllele::get_instance()->_no_call.get(), s_const->k_empty_allele});
    for (size_t i = 0, len = alleles.size(); i < len; ++i) {
        mm.insert({alleles.at(i), '0' + i});
    }
    return mm;
}

std::pmr::map<pAllele, int32_t> DebugUtils::bcf_build_allele_map(const AlleleVector& alleles)
{
    std::pmr::map<pAllele, int32_t> mm{alleles.get_allocator()};
    for (int32_t i = 0, len = int32_t(alleles.size()); i < len; ++i) {
        mm.insert({alleles.at(i), i});
    }
    return mm;
}

const char* DebugUtils::chr_id2name(int32_t id, bam_hdr_t* header) { return header->target_name[id]; }

int32_t DebugUtils::chr_name2id(const char* name, bam_hdr_t* header)
{
    if (strcmp(name, "null") == 0) {
        return -1;
    }

    for (int32_t i = 0; i < header->n_targets; ++i) {
        if (strcmp(header->target_name[i], name) == 0) {
            return i;
        }
    }
    RovacaLogger::error("tid == -1");
    exit(EXIT_FAILURE);
}

void DebugUtils::print_reads(const ReadVector& rs, bam_hdr_t* hdr, pMemoryPool pool)
{
    setvbuf(stdout, nullptr, _IONBF, 0);
    for (pReadRecord r : rs) {
        std::string s = read2str(r, hdr, pool);
        fprintf(stdout, "%s", s.c_str());
    }
    fflush(stdout);
    printf("\n");
}

void DebugUtils::print_reads(const ReadHashSet& rs, bam_hdr_t* hdr, pMemoryPool pool)
{
    print_reads(ReadVector{{rs.begin(), rs.end()}, pool}, hdr, pool);
}

const char* debug_file = "/data/rovaca-dev/zhangtiefeng/script/python/result/rovaca.log";
void DebugUtils::print_reads_to_file(const ReadVector& rs, bam_hdr_t* hdr, pMemoryPool pool)
{
    FILE* file = fopen(debug_file, "w");
    if (file == nullptr) {
        printf("open debug_file error: %s\n", strerror(errno));
        return;
    }

    for (pReadRecord r : rs) {
        std::string s = read2str(r, hdr, pool);
        fprintf(file, "%s", s.c_str());
    }

    fclose(file);
    printf("内容已成功写入文件。\n");
}

void DebugUtils::print_reads(const std::vector<bam1_t*>& rs, bam_hdr_t* hdr, pMemoryPool pool)
{
    ReadVector rr{pool};
    for (bam1_t* b : rs) {
        b->core.pos += 1;
        b->core.mpos += 1;
        rr.push_back(ReadRecord::create(pool, hdr, b));
    }
    print_reads(rr, hdr, pool);
}

void DebugUtils::print_reads_to_file(const ReadHashSet& rs, bam_hdr_t* hdr, pMemoryPool pool)
{
    print_reads_to_file(ReadVector{{rs.begin(), rs.end()}, pool}, hdr, pool);
}

void DebugUtils::print_reads_to_file(const std::vector<bam1_t*>& rs, bam_hdr_t* hdr, pMemoryPool pool)
{
    ReadVector rr{pool};
    for (bam1_t* b : rs) {
        b->core.pos += 1;
        b->core.mpos += 1;
        rr.push_back(ReadRecord::create(pool, hdr, b));
    }
    print_reads_to_file(rr, hdr, pool);
}

void DebugUtils::print_haplotypes(const HaplotypeVector& hs, bam_hdr_t* hdr)
{
    setvbuf(stdout, nullptr, _IONBF, 0);
    for (pHaplotype h : hs) {
        std::string s = haplotype2str(h, hdr);
        fprintf(stdout, "%s", s.c_str());
    }
    fflush(stdout);
    printf("\n");
}

void DebugUtils::print_haplotypes_to_file(const HaplotypeVector& hs, bam_hdr_t* hdr)
{
    FILE* file = fopen(debug_file, "w");
    if (file == nullptr) {
        printf("open debug_file error: %s\n", strerror(errno));
        return;
    }

    for (pHaplotype h : hs) {
        std::string s = haplotype2str(h, hdr);
        fprintf(file, "%s", s.c_str());
    }

    fclose(file);
    printf("内容已成功写入文件。\n");
}

}  // namespace rovaca