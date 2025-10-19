#include "writer.h"

#include <chrono>
#include <cstddef>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <exception>
#include <fstream>
#include <map>
#include <sstream>
#include <string>
#include <unordered_map>

#include "common/utils/BlockingQueue.h"
#include "genotype/genotype_enum.h"
#include "genotype/genotype_macors.h"
#include "htslib/bgzf.h"
#include "htslib/hfile.h"
#include "htslib/hts.h"
#include "htslib/kstring.h"
#include "htslib/sam.h"
#include "htslib/vcf.h"
#include "rovaca_logger.h"

namespace rovaca
{

#define BCF_CACHE_SIZE      (128 * 1024)
#define ROVACA_HTS_CACHE_SIZE (64 * 1024 * 1024)
#define RESULT_CACHE_COUNT  (128 * 1024 * 1024 / 8)

static inline bool is_non_variant(kstring_t *s)
{
    int32_t alt_pos = 4;
    const char *a = s->s;
    while (alt_pos--) {
        for (; *a != '\t'; ++a) {}
        ++a;
    }
    return *a == '<';
}

static inline size_t number_of_digits(int64_t n)
{
    size_t count = 0;
    do {
        n = n / 10;
        ++count;
    } while (n != 0);
    return count;
}

Writer::Writer(bool index, int32_t compression_level, pHCArgs args, bam_hdr_t *bam_hdr, htsThreadPool *tp,
               BlockingQueue<pWriterTask> *result_queue)
    : index_(index)
    , compression_level_(compression_level)
    , args_{args}
    , bam_hdr_{bam_hdr}
    , hts_pool_{tp}
    , result_queue_{result_queue}
    , out_{}
    , last_id_back_{INVALID_INT}
    , merge_thread_{nullptr}
    , ks_pool_{}
{}

kstring_t *Writer::pop_bcf_cache()
{
    kstring_t *ret = ks_pool_.pop();
    if (nullptr == ret) {
        ret = new kstring_t{0, 0, nullptr};
        int32_t r = ks_resize(ret, BCF_CACHE_SIZE);
        CHECK_CONDITION_EXIT(r != 0, "error: ks_resize");
    }
    return ret;
}

void Writer::push_bcf_cache(kstring_t *s)
{
    ks_clear(s);
    ks_pool_.push(s);
}

bool Writer::start()
{
    if (!args_->output) {
        return false;
    }

    // 创建写出文件
    const char *mode = "w";
    if (is_gz_file(args_->output)) {
        const char *level_to_model[10]{"wz0", "wz1", "wz2", "wz3", "wz4", "wz5", "wz", "wz7", "wz8", "wz9"};
        mode = level_to_model[compression_level_];
        if (index_) {
            args_->idx_name = std::string(args_->output) + ".tbi";
        }
    }
    else {
        index_ = false;  // 非 gz 不建立索引
    }
    out_.out_file = hts_open(args_->output, mode);
    if (!out_.out_file) {
        RovacaLogger::error("failed to open: {}", args_->output);
        exit(EXIT_FAILURE);
    }

    // 创建 初始化 header内容
    if (args_->reference_confidence_mode == ReferenceConfidenceMode::GVCF) {
        out_.hdr = init_gvcf_header(args_, bam_hdr_);
    }
    else if (args_->reference_confidence_mode == ReferenceConfidenceMode::NONE) {
        out_.hdr = init_vcf_header(args_, bam_hdr_);
    }
    else if (args_->reference_confidence_mode == ReferenceConfidenceMode::BP_RESOLUTION) {
        RovacaLogger::error("BP_RESOLUTION is not currently supported");
        exit(EXIT_FAILURE);
    }
    else {
        RovacaLogger::error("invalid ReferenceConfidenceMode/ERC type");
        exit(EXIT_FAILURE);
    }
    CHECK_CONDITION_EXIT(out_.hdr == nullptr, "can not create vcf header");

    merge_thread_ = new std::thread{&Writer::merge_task_thread_call, this};
    pthread_setname_np(merge_thread_->native_handle(), "Writer");

    return true;
}

void Writer::stop()
{
    merge_thread_->join();
    delete merge_thread_;

    kstring_t *s;
    while ((s = ks_pool_.pop()) != nullptr) {
        ks_free(s);
        delete s;
    }

    RovacaLogger::info("writer done");
}

void Writer::merge_task_thread_call()
{
    pWriterTask task{nullptr};
    int32_t ret, expected_id = 0;

    hts_set_thread_pool(out_.out_file, hts_pool_);

    // 存储del后面特殊位点的cache
    ks_resize(&del_cache_, KS_CACHE_SIZE);

    ret = bcf_hdr_write(out_.out_file, out_.hdr);
    CHECK_CONDITION_EXIT(ret != 0, "write header: {}", std::strerror(errno));

    if (!args_->idx_name.empty()) {
        ret = bcf_idx_init(out_.out_file, out_.hdr, 0, args_->idx_name.c_str());
        CHECK_CONDITION_EXIT(ret != 0, "failed to initialise index for current output");
    }

    std::vector<pWriterTask> discontinuous_result{RESULT_CACHE_COUNT, nullptr};

    while (expected_id != last_id_back_) {
        task = result_queue_->pop(1);
        if (nullptr == task) {
            continue;
        }

        if (task->source_id == expected_id) {
            if (nullptr != task->cache) {
                if (task->cache->l > 0) {
                    resolve_overrides_and_write_task(task);
                }
                push_bcf_cache(task->cache);
            }
            delete task;
            ++expected_id;
        }
        else {
            discontinuous_result[task->source_id] = task;
        }

        while ((task = discontinuous_result[expected_id]) != nullptr) {
            if (nullptr != task->cache) {
                if (task->cache->l > 0) {
                    resolve_overrides_and_write_task(task);
                }
                push_bcf_cache(task->cache);
            }
            delete task;
            expected_id++;
        }
    }

    ks_free(&del_cache_);

    if (!args_->idx_name.empty()) {
        ret = bcf_idx_save(out_.out_file);
        CHECK_CONDITION_EXIT(ret != 0, "failed to build index");
    }
    hts_close(out_.out_file);
    bcf_hdr_destroy(out_.hdr);
}

void Writer::close_file()
{
    if (!args_->idx_name.empty()) {
        int ret = bcf_idx_save(out_.out_file);
        CHECK_CONDITION_EXIT(ret != 0, "failed to build index");
    }
    hts_close(out_.out_file);
}

bool Writer::check_deletion_variant(const std::shared_ptr<char> &ref, kstring_t *s, size_t *current_offset)
{
    char *c, *p, *i, *r, *q, *e, *o;
    bool punctured{false};
    int32_t tid;
    int64_t start, stop;
    kstring_t sub_str{0, 0, nullptr};
    size_t offset{0}, current_line_len, start_len, stop_len;
    for (offset = 0; offset < s->l; offset += current_line_len) {
        current_line_len = 0;
        for (q = s->s + offset; *q != '\n'; ++q, ++current_line_len) {}
        current_line_len += 1;

        sub_str = {current_line_len, current_line_len, s->s + offset};

        // 变异位点，直接从此位置输出
        if (!is_non_variant(&sub_str)) {
            *current_offset = offset;
            break;
        }

        for (c = sub_str.s, q = c; *q != '\t'; ++q) {}  // chr
        strncpy(contig_cache, c, q - c);
        contig_cache[q - c] = '\0';  // 追加结尾
        tid = bcf_hdr_name2id(out_.hdr, contig_cache);

        for (p = q + 1, q = p; *q != '\t'; ++q) {}  // pos

        start_len = q - p;
        strncpy(pos_cache, p, start_len);
        pos_cache[start_len] = '\0';  // 追加结尾
        start = atol(pos_cache);

        // start < next_available_start_成立时要看 end_key 是否大于 next_available_start_
        if (start < next_available_start_) {
            // 走到这里的必是非变异，stop取end_key
            for (i = q + 1, q = i; *q != '\t'; ++q) {}  // id
            for (r = q + 1, q = r; *q != '\t'; ++q) {}  // ref
            for (o = q + 1, q = o; *q != '\t'; ++q) {}  // alt
            for (o = q + 1, q = o; *q != '\t'; ++q) {}  // filter
            for (o = q + 1, q = o; *q != '\t'; ++q) {}  // qual
            for (e = q + 1, q = e; *q != '='; ++q) {}   // =
            for (e = q + 1, q = e; *q != '\t'; ++q) {}  // end
            stop_len = q - e;
            strncpy(pos_cache, e, stop_len);
            pos_cache[stop_len] = '\0';  // 追加结尾
            stop = atol(pos_cache);
            if (stop < next_available_start_) {
                continue;  // 被del跨越
            }

            /*
            以下逻辑需要适配 pos ref 到 next_available_start_
            当 pos=999, next_available_start_为1000时, 多一个字符, 为了避免对原字符串更改, 此处需要特殊处理
            */

            r[0] = ref.get()[next_available_start_ - 1];
            if (start_len == number_of_digits(next_available_start_)) {
                memset(del_cache_.s, 0, (start_len + 8) & ~7);
                sprintf(del_cache_.s, "%ld", next_available_start_);
                strncpy(p, del_cache_.s, start_len);
                *current_offset = offset;
                break;
            }
            else {
                del_cache_.l = current_line_len + number_of_digits(next_available_start_) - start_len;
                size_t current_len = 0;
                size_t len = 0;
                for (char *cc = c; *cc != '\t'; ++cc, ++len) {}
                strncpy(del_cache_.s + current_len, c, len + 1);
                current_len = len + 1;
                sprintf(del_cache_.s + current_len, "%ld\t", next_available_start_);
                current_len += number_of_digits(next_available_start_) + 1;
                len = 0;
                for (char *ii = i; *ii != '\n'; ++ii, ++len) {}
                len += 1;
                strncpy(del_cache_.s + current_len, i, len);
                current_len += len;
                CHECK_CONDITION_EXIT(current_len != del_cache_.l, "resolve del overrides error");
                if (out_.out_file->format.compression == bgzf) {
                    int32_t ret = int32_t(bgzf_write(out_.out_file->fp.bgzf, del_cache_.s, del_cache_.l));
                    CHECK_CONDITION_EXIT(size_t(ret) != del_cache_.l, "error: vcf_write_line");

                    if (index_) {
                        int32_t btid = hts_idx_tbi_name(out_.out_file->idx, tid, contig_cache);
                        CHECK_CONDITION_EXIT(btid == -1, "error: hts_idx_tbi_name");

                        ret = bgzf_idx_push(out_.out_file->fp.bgzf, out_.out_file->idx, btid, start - 1, stop - 1,
                                            bgzf_tell(out_.out_file->fp.bgzf), 1);
                        CHECK_CONDITION_EXIT(ret != 0, "error: hts_idx_push");
                    }
                }
                else {
                    int32_t ret = int32_t(hwrite(out_.out_file->fp.hfile, del_cache_.s, del_cache_.l));
                    CHECK_CONDITION_EXIT(size_t(ret) != del_cache_.l, "error: vcf_write_line");
                }
                *current_offset = offset + current_line_len;
                break;
            }
        }
        else {
            *current_offset = offset;
            break;
        }
    }

    if (offset == s->l && 0 == *current_offset) {
        punctured = true;
    }

    return punctured;
}

void Writer::resolve_overrides_and_write_task(pWriterTask task)
{
    const char *q, *w;
    int32_t tid, btid, ret = 0;
    int64_t start, stop;
    bool punctured{false};
    size_t current_line_len, current_offset = 0;

    kstring_t *str = task->cache;
    kstring_t sub_str{0, 0, nullptr};

    bool gvcf_model = args_->emit_reference_confidence();

    if (gvcf_model) {
        if (pre_tid_ == task->tid) {
            if (pre_is_del_) {
                punctured = check_deletion_variant(task->ref, str, &current_offset);
            }

            pre_is_del_ = punctured ? pre_is_del_ : task->del_line_offset;
            next_available_start_ = punctured ? next_available_start_ : task->next_available_start;
        }
        else {
            pre_tid_ = task->tid;
            pre_is_del_ = task->del_line_offset;
            next_available_start_ = task->next_available_start;
        }
    }

    // gvcf 模式击穿不输出
    if (ROVACA_UNLIKELY(punctured)) {
        return;
    }

    if (out_.out_file->format.compression == bgzf) {
        while (current_offset < str->l) {
            current_line_len = 0;
            for (q = str->s + current_offset; *q != '\n'; ++q, ++current_line_len) {}
            current_line_len += 1;

            sub_str = {current_line_len, current_line_len, str->s + current_offset};

            ret = int32_t(bgzf_write(out_.out_file->fp.bgzf, sub_str.s, sub_str.l));
            CHECK_CONDITION_EXIT(size_t(ret) != sub_str.l, "error: vcf_write_line");

            if (index_) {
                // chr
                for (q = sub_str.s, w = q; *w != '\t'; ++w) {}
                strncpy(contig_cache, q, w - q);
                contig_cache[w - q] = '\0';  // 追加结尾
                tid = bcf_hdr_name2id(out_.hdr, contig_cache);

                // pos
                for (q = w + 1, w = q; *w != '\t'; ++w) {}
                strncpy(pos_cache, q, w - q);
                pos_cache[w - q] = '\0';  // 追加结尾
                start = atol(pos_cache);

                // id
                for (q = w + 1, w = q; *w != '\t'; ++w) {}

                // ref
                for (q = w + 1, w = q; *w != '\t'; ++w) {}
                stop = start + int64_t(w - q);
                if (gvcf_model) {
                    // alt
                    for (q = w + 1, w = q; *w != '\t'; ++w) {}
                    // qual
                    for (q = w + 1, w = q; *w != '\t'; ++w) {}

                    // non_ref 位点
                    if (*q == '.') {
                        // filter
                        for (q = w + 1, w = q; *w != '\t'; ++w) {}
                        // info
                        for (q = w + 5, w = q; *w != '\t'; ++w) {}
                        strncpy(pos_cache, q, w - q);
                        pos_cache[w - q] = '\0';  // 追加结尾
                        stop = atol(pos_cache) + 1;
                    }
                }

                btid = hts_idx_tbi_name(out_.out_file->idx, tid, contig_cache);
                CHECK_CONDITION_EXIT(btid == -1, "error: hts_idx_tbi_name");

                ret = bgzf_idx_push(out_.out_file->fp.bgzf, out_.out_file->idx, btid, start - 1, stop - 1,
                                    bgzf_tell(out_.out_file->fp.bgzf), 1);
                
                CHECK_CONDITION_EXIT(ret != 0, "error: hts_idx_push");
            }

            current_offset += current_line_len;
        }
    }
    else {
        ret = int32_t(hwrite(out_.out_file->fp.hfile, task->cache->s + current_offset, task->cache->l - current_offset));
        CHECK_CONDITION_EXIT(size_t(ret) != task->cache->l - current_offset, "error: vcf_write_line");
    }
}

bool Writer::is_gz_file(const char *filename)
{
    size_t len = strlen(filename);
    if (len < 3) {
        return false;
    }
    const char *suffix = filename + len - 3;
    return strcmp(suffix, ".gz") == 0;
}

std::string Writer::generate_gvcf_block_line(int32_t min_gq, int32_t max_gq)
{
    std::ostringstream oss;
    oss << "##GVCFBlock" << min_gq << "-" << max_gq << "=minGQ=" << min_gq << "(inclusive),maxGQ=" << max_gq << "(exclusive)";
    return oss.str();
}

std::string Writer::generate_contig_line(const char *chr_name, uint32_t len)
{
    std::ostringstream oss;
    oss << "##contig=<ID=" << std::string(chr_name) << ",length=" << len << ">";
    return oss.str();
}

std::string Writer::generate_command_line(const char *tool_name, const char *command_line)
{
    CHECK_CONDITION_EXIT(nullptr == tool_name, "args->tool_name is nullptr");
    CHECK_CONDITION_EXIT(nullptr == command_line, "args->command_line is empty");
    std::ostringstream oss;
    oss << "##LUSHCommandLine=<ID=" << std::string(tool_name) << ",CommandLine=\"" << std::string(command_line) << "\">";
    return oss.str();
}

bcf_hdr_t *Writer::init_vcf_header(pHCArgs args, bam_hdr_t *bam_hdr)
{
    bcf_hdr_t *hdr = bcf_hdr_init("w");
    if (!hdr) {
        RovacaLogger::error("failed to init bcf header");
        exit(EXIT_FAILURE);
    }
    bcf_hdr_remove(hdr, BCF_HL_FLT, "PASS");

    int32_t ret;
    // clang-format off
    ret = bcf_hdr_append(hdr, "##FILTER=<ID=LowQual,Description=\"Low quality\">");
    CHECK_CONDITION_EXIT(ret != 0, "bcf_hdr_append: FILTER");

    ret = bcf_hdr_append(hdr, "##FORMAT=<ID=AD,Number=R,Type=Integer,Description=\"Allelic depths for the ref and alt alleles in the order listed\">");
    CHECK_CONDITION_EXIT(ret != 0, "bcf_hdr_append: FORMAT=AD");
    ret = bcf_hdr_append(hdr, "##FORMAT=<ID=DP,Number=1,Type=Integer,Description=\"Approximate read depth (reads with MQ=255 or with bad mates are filtered)\">");
    CHECK_CONDITION_EXIT(ret != 0, "bcf_hdr_append: FORMAT=DP");
    ret = bcf_hdr_append(hdr, "##FORMAT=<ID=GQ,Number=1,Type=Integer,Description=\"Genotype Quality\">");
    CHECK_CONDITION_EXIT(ret != 0, "bcf_hdr_append: FORMAT=GQ");
    ret = bcf_hdr_append(hdr, "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">");
    CHECK_CONDITION_EXIT(ret != 0, "bcf_hdr_append: FORMAT=GT");
    ret = bcf_hdr_append(hdr, "##FORMAT=<ID=PL,Number=G,Type=Integer,Description=\"Normalized, Phred-scaled likelihoods for genotypes as defined in the VCF specification\">");
    CHECK_CONDITION_EXIT(ret != 0, "bcf_hdr_append: FORMAT=PL");
    // clang-format on

    std::string command_line = generate_command_line(args->tool_name, args->command_line.c_str());
    ret = bcf_hdr_append(hdr, command_line.c_str());
    CHECK_CONDITION_EXIT(ret != 0, "bcf_hdr_append: LUSHCommandLine=%s", command_line.c_str());

    // clang-format off
    ret = bcf_hdr_append(hdr, "##INFO=<ID=AC,Number=A,Type=Integer,Description=\"Allele count in genotypes, for each ALT allele, in the same order as listed\">");
    CHECK_CONDITION_EXIT(ret != 0, "bcf_hdr_append: INFO=AC");
    ret = bcf_hdr_append(hdr, "##INFO=<ID=AF,Number=A,Type=Float,Description=\"Allele Frequency, for each ALT allele, in the same order as listed\">");
    CHECK_CONDITION_EXIT(ret != 0, "bcf_hdr_append: INFO=AF");
    ret = bcf_hdr_append(hdr, "##INFO=<ID=AN,Number=1,Type=Integer,Description=\"Total number of alleles in called genotypes\">");
    CHECK_CONDITION_EXIT(ret != 0, "bcf_hdr_append: INFO=AN");
    ret = bcf_hdr_append(hdr, "##INFO=<ID=BaseQRankSum,Number=1,Type=Float,Description=\"Z-score from Wilcoxon rank sum test of Alt Vs. Ref base qualities\">");
    CHECK_CONDITION_EXIT(ret != 0, "bcf_hdr_append: INFO=BaseQRankSum");
    ret = bcf_hdr_append(hdr, "##INFO=<ID=DP,Number=1,Type=Integer,Description=\"Approximate read depth; some reads may have been filtered\">");
    CHECK_CONDITION_EXIT(ret != 0, "bcf_hdr_append: INFO=DP");
    ret = bcf_hdr_append(hdr, "##INFO=<ID=ExcessHet,Number=1,Type=Float,Description=\"Phred-scaled p-value for exact test of excess heterozygosity\">");
    CHECK_CONDITION_EXIT(ret != 0, "bcf_hdr_append: INFO=ExcessHet");
    ret = bcf_hdr_append(hdr, "##INFO=<ID=FS,Number=1,Type=Float,Description=\"Phred-scaled p-value using Fisher's exact test to detect strand bias\">");
    CHECK_CONDITION_EXIT(ret != 0, "bcf_hdr_append: INFO=FS");
    ret = bcf_hdr_append(hdr, "##INFO=<ID=InbreedingCoeff,Number=1,Type=Float,Description=\"Inbreeding coefficient as estimated from the genotype likelihoods per-sample when compared against the Hardy-Weinberg expectation\">");
    CHECK_CONDITION_EXIT(ret != 0, "bcf_hdr_append: INFO=InbreedingCoeff");
    ret = bcf_hdr_append(hdr, "##INFO=<ID=MLEAC,Number=A,Type=Integer,Description=\"Maximum likelihood expectation (MLE) for the allele counts (not necessarily the same as the AC), for each ALT allele, in the same order as listed\">");
    CHECK_CONDITION_EXIT(ret != 0, "bcf_hdr_append: INFO=MLEAC");
    ret = bcf_hdr_append(hdr, "##INFO=<ID=MLEAF,Number=A,Type=Float,Description=\"Maximum likelihood expectation (MLE) for the allele frequency (not necessarily the same as the AF), for each ALT allele, in the same order as listed\">");
    CHECK_CONDITION_EXIT(ret != 0, "bcf_hdr_append: INFO=MLEAF");
    ret = bcf_hdr_append(hdr, "##INFO=<ID=MQ,Number=1,Type=Float,Description=\"RMS Mapping Quality\">");
    CHECK_CONDITION_EXIT(ret != 0, "bcf_hdr_append: INFO=MQ");
    ret = bcf_hdr_append(hdr, "##INFO=<ID=MQRankSum,Number=1,Type=Float,Description=\"Z-score From Wilcoxon rank sum test of Alt vs. Ref read mapping qualities\">");
    CHECK_CONDITION_EXIT(ret != 0, "bcf_hdr_append: INFO=MQRankSum");
    ret = bcf_hdr_append(hdr, "##INFO=<ID=QD,Number=1,Type=Float,Description=\"Variant Confidence/Quality by Depth\">");
    CHECK_CONDITION_EXIT(ret != 0, "bcf_hdr_append: INFO=QD");
    ret = bcf_hdr_append(hdr, "##INFO=<ID=ReadPosRankSum,Number=1,Type=Float,Description=\"Z-score from Wilcoxon rank sum test of Alt vs. Ref read position bias\">");
    CHECK_CONDITION_EXIT(ret != 0, "bcf_hdr_append: INFO=ReadPosRankSum");
    ret = bcf_hdr_append(hdr, "##INFO=<ID=SOR,Number=1,Type=Float,Description=\"Symmetric Odds Ratio of 2x2 contingency table to detect strand bias\">");
    CHECK_CONDITION_EXIT(ret != 0, "bcf_hdr_append: INFO=SOR");
    // clang-format on

    generate_samples(hdr, bam_hdr);
    generate_contig(hdr, bam_hdr);

    ret = bcf_hdr_append(hdr, "##source=HaplotypeCaller");
    CHECK_CONDITION_EXIT(ret != 0, "bcf_hdr_append: source=HaplotypeCaller");

    return hdr;
}

bcf_hdr_t *Writer::init_gvcf_header(pHCArgs args, bam_hdr_t *bam_hdr)
{
    bcf_hdr_t *hdr = bcf_hdr_init("w");
    if (!hdr) {
        RovacaLogger::error("failed to init bcf header");
        exit(EXIT_FAILURE);
    }
    bcf_hdr_remove(hdr, BCF_HL_FLT, "PASS");

    int32_t ret;
    // clang-format off
    ret = bcf_hdr_append(hdr, "##ALT=<ID=NON_REF,Description=\"Represents any possible alternative allele not already represented at this location by REF and ALT\">");
    CHECK_CONDITION_EXIT(ret != 0, "bcf_hdr_append: ALT");
    
    ret = bcf_hdr_append(hdr, "##FILTER=<ID=LowQual,Description=\"Low quality\">");
    CHECK_CONDITION_EXIT(ret != 0, "bcf_hdr_append: FILTER");

    ret = bcf_hdr_append(hdr, "##FORMAT=<ID=AD,Number=R,Type=Integer,Description=\"Allelic depths for the ref and alt alleles in the order listed\">");
    CHECK_CONDITION_EXIT(ret != 0, "bcf_hdr_append: FORMAT=AD");
    ret = bcf_hdr_append(hdr, "##FORMAT=<ID=DP,Number=1,Type=Integer,Description=\"Approximate read depth (reads with MQ=255 or with bad mates are filtered)\">");
    CHECK_CONDITION_EXIT(ret != 0, "bcf_hdr_append: FORMAT=DP");
    ret = bcf_hdr_append(hdr, "##FORMAT=<ID=GQ,Number=1,Type=Integer,Description=\"Genotype Quality\">");
    CHECK_CONDITION_EXIT(ret != 0, "bcf_hdr_append: FORMAT=GQ");
    ret = bcf_hdr_append(hdr, "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">");
    CHECK_CONDITION_EXIT(ret != 0, "bcf_hdr_append: FORMAT=GT");
    ret = bcf_hdr_append(hdr, "##FORMAT=<ID=MIN_DP,Number=1,Type=Integer,Description=\"Minimum DP observed within the GVCF block\">");
    CHECK_CONDITION_EXIT(ret != 0, "bcf_hdr_append: FORMAT=MIN_DP");
    ret = bcf_hdr_append(hdr, "##FORMAT=<ID=PGT,Number=1,Type=String,Description=\"Physical phasing haplotype information, describing how the alternate alleles are phased in relation to one another; will always be heterozygous and is not intended to describe called alleles\">");
    CHECK_CONDITION_EXIT(ret != 0, "bcf_hdr_append: FORMAT=PGT");
    ret = bcf_hdr_append(hdr, "##FORMAT=<ID=PID,Number=1,Type=String,Description=\"Physical phasing ID information, where each unique ID within a given sample (but not across samples) connects records within a phasing group\">");
    CHECK_CONDITION_EXIT(ret != 0, "bcf_hdr_append: FORMAT=PID");
    ret = bcf_hdr_append(hdr, "##FORMAT=<ID=PL,Number=G,Type=Integer,Description=\"Normalized, Phred-scaled likelihoods for genotypes as defined in the VCF specification\">");
    CHECK_CONDITION_EXIT(ret != 0, "bcf_hdr_append: FORMAT=PL");
    ret = bcf_hdr_append(hdr, "##FORMAT=<ID=PS,Number=1,Type=Integer,Description=\"Phasing set (typically the position of the first variant in the set)\">");
    CHECK_CONDITION_EXIT(ret != 0, "bcf_hdr_append: FORMAT=PS");
    ret = bcf_hdr_append(hdr, "##FORMAT=<ID=SB,Number=4,Type=Integer,Description=\"Per-sample component statistics which comprise the Fisher's Exact Test to detect strand bias.\">");
    CHECK_CONDITION_EXIT(ret != 0, "bcf_hdr_append: FORMAT=SB");
    // clang-format on

    std::string command_line = generate_command_line(args->tool_name, args->command_line.c_str());
    ret = bcf_hdr_append(hdr, command_line.c_str());
    CHECK_CONDITION_EXIT(ret != 0, "bcf_hdr_append: LUSHCommandLine=%s", command_line.c_str());

    generate_gvcf_block(hdr, args);

    // clang-format off
    ret = bcf_hdr_append(hdr, "##INFO=<ID=BaseQRankSum,Number=1,Type=Float,Description=\"Z-score from Wilcoxon rank sum test of Alt Vs. Ref base qualities\">");
    CHECK_CONDITION_EXIT(ret != 0, "bcf_hdr_append: INFO=BaseQRankSum");
    ret = bcf_hdr_append(hdr, "##INFO=<ID=DP,Number=1,Type=Integer,Description=\"Approximate read depth; some reads may have been filtered\">");
    CHECK_CONDITION_EXIT(ret != 0, "bcf_hdr_append: INFO=DP");
    ret = bcf_hdr_append(hdr, "##INFO=<ID=END,Number=1,Type=Integer,Description=\"Stop position of the interval\">");
    CHECK_CONDITION_EXIT(ret != 0, "bcf_hdr_append: INFO=END");
    ret = bcf_hdr_append(hdr, "##INFO=<ID=ExcessHet,Number=1,Type=Float,Description=\"Phred-scaled p-value for exact test of excess heterozygosity\">");
    CHECK_CONDITION_EXIT(ret != 0, "bcf_hdr_append: INFO=ExcessHet");
    ret = bcf_hdr_append(hdr, "##INFO=<ID=InbreedingCoeff,Number=1,Type=Float,Description=\"Inbreeding coefficient as estimated from the genotype likelihoods per-sample when compared against the Hardy-Weinberg expectation\">");
    CHECK_CONDITION_EXIT(ret != 0, "bcf_hdr_append: INFO=InbreedingCoeff");
    ret = bcf_hdr_append(hdr, "##INFO=<ID=MLEAC,Number=A,Type=Integer,Description=\"Maximum likelihood expectation (MLE) for the allele counts (not necessarily the same as the AC), for each ALT allele, in the same order as listed\">");
    CHECK_CONDITION_EXIT(ret != 0, "bcf_hdr_append: INFO=MLEAC");
    ret = bcf_hdr_append(hdr, "##INFO=<ID=MLEAF,Number=A,Type=Float,Description=\"Maximum likelihood expectation (MLE) for the allele frequency (not necessarily the same as the AF), for each ALT allele, in the same order as listed\">");
    CHECK_CONDITION_EXIT(ret != 0, "bcf_hdr_append: INFO=MLEAF");
    ret = bcf_hdr_append(hdr, "##INFO=<ID=MQRankSum,Number=1,Type=Float,Description=\"Z-score From Wilcoxon rank sum test of Alt vs. Ref read mapping qualities\">");
    CHECK_CONDITION_EXIT(ret != 0, "bcf_hdr_append: INFO=MQRankSum");
    ret = bcf_hdr_append(hdr, "##INFO=<ID=RAW_MQandDP,Number=2,Type=Integer,Description=\"Raw data (sum of squared MQ and total depth) for improved RMS Mapping Quality calculation. Incompatible with deprecated RAW_MQ formulation.\">");
    CHECK_CONDITION_EXIT(ret != 0, "bcf_hdr_append: INFO=RAW_MQandDP");
    ret = bcf_hdr_append(hdr, "##INFO=<ID=ReadPosRankSum,Number=1,Type=Float,Description=\"Z-score from Wilcoxon rank sum test of Alt vs. Ref read position bias\">");
    CHECK_CONDITION_EXIT(ret != 0, "bcf_hdr_append: INFO=ReadPosRankSum");
    // clang-format on

    generate_samples(hdr, bam_hdr);
    generate_contig(hdr, bam_hdr);

    ret = bcf_hdr_append(hdr, "##source=HaplotypeCaller");
    CHECK_CONDITION_EXIT(ret != 0, "bcf_hdr_append: source=HaplotypeCaller");

    return hdr;
}

void Writer::generate_gvcf_block(bcf_hdr_t *hdr, pHCArgs args)
{
    int32_t pre = 0, ret;
    std::string block_line;
    for (int32_t b : args->gvcf_gq_bands) {
        block_line = generate_gvcf_block_line(pre, b);
        ret = bcf_hdr_append(hdr, block_line.c_str());
        CHECK_CONDITION_EXIT(ret != 0, "bcf_hdr_append: GVCFBlock={}", block_line.c_str());
        pre = b;
    }
    if (pre == MAX_GENOTYPE_QUAL) {
        block_line = generate_gvcf_block_line(pre, pre + 1);
        ret = bcf_hdr_append(hdr, block_line.c_str());
        CHECK_CONDITION_EXIT(ret != 0, "bcf_hdr_append: GVCFBlock={}", block_line.c_str());
    }
}

void Writer::generate_contig(bcf_hdr_t *hdr, bam_hdr_t *bam_hdr)
{
    int32_t ret;
    std::string contig_line;
    for (int32_t i = 0; i < bam_hdr->n_targets; ++i) {
        contig_line = generate_contig_line(bam_hdr->target_name[i], bam_hdr->target_len[i]);
        ret = bcf_hdr_append(hdr, contig_line.c_str());
        CHECK_CONDITION_EXIT(ret != 0, "bcf_hdr_append: contig={}", contig_line.c_str());
    }
}

void Writer::generate_samples(bcf_hdr_t *hdr, bam_hdr_t *bam_hdr)
{
    const char *q, *p = strstr(strstr(bam_hdr->text, "@RG"), "SM");
    for (p += 3, q = p; *q && *q != '\t' && *q != '\n'; ++q) {}
    std::string sample_name(p, q - p);
    int32_t ret = bcf_hdr_add_sample(hdr, sample_name.c_str());
    CHECK_CONDITION_EXIT(ret != 0, "bcf_hdr_add_sample");
}

}  // namespace rovaca