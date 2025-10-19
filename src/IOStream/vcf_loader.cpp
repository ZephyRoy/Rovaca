#include "vcf_loader.h"

#include "../rovaca_logger/rovaca_logger.h"

constexpr uint32_t s_init_capacity = 1000;

void VcfLoader::VcfPool::initialization()
{
    static constexpr uint32_t s_init_capacity = 100000;
    data_.reserve(s_init_capacity);
}

void VcfLoader::VcfPool::terminalization()
{
    for (bcf1_t* b : data_) {
        bcf_clear(b);
        bcf_destroy(b);
    }
    data_.clear();
}

bool VcfLoader::VcfPool::push(bcf1_t* b)
{
    if (nullptr == b) {
        return false;
    }

    {
        std::lock_guard lock{push_mutex_};
        data_.push_back(b);
    }

    return true;
}

bool VcfLoader::VcfPool::pop(bcf1_t** b)
{
    {
        std::lock_guard lock{push_mutex_};
        if (!data_.empty()) {
            *b = data_.back();
            data_.pop_back();
            return true;
        }
    }

    bcf1_t* result = bcf_init();
    if (nullptr == result) {
        RovacaLogger::error("alloc bcf1_t error.");
        return false;
    }

    *b = result;
    return true;
}

void VcfLoader::initialization(const std::string& vcf_db, htsThreadPool* tp, contig_info_t* info)
{
    fd_ = bcf_open(vcf_db.c_str(), "r");
    if (nullptr == fd_) {
        RovacaLogger::error("open dbsnp file error: {}.", vcf_db);
        exit(-1);
    }

    hts_set_thread_pool(fd_, tp);

    hdr_ = bcf_hdr_read(fd_);
    if (nullptr == hdr_) {
        RovacaLogger::error("read dbsnp header error.");
        exit(-1);
    }

    // Check contig matching status
    for (int32_t id{0}; id < hdr_->n[BCF_DT_CTG]; ++id) {
        if (0 != strcmp(info->key.at(id).c_str(), hdr_->id[BCF_DT_CTG][id].key)) {
            RovacaLogger::error("dbsnp contigs not matching fasta.");
            exit(-1);
        }
    }

    const std::string idx_name = vcf_db + ".tbi";
    tbx_ = tbx_index_load(idx_name.c_str());
    if (nullptr == tbx_) {
        RovacaLogger::error("open dbsnp index error: {}.", idx_name);
        exit(-1);
    }

    pool_.initialization();

    ks_resize(&cache_line_, s_init_capacity);

    seqnames_ = bcf_hdr_seqnames(hdr_, &n_seqname_);

    for (int32_t i{0}; i < n_seqname_; ++i) {
        seqnames_set_.insert(std::string(seqnames_[i]));
    }
}

void VcfLoader::terminalization()
{
    free(seqnames_);
    pool_.terminalization();
    ks_free(&cache_line_);
    tbx_destroy(tbx_);
    bcf_hdr_destroy(hdr_);
    bcf_close(fd_);
}

std::vector<bcf1_t*> VcfLoader::load_vcf(const char* chr, p_bed_intervals bed)
{
    if (!seqnames_set_.count(chr)) {
        return {};
    }

    std::vector<bcf1_t*> result;
    result.reserve(s_init_capacity);

    int32_t ret = 0;
    bcf1_t* b = nullptr;
    hts_itr_t* itr = nullptr;

    for (int i{0}; i < (nullptr == bed ? 1 : bed->n); ++i) {
        itr = (nullptr == bed ? tbx_itr_querys(tbx_, chr) : tbx_itr_queryi(tbx_, bcf_hdr_name2id(hdr_, chr), bed->start[i], bed->end[i]));
        if (nullptr == itr) {
            RovacaLogger::error("query dbsnp error: {}", chr);
            exit(-1);
        }

        if (ret < -1) {
            RovacaLogger::error("error reading DB file: {}", chr);
            exit(-1);
        }

        while ((ret = tbx_itr_next(fd_, tbx_, itr, &cache_line_)) >= 0) {
            if (!pool_.pop(&b)) {
                RovacaLogger::error("pop bcf1_t error");
                exit(-1);
            }

            if ((ret = vcf_parse1(&cache_line_, hdr_, b)) < -1) {
                RovacaLogger::error("vcf_parse1 error: {}", chr);
                exit(-1);
            }

            result.push_back(b);
        }

        bcf_itr_destroy(itr);
    }

    return result;
}
