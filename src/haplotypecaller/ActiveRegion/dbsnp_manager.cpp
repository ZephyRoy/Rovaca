#include "dbsnp_manager.h"

#include <algorithm>

#include "bed_loader.h"
#include "rovaca_logger.h"

void DbsnpManager::initialization(const std::string& vcf_file, BedLoader* bl, contig_info_t* ci, htsThreadPool* tp, int32_t prefetch)
{
    bed_loader_ = bl;
    contig_info_ = ci;

    prefetch_num_ = prefetch;

    vcf_loader_.initialization(vcf_file, tp, ci);
}

void DbsnpManager::terminalization() { vcf_loader_.terminalization(); }

void DbsnpManager::run()
{
    RovacaLogger::info("start reading dbsnp file");
    p_bed_intervals bi;
    int32_t len = static_cast<int32_t>(contig_info_->key.size());

    for (int32_t id{0}; id < len; ++id) {
        {
            // 缩小上锁范围
            std::unique_lock<std::mutex> lock{mutex_};
            while (undo_num_ == prefetch_num_ && !finish_) {
                not_full_.wait(lock);
            }
        }

        if (finish_) {
            break;
        }

        // bed中没有此chr，构造一个空数据
        if (bed_loader_ && !bed_loader_->get_bed_intervals().count(contig_info_->key[id])) {
            std::shared_ptr<std::vector<bcf1_t*>> data{new std::vector<bcf1_t*>{}};
            std::unique_lock<std::mutex> lock{mutex_};
            ++undo_num_;
            vcf_.insert({id, data});
            has_data_.notify_one();

            continue;
        }

        // 自定义删除器，回收 bcf1_t
        bi = bed_loader_ ? bed_loader_->get_bed_intervals().at(contig_info_->key[id]) : nullptr;
        std::vector<bcf1_t*> result = vcf_loader_.load_vcf(contig_info_->key.at(id).c_str(), bi);
        std::shared_ptr<std::vector<bcf1_t*>> data{new std::vector<bcf1_t*>{std::move(result)}, [&](const std::vector<bcf1_t*>* d) {
                                                       std::for_each(d->begin(), d->end(), [&](bcf1_t* b) { vcf_loader_.recover_bcf(b); });
                                                   }};

        std::unique_lock<std::mutex> lock{mutex_};
        ++undo_num_;
        vcf_.insert({id, data});
        has_data_.notify_one();
    }

    RovacaLogger::info("end reading dbsnp file");
}

std::shared_ptr<std::vector<bcf1_t*>> DbsnpManager::get(int32_t tid)
{
    if (!vcf_.count(tid)) {
        std::unique_lock<std::mutex> lock{mutex_};
        while (undo_num_ == 0 || !vcf_.count(tid)) {
            has_data_.wait(lock);
        }
    }

    return vcf_.at(tid);
}

void DbsnpManager::notify(int32_t tid)
{
    assert(std::size_t(1) == vcf_.count(tid));

    RovacaLogger::info("db tid {} finished", tid);

    std::unique_lock<std::mutex> lock{mutex_};
    vcf_.erase(tid);
    --undo_num_;
    not_full_.notify_one();
}
