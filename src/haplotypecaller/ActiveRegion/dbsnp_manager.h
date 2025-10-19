#ifndef DBSNP_MANAGER_H
#define DBSNP_MANAGER_H
#include <condition_variable>
#include <map>
#include <memory>
#include <mutex>

#include "vcf_loader.h"

class BedLoader;

class DbsnpManager
{
private:
    int32_t undo_num_;
    int32_t prefetch_num_;
    volatile int32_t finish_;
    BedLoader *bed_loader_;
    VcfLoader vcf_loader_;
    contig_info_t *contig_info_;
    std::map<int32_t, std::shared_ptr<std::vector<bcf1_t *>>> vcf_;

    std::mutex mutex_;
    std::condition_variable not_full_;
    std::condition_variable has_data_;

public:
    DbsnpManager()
        : undo_num_{0}
        , prefetch_num_{0}
        , finish_{0}
        , bed_loader_{nullptr}
        , vcf_loader_{}
        , contig_info_{nullptr}
        , vcf_{}
        , mutex_{}
        , not_full_{}
        , has_data_{}
    {}

    ~DbsnpManager() = default;

    /// @brief
    /// @param vcf_file
    /// @param bl 此参数可以为nullptr
    /// @param ci ReferenceManager中的contig_info_t，与region相同的顺序获取vcf
    /// @param tp 公用线程池，用于vcf文件读取
    /// @param prefetch run函数提前读取的染色体数量
    void initialization(const std::string &vcf_file, BedLoader *bl, contig_info_t *ci, htsThreadPool *tp, int32_t prefetch);

    /// @brief 资源释放，work结束后调用，否则会提前释放资源，调用此函数时，理论上所有发出的bcf1_t*都应该返回
    void terminalization();

    /// @brief 单独线程运行，持续读取直到
    void run();

    /// @brief 根据tid获取对应的数据
    /// @param tid
    /// @return
    std::shared_ptr<std::vector<bcf1_t *>> get(int32_t tid);

    /// @brief region携带shared_ptr发出数据，此处的备份可以删除
    /// @param tid
    void notify(int32_t tid);

    /// @brief 有bed的模式，读完bed对应的数据，通知db读取线程结束
    void update_finish()
    {
        finish_ = 1;
        not_full_.notify_one();
    }
};

#endif  // DBSNP_MANAGER_H