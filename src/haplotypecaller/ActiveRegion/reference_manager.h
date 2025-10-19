#ifndef __REFERENCE_MANAGER_H
#define __REFERENCE_MANAGER_H
#include <condition_variable>
#include <map>
#include <memory>
#include <mutex>
#include <thread>

#include "bed_loader.h"
#include "fasta_loader.h"

const int default_wgs_prefetch_number = 10;
const int default_wes_prefetch_number = 10;

class ReferenceManager
{
    static constexpr int padding = 500;

public:
    ReferenceManager(const std::string &reference, BedLoader *bed, int num)
        : fasta_loader(reference)
        , bed_loader(bed)
        , prefetch_number(num)
        , current_tid(0)
        , last_tid(0)
    {
        FastaLoader::get_fasta_dict(reference, &contig);
    }

    ReferenceManager(const std::string &reference, int num)
        : ReferenceManager(reference, nullptr, num)
    {}
    ReferenceManager(std::string &reference)
        : ReferenceManager(reference, nullptr, default_wgs_prefetch_number)
    {}

    ReferenceManager(std::string &reference, BedLoader *bed)
        : ReferenceManager(reference, bed, default_wes_prefetch_number)
    {}

    void run();

    contig_info_t &get_contig() { return contig; }

    std::shared_ptr<char> get(int tid, int &length);

    void pop();

    void assign_finish();

private:
    bool force_finish_{false};
    FastaLoader fasta_loader;
    BedLoader *bed_loader;
    size_t prefetch_number;
    std::map<int, std::shared_ptr<char>> ref_map;
    int current_tid;
    int last_tid;
    contig_info_t contig;
    std::condition_variable m_not_full_condition_;
    std::condition_variable m_new_comming_condition_;
    std::mutex m_mutex_;
};
#endif