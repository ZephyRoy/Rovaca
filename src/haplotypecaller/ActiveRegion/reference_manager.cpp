#include "reference_manager.h"

#include <algorithm>
#include <iostream>

#include "rovaca_logger.h"
void ReferenceManager::run()
{
    while (current_tid < (int)contig.key.size() && !force_finish_) {
        std::unique_lock<std::mutex> lock(m_mutex_);
        while (prefetch_number == ref_map.size() && !force_finish_) {
            m_not_full_condition_.wait(lock);
        }
        if (force_finish_) {
            break;
        }
        if (!bed_loader) {
            char *raw_chr_bases = fasta_loader.fetch_target_seq(current_tid, 0, contig.idict[current_tid]);
            std::shared_ptr<char> chr_bases(raw_chr_bases, free);
            ref_map.insert({current_tid, chr_bases});
        }
        else {
            const std::vector<std::string> &keys = bed_loader->get_bed_keys();
            const std::string &chrom = contig.key[current_tid];
            if (std::find(keys.begin(), keys.end(), chrom) == keys.end()) {
                // ref_map.insert({current_tid, nullptr});
            }
            else {
                //
                const std::map<std::string, p_bed_intervals> &intervals = bed_loader->get_bed_intervals();
                p_bed_intervals bi = intervals.at(chrom);
                p_bed_intervals new_bi = BedLoader::init_bed_intervals();
                char *raw_chr_bases = (char *)malloc(sizeof(char) * contig.idict[current_tid]);

                for (int i = 0; i < bi->n; i++) {
                    if (i == 0) {
                        new_bi->start[0] = bi->start[i] - padding < 0 ? 0 : bi->start[i] - padding;
                        new_bi->end[0] =
                            bi->end[i] + padding > (hts_pos_t)contig.idict[current_tid] ? contig.idict[current_tid] : bi->end[i] + padding;
                        new_bi->n++;
                    }
                    else if (new_bi->end[new_bi->n - 1] > bi->start[i] - padding) {
                        new_bi->end[new_bi->n - 1] = bi->end[i] + padding;
                    }
                    else {
                        if (new_bi->m == new_bi->n - 1) {
                            new_bi->m = new_bi->m << 1;
                            new_bi->start = (hts_pos_t *)realloc(new_bi->start, sizeof(hts_pos_t) * new_bi->m);
                            new_bi->end = (hts_pos_t *)realloc(new_bi->end, sizeof(hts_pos_t) * new_bi->m);
                        }
                        new_bi->start[new_bi->n] = bi->start[i] - padding;
                        new_bi->end[new_bi->n] =
                            bi->end[i] + padding > (hts_pos_t)contig.idict[current_tid] ? contig.idict[current_tid] : bi->end[i] + padding;
                        ;
                        new_bi->n++;
                    }
                }
                for (int i = 0; i < new_bi->n; i++) {
                    char *refbase = fasta_loader.fetch_target_seq(current_tid, new_bi->start[i], new_bi->end[i]);
                    memcpy(raw_chr_bases + new_bi->start[i], refbase, (new_bi->end[i] - new_bi->start[i]) * sizeof(char));
                    free(refbase);
                }
                std::shared_ptr<char> chr_bases(raw_chr_bases, free);
                ref_map.insert({current_tid, chr_bases});
                free(new_bi->start);
                free(new_bi->end);
                free(new_bi);
            }
        }
        current_tid++;
        m_new_comming_condition_.notify_one();
    }
    return;
}

/*当前框架不用考虑退出问题*/
std::shared_ptr<char> ReferenceManager::get(int tid, int &length)
{
    if (current_tid > tid) {
        std::shared_ptr<char> chr_bases = ref_map[tid];
        length = contig.idict[tid];
        return chr_bases;
    }
    std::unique_lock<std::mutex> lock(m_mutex_);
    while (current_tid <= tid) {
        m_new_comming_condition_.wait(lock);
    }
    std::shared_ptr<char> chr_bases = ref_map[tid];
    length = contig.idict[tid];
    return chr_bases;
}
/*当前框架不用考虑退出问题*/
void ReferenceManager::pop()
{
    std::unique_lock<std::mutex> lock(m_mutex_);
    while (current_tid <= last_tid) {
        m_new_comming_condition_.wait(lock);
    }
    // std::shared_ptr<char> chr_bases = ref_map[last_tid];
    ref_map.erase(last_tid);
    last_tid++;
    m_not_full_condition_.notify_one();
}

void ReferenceManager::assign_finish()
{
    std::unique_lock<std::mutex> lock(m_mutex_);
    force_finish_ = true;
    m_not_full_condition_.notify_one();
}