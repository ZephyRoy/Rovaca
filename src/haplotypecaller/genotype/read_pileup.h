/**
 * @brief Represents a pileup of reads at a given position.
 * @date 2022-06-12
 */
#ifndef ROVACA_HC_READ_PLIEUP_H_
#define ROVACA_HC_READ_PLIEUP_H_

#include <list>
#include <memory_resource>

#include "forward.h"
#include "genotype_macors.h"
#include "interface/interface_locatable.hpp"

namespace rovaca
{

#define FLAG_STATUS   2
#define HC_BASE_Q_MAX 64
enum FlagType { FS_REF = 0, FS_NON_REF };

class ReadPileup : public InterfaceLocatable
{
    /* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
private:
    int32_t _tid{INVALID_INT};
    int64_t _start{INVALID_INT}, _stop{INVALID_INT};
    std::pmr::list<pPileupElement> _pileup_elements;

    /* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
public:
    uint32_t qual_his[FLAG_STATUS][HC_BASE_Q_MAX]{};
    uint8_t qual_max[FLAG_STATUS]{};

    /* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
public:
    static pReadPileup create(int32_t contig, int64_t start, int64_t stop, pMemoryPool pool);

    DISALLOW_COPY_AND_ASSIGN(ReadPileup);
    ~ReadPileup() override = default;

    int32_t get_tid() const override { return _tid; }
    int64_t get_start() const override { return _start; }
    int64_t get_stop() const override { return _stop; }
    void set_tid(int32_t tid) override { _tid = tid; }
    void set_start(int64_t start) override { _start = start; }
    void set_stop(int64_t stop) override { _stop = stop; }

    void add(pPileupElement ele) { _pileup_elements.push_back(ele); }
    std::pmr::list<pPileupElement>& get_pileup_element() { return _pileup_elements; }

    /* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
private:
    ReadPileup(int32_t contig, int64_t start, int64_t stop, pMemoryPool pool);

};  // ReadPileup end

}  // namespace rovaca

#endif  // ROVACA_HC_READ_PLIEUP_H_
