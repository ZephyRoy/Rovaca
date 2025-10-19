#ifndef ROVACA_HC_BAM_DATA_POOL_H_
#define ROVACA_HC_BAM_DATA_POOL_H_
#include <cstdint>
#include <cstdlib>

#include "htslib/sam.h"

namespace rovaca
{

/*!
 * 此类给bam对bam1_t中data做内存池管理，bam1结构管理由调用者管理
 * alloc_bam_struct_pool 为从内存池里面申请空间，如果内存池空间不够调用 bam1_init()
 * 其中assign_bam_mem_pool指定bam1_t中data对应内存池未使用的空间，若后续操作sam_read1,bam1_set1函数需要的空间大于剩余内存空间，htslib内部则会堆上申请。
 * 当然，bam_destroy1释放bam,他不会释放内存池的空间。
 * 注意：调用顺序 alloc_bam_struct_pool -> assign_bam_mem_pool -> 多次 sam_read1 | bam1_set1 | bam_copy1 -> peek_bam ;
 * 调用 sam_read1 | bam1_set1 | bam_copy1 函数后，不能调用assign_bam_mem_pool 函数
 */
class BamDataPool
{
    /* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
private:
    uint32_t _max;
    uint32_t _used;
    char *_data;

    /* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
public:
    explicit BamDataPool(uint32_t n)
        : _max(n)
        , _used(0)
        , _data((char *)malloc(n))
    {}
    ~BamDataPool() { free(_data); }

    // 从申请内存池申请bam1的struct空间,同时将bam的data指定内存池的剩余空间,设置其m_data为内存池剩余空间
    bam1_t *alloc_bam_struct_pool()
    {
        if (_used + sizeof(bam1_t) > _max) {
            return bam_init1();
        }

        // 分配 struct
        auto *ret = (bam1_t *)(_data + _used);
        _used += sizeof(bam1_t);

        // 赋值 data
        ret->data = (uint8_t *)(_data + _used);
        ret->m_data = _max - _used;

        bam_set_mempolicy(ret, BAM_USER_OWNS_STRUCT | BAM_USER_OWNS_DATA);

        return ret;
    }

    // 将bam的data指定内存池的剩余空间
    // int assign_bam_mem_pool(bam1_t *bam);

    // 调用htslib相关接口后，计算内存池还剩多少空间
    void peek_bam(bam1_t *bam)
    {
        if ((bam_get_mempolicy(bam) & BAM_USER_OWNS_DATA) == 0) {
            return;
        }
        uint32_t used_len = ((uint32_u)bam->l_data + 7) & (~7U);
        if (used_len < bam->m_data) {
            bam->m_data = used_len;
        }
        _used += used_len;
    }
    uint32_t get_unsed_data_size() { return _max - _used; }

    /*!
     * @brief 在一个循环内重复利用某一段空间, 循环开始记录 _used, 循环结束还原
     */
    uint32_t get_used_data_size() const { return _used; }
    void reset(uint32_t used) { _used = used; }

    // 回收内存使用空间，注意那些要确认引用内存池的bam1_t不能再使用了
    void finalize() { _used = 0; }
};

}  // namespace rovaca

#endif  // ROVACA_HC_BAM_DATA_POOL_H_
