#include "ActiveBase.h"

bool ActiveBase::compulte_all_likelihood()
{
    for (hts_pos_t iter = m_actual_start; iter < m_actual_end; iter++) {
        const int offset = get_actual_offset(m_tid, iter);
        if (offset == -1) continue;
        compute_slot_likelihood(offset);
    }
    return true;
}

int ActiveBase::get_offset(int tid, hts_pos_t pos)
{
    if (m_tid != tid) return -1;
    if (pos < m_start || pos >= m_end) return -1;
    if (!m_bit->test(pos - m_start)) return -1;
    return pos - m_start;
}

int ActiveBase::get_actual_offset(int tid, hts_pos_t pos)
{
    hts_pos_t ret = (pos - m_actual_start);
    // if ((uint64_t)ret >= (uint64_t)(m_actual_end - m_actual_start)) return -1;
    if (m_tid != tid || pos < m_actual_start || pos >= m_actual_end) return -1;
    if (m_bit && !m_bit->test(ret)) return -1;
    return ret;
}

int ActiveBase::offset_transform_actual(int offset)
{
    hts_pos_t pos = m_start + offset;
    if (pos < m_actual_start || pos >= m_actual_end) return -1;
    int real_offset = pos - m_actual_start;
    if (m_bit && !m_bit->test(real_offset)) return -1;
    return real_offset;
}