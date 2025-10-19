#ifndef ROVACA_HC_INTERFACE_LOCATABLE_H_
#define ROVACA_HC_INTERFACE_LOCATABLE_H_
#include <cstdint>

namespace rovaca
{

class InterfaceLocatable
{
    /* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
public:
    virtual ~InterfaceLocatable() = default;
    virtual int32_t get_tid() const = 0;
    virtual int64_t get_start() const = 0;
    virtual int64_t get_stop() const = 0;
    virtual void set_tid(int32_t tid) = 0;
    virtual void set_start(int64_t start) = 0;
    virtual void set_stop(int64_t stop) = 0;

    /*! @brief Determines whether this interval overlaps the provided interval */
    virtual bool overlaps(const InterfaceLocatable& other) const { return within_distance_of(other, 0); }

    int64_t get_length_on_reference() const { return get_stop() - get_start() + 1; }

    /*! @brief Determines whether this interval comes within of overlapping the provided interval */
    bool within_distance_of(const InterfaceLocatable& other, int64_t distance) const
    {
        int64_t this_start = get_start(), this_stop = get_stop();
        int64_t other_start = other.get_start(), other_stop = other.get_stop();
        return contigs_match(other) && overlaps(this_start, this_stop, other_start - distance, other_stop + distance);
    }

    /** @brief Determines whether this interval contains provided interval */
    bool contains(const InterfaceLocatable& other) const
    {
        const int64_t &this_start = get_start(), this_stop = get_stop();
        const int64_t &other_start = other.get_start(), other_stop = other.get_stop();
        return contigs_match(other) && encloses(this_start, this_stop, other_start, other_stop);
    }

    /*! @brief Determine if this is on the same chr as other */
    bool contigs_match(const InterfaceLocatable& other) const { return this->get_tid() == other.get_tid(); }

    /* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
private:
    /*! @brief Checks to see if the two sets of coordinates have any overlap */
    static bool overlaps(int64_t start, int64_t end, int64_t start2, int64_t end2)
    {
        return (start2 >= start && start2 <= end) || (end2 >= start && end2 <= end) || encloses(start2, end2, start, end);
    }

    /*! @brief outter contains inner */
    static bool encloses(int64_t outer_start, int64_t outer_end, int64_t inner_start, int64_t inner_end)
    {
        return inner_start >= outer_start && inner_end <= outer_end;
    }
};

}  // namespace rovaca

namespace std
{

template <>
struct less<rovaca::pInterfaceLocatable>
{
    bool operator()(rovaca::pInterfaceLocatable l, rovaca::pInterfaceLocatable r) const { return l->get_stop() < r->get_stop(); }
};

}  // namespace std

#endif  // ROVACA_HC_INTERFACE_LOCATABLE_H_
