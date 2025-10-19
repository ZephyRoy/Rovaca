#ifndef BAM_POO_H
#define BAM_POO_H

#include <cassert>
#include <cstddef>
#include <iostream>
#include <stdexcept>
#include <vector>

#include "htslib/sam.h"

class Bam1Pool
{
public:
    Bam1Pool(size_t object_count, size_t data_size)
        : object_count_(object_count)
        , object_size_(data_size)
    {
        assert(object_count_ > 0);
        assert(data_size > 0);

        // Initialize the objects and add them to the free_list_
        free_list_.reserve(object_count_);
        for (size_t i = 0; i < object_count_; i++) {
            try {
                bam1_t* object = bam_init1();
                object->mempolicy = 0;
                object->m_data = object_size_;
                object->data = (uint8_t*)malloc(object_size_);
                free_list_.push_back(object);
            }
            catch (...) {
                for (auto item : free_list_) {
                    bam_destroy1(item);
                }
                throw std::runtime_error("Bam1Pool alloc error when init.");
            }
        }
    }

    ~Bam1Pool()
    {
        // Destroy the objects, ensure all source are recycled.
        for (size_t i = 0; i < free_list_.size(); i++) {
            free(free_list_[i]->data);
            free(free_list_[i]);
        }
        free_list_.clear();
    }

    bam1_t* allocate() noexcept
    {
        if (free_list_.empty()) {
            object_count_++;
            bam1_t* object = bam_init1();
            object->mempolicy = 0;
            object->m_data = object_size_;
            object->data = (uint8_t*)malloc(object_size_);
            return object;
        }

        bam1_t* object = free_list_.back();
        free_list_.pop_back();
        return object;
    }

    void deallocate(bam1_t* object) noexcept
    {
        if (object == NULL) {
            std::cerr << "deallocater get a nullptr in Bam1Pool" << std::endl;
            object = bam_init1();
            object->mempolicy = 0;
            object->m_data = object_size_;
            object->data = (uint8_t*)malloc(object_size_);
        }

        // uint32_t m_data = object->m_data;
        // uint8_t* data_ptr = object->data;
        // memset(data_ptr, 0, object->l_data * sizeof(uint8_t));
        // memset(object, 0, sizeof(bam1_t));
        // object->data = data_ptr;
        // object->m_data = m_data;
        free_list_.push_back(object);
    }

private:
    size_t object_count_;
    size_t object_size_;
    std::vector<bam1_t*> free_list_;

    Bam1Pool(const Bam1Pool&) = delete;
    Bam1Pool& operator=(const Bam1Pool&) = delete;
    Bam1Pool(Bam1Pool&&) = delete;
    Bam1Pool& operator=(Bam1Pool&&) = delete;
};

#endif  // BAM_POO_H