/******************************************************************************************************************
* @copyright  : Boost Software License
* @file       : temporal_vertex_index.h
* @brief      : a index for temporal edge
* @version    : 1.0
* @date       : 2021/10/24
******************************************************************************************************************/

#pragma once
#include "multiple_core_utility.h"

namespace scnu{
    class multiple_core_number_set_index {
    public:
        multiple_core_number_set_index();

        ~multiple_core_number_set_index() = default;

        bool compare(const shared_ptr<multiple_core_number_set_index>& other_multiple_core_set_index);

        bool count(uint32_t k, uint32_t h);

        shared_ptr<unordered_set<pair<uint32_t, uint32_t>, hash_pair>> get_multiple_core_map();

        uint32_t get_memory_cost();

        void insert(uint32_t k, uint32_t h);

    private:
        shared_ptr<unordered_set<pair<uint32_t, uint32_t>, hash_pair>> core_number_map;
    };
}



