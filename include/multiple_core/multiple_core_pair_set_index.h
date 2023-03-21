/******************************************************************************************************************
* @copyright  : Boost Software License
* @file       : multiple_core_pair_set_index.h
* @brief      : a index for temporal vertex
* @version    : 1.0
* @date       : 2022/03/07
******************************************************************************************************************/

#pragma once
#include "multiple_core_utility.h"

namespace scnu{
    class multiple_core_pair_set_index {
    public:
        multiple_core_pair_set_index();

        bool compare(const shared_ptr<multiple_core_pair_set_index>& other_multiple_core_pair_index);

        bool count(uint32_t k, uint32_t h);

        shared_ptr<unordered_set<pair<uint32_t, uint32_t>, hash_pair>> get_core_pair_set();

        uint32_t get_memory_cost();

        void insert(uint32_t k,uint32_t h);

    private:
        shared_ptr<unordered_set<pair<uint32_t, uint32_t>, hash_pair>> core_pair_set;

    };
}




