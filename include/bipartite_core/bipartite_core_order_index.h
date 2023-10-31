/******************************************************************************************************************
* @copyright  : Boost Software License
* @file       : bipartite_core_order.h
* @details    : common head files for bipartite core
* @version    : 1.0
* @date       : 2023/06/13
******************************************************************************************************************/
#pragma once
#include "bipartite_core/bipartite_core_utility.h"

namespace scnu{
    class bipartite_core_order_index {
    public:
        bipartite_core_order_index();

        explicit bipartite_core_order_index(const shared_ptr<bipartite_core_order_index>& other_core_order_index);

        shared_ptr<unordered_map<uint32_t, shared_ptr<unordered_map<uint32_t, shared_ptr<extend_list<int, uint32_t>>>>>>
        get_left_map();

        shared_ptr<unordered_map<uint32_t, shared_ptr<extend_list<int, uint32_t>>>> get_left_map(uint32_t k);

        shared_ptr<unordered_map<uint32_t, shared_ptr<extend_list<int, uint32_t>>>> get_middle_map();

        shared_ptr<unordered_map<uint32_t, shared_ptr<unordered_map<uint32_t , shared_ptr<extend_list<int, uint32_t>>>>>>
        get_right_map();

        shared_ptr<unordered_map<uint32_t, shared_ptr<extend_list<int, uint32_t>>>> get_right_map(uint32_t k);

        void insert_left(uint32_t k,
                         const shared_ptr<unordered_map<uint32_t, shared_ptr<extend_list<int, uint32_t>>>> &order_map);

        void insert_right(uint32_t k,
                          const shared_ptr<unordered_map<uint32_t, shared_ptr<extend_list<int, uint32_t>>>> &order_map);

    private:

        shared_ptr<unordered_map<uint32_t, shared_ptr<extend_list<int, uint32_t>>>> middle_map;

        shared_ptr<unordered_map<uint32_t, shared_ptr<unordered_map<uint32_t, shared_ptr<extend_list<int, uint32_t>>>>>> left_map;

        shared_ptr<unordered_map<uint32_t, shared_ptr<unordered_map<uint32_t, shared_ptr<extend_list<int, uint32_t>>>>>> right_map;
    };
}

