/******************************************************************************************************************
* @copyright  : Boost Software License
* @file       : multiple_core_pair_map_index.h
* @brief      : a index for temporal vertex
* @version    : 1.0
* @date       : 2021/08/09
******************************************************************************************************************/

#pragma once
#include "multiple_core_utility.h"

namespace scnu{
    class multiple_core_pair_map_index {
    public:
        explicit multiple_core_pair_map_index();

        explicit multiple_core_pair_map_index(const shared_ptr<multiple_core_pair_map_index>& other_temporal_vertex_index);

        ~multiple_core_pair_map_index() = default;

        void clear();

        bool compare(const shared_ptr<multiple_core_pair_map_index>& other_temporal_vertex_index);

        bool count(uint32_t k, uint32_t h);

        bool empty();

        bool equal(uint32_t k, uint32_t h);

        shared_ptr<unordered_map<uint32_t, uint32_t>> get_left_map();

        shared_ptr<unordered_map<uint32_t, uint32_t>> get_right_map();

        uint32_t get_memory_cost();

        uint32_t get_core_size();

        uint32_t get_delta();

        uint32_t get_k(uint32_t h);

        uint32_t get_h(uint32_t k);

        void merge_insert(const shared_ptr<multiple_core_pair_map_index>& other_index);

        void merge_remove(const shared_ptr<multiple_core_pair_map_index>& other_index);

        void insert(uint32_t k, uint32_t h);

        void left_insert(uint32_t key, uint32_t value);

        void right_insert(uint32_t key, uint32_t value);

        void remove(uint32_t k, uint32_t h);

        void left_remove(uint32_t key);

        void left_remove(uint32_t key, uint32_t value);

        void right_remove(uint32_t key);

        void right_remove(uint32_t key, uint32_t value);

        void set(uint32_t key, uint32_t value);

        void set_left(uint32_t key, uint32_t value);

        void set_right(uint32_t key, uint32_t value);



    private:
        shared_ptr<unordered_map<uint32_t, uint32_t>> left_map;
        shared_ptr<unordered_map<uint32_t, uint32_t>> right_map;
    };
}


