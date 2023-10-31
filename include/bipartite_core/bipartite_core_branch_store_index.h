/******************************************************************************************************************
* @copyright  : Boost Software License
* @file       : bipartite_vertex_index.h
* @details    : common head files for bipartite core
* @version    : 1.0
* @date       : 2023/06/13
******************************************************************************************************************/

#pragma once

#include "bipartite_core/bipartite_core_utility.h"

namespace scnu {
    class bipartite_core_branch_store_index {
    public:
        explicit bipartite_core_branch_store_index();

        explicit bipartite_core_branch_store_index(const shared_ptr<bipartite_core_branch_store_index> &other_temporal_vertex_index);

        ~bipartite_core_branch_store_index() = default;

        void clear();

        bool compare(const shared_ptr<bipartite_core_branch_store_index> &other_temporal_vertex_index);

        bool count(uint32_t i, uint32_t j);

        bool empty();

        bool equal(uint32_t i, uint32_t j);

        shared_ptr<unordered_map<uint32_t, uint32_t>> get_left_map();

        shared_ptr<unordered_map<uint32_t, uint32_t>> get_right_map();

        uint32_t get_memory_cost();

        uint32_t get_k();

        uint32_t get_i(uint32_t j);

        uint32_t get_j(uint32_t i);

        void merge_insert(const shared_ptr<bipartite_core_branch_store_index> &other_index);

        void merge_remove(const shared_ptr<bipartite_core_branch_store_index> &other_index);

        void insert(uint32_t i, uint32_t j);

        void left_insert(uint32_t key, uint32_t value);

        void right_insert(uint32_t key, uint32_t value);

        void remove(uint32_t i, uint32_t j);

        void left_remove(uint32_t key);

        void left_remove(uint32_t key, uint32_t value);

        void right_remove(uint32_t key);

        void right_remove(uint32_t key, uint32_t value);

        void set(uint32_t i, uint32_t j);

        void set_left(uint32_t key, uint32_t value);

        void set_right(uint32_t key, uint32_t value);

    private:
        shared_ptr<unordered_map<uint32_t, uint32_t>> left_map;
        shared_ptr<unordered_map<uint32_t, uint32_t>> right_map;
    };
}

