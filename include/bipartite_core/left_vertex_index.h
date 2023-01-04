/******************************************************************************************************************
* @copyright  : Boost Software License
* @file       : bipartite_core_utility.hpp
* @details    : common head files for bipartite core
* @version    : 1.0
* @date       : 2021/02/01
******************************************************************************************************************/

#pragma once
#include "bipartite_core/bipartite_core_utility.h"

namespace scnu
{
    class left_vertex_index {
    public:
        left_vertex_index();

        explicit left_vertex_index(const shared_ptr<left_vertex_index>& other_left_node_index);

        bool compare(shared_ptr<left_vertex_index> &another_left_node_index);

        void clear();

        bool count(uint32_t i);

        bool count(uint32_t i, uint32_t j);

        bool empty();
        shared_ptr<map<uint32_t , uint32_t>> get_index_map();

        uint32_t get_j(uint32_t i);

        uint32_t get_maximal_i(uint32_t j);

        void insert(uint32_t i, uint32_t j);

        void remove(uint32_t i);

        void remove(uint32_t i, uint32_t j);

        void set(uint32_t i, uint32_t j);

        bool operator==(shared_ptr<left_vertex_index> &another_left_node_index);

        bool operator!=(shared_ptr<left_vertex_index> &another_left_node_index);

    private:
        shared_ptr<map<uint32_t , uint32_t>> index_map;
    };
}
