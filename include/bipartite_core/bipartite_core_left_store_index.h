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
    class bipartite_core_left_store_index {
    public:
        bipartite_core_left_store_index();

        explicit bipartite_core_left_store_index(const shared_ptr<bipartite_core_left_store_index>& other_left_node_index);

        bool compare(shared_ptr<bipartite_core_left_store_index> &another_left_node_index);

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

        bool operator==(shared_ptr<bipartite_core_left_store_index> &another_left_node_index);

        bool operator!=(shared_ptr<bipartite_core_left_store_index> &another_left_node_index);

    private:
        shared_ptr<map<uint32_t , uint32_t>> index_map;
    };
}
