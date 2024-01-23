/******************************************************************************************************************
* @copyright  : Boost Software License
* @file       : bipartite_core.h
* @details    : auxiliary data structure for bipartite core
* @version    : 1.0
* @date       : 2021/02/01
******************************************************************************************************************/

#pragma once
#include "bipartite_core/bipartite_core_utility.h"

namespace scnu{
    class bipartite_core
    {
    public:
        bipartite_core();

        explicit bipartite_core(const shared_ptr<bipartite_core>& other_bipartite_core);

        bipartite_core(const shared_ptr<unordered_set<uint32_t>>& other_left_vertex_set,
                       const shared_ptr<unordered_set<uint32_t>>& other_right_vertex_set);

        bool count_left_vertex(uint32_t l);


        shared_ptr<unordered_set<uint32_t>> get_left_vertex_set();

        shared_ptr<unordered_set<uint32_t>> get_right_vertex_set();

        void insert_left_vertex(uint32_t l);

        void insert_right_vertex(uint32_t r);

        void remove_left_vertex(uint32_t l);

        void remove_right_vertex(uint32_t r);

        bool count_right_vertex(uint32_t r);

        bool empty();

    private:
        shared_ptr<unordered_set<uint32_t>> left_vertex_set;
        shared_ptr<unordered_set<uint32_t>> right_vertex_set;
    };
}


