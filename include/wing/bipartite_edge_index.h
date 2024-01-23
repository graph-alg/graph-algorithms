/******************************************************************************************************************
* @copyright  : Boost Software License
* @file       : bipartite_edge_index.h
* @brief      : An bipartite edge index
* @version    : 1.0
* @date       : 2021/7/9
******************************************************************************************************************/

#pragma once
#include "wing_utility.h"
namespace scnu{
    class bipartite_edge_index {
    public:
        bipartite_edge_index(const shared_ptr<abstract_bipartite_edge>& other_edge);

        shared_ptr<abstract_bipartite_edge> get_edge();

        shared_ptr<unordered_set<uint32_t>> get_left_vertex_id_set();

        shared_ptr<unordered_set<uint32_t>> get_right_vertex_id_set();

        void insert_dual_edge(const shared_ptr<abstract_bipartite_edge>& dual_edge);

        void remove_dual_edge(const shared_ptr<abstract_bipartite_edge>& dual_edge);

    private:
        shared_ptr<abstract_bipartite_edge> edge;
        shared_ptr<unordered_set<uint32_t>> left_vertex_id_set;
        shared_ptr<unordered_set<uint32_t>> right_vertex_id_set;
    };
}



