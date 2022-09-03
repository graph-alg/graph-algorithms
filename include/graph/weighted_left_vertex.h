/******************************************************************************************************************
* @copyright  : Boost Software License
* @file       : weighted_left_vertex.h
* @brief      : a left vertex with weighted edges
* @version    : 1.0
* @date       : 2022/06/17
******************************************************************************************************************/

#pragma once
#include "graph/weighted_bipartite_edge.h"

namespace  scnu
{
    class weighted_left_vertex
    {
    public:
        explicit weighted_left_vertex(uint32_t other_left_vertex_id);

        explicit weighted_left_vertex(const shared_ptr<weighted_left_vertex>& left_vertex);

        virtual ~weighted_left_vertex() = default;

        [[nodiscard]] uint32_t get_left_vertex_id() const;

        shared_ptr<weighted_bipartite_edge> get_edge(uint32_t right_vertex_id);

        shared_ptr<unordered_map<uint32_t,shared_ptr<weighted_bipartite_edge>>> get_edge_map();

        uint32_t get_degree();

        void insert_edge(uint32_t right_vertex_id, const shared_ptr<weighted_bipartite_edge>& edge);

        void remove_edge(uint32_t right_vertex_id);

    private:
        uint32_t left_vertex_id;
        shared_ptr<unordered_map<uint32_t,shared_ptr<weighted_bipartite_edge>>> edge_map;
    };
}
