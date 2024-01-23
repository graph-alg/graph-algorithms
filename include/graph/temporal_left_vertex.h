/******************************************************************************************************************
* @copyright  : Boost Software License
* @file       : temporal_left_vertex.h
* @brief      : a left vertex with temporal edges
* @version    : 1.0
* @date       : 2022/06/17
******************************************************************************************************************/

#pragma once
#include "graph/temporal_bipartite_edge.h"

namespace  scnu
{
    class temporal_left_vertex
    {
    public:
        explicit temporal_left_vertex(uint32_t other_l);

        explicit temporal_left_vertex(const shared_ptr<temporal_left_vertex>& l_vertex);

        virtual ~temporal_left_vertex() = default;

        [[nodiscard]] uint32_t get_left_vertex_id() const;

        shared_ptr<unordered_set<shared_ptr<temporal_bipartite_edge>>> get_edge_set(uint32_t r);

        shared_ptr<unordered_map<uint32_t,shared_ptr<unordered_set<shared_ptr<temporal_bipartite_edge>>>>> get_edge_map();

        uint32_t get_neighbor_size();

        void insert_edge(uint32_t r, const shared_ptr<temporal_bipartite_edge>& e);

        void insert_edge_set(uint32_t r, const shared_ptr<unordered_set<shared_ptr<temporal_bipartite_edge>>>& e);

        void remove_edge(uint32_t r, const shared_ptr<temporal_bipartite_edge>& e);

        void remove_edge_set(uint32_t r);

        void remove_edge_set(uint32_t r,
                             const shared_ptr<unordered_set<shared_ptr<temporal_bipartite_edge>>> &edge_set);

    private:
        uint32_t left_vertex_id;
        shared_ptr<unordered_map<uint32_t,shared_ptr<unordered_set<shared_ptr<temporal_bipartite_edge>>>>> edge_map;
    };
}

