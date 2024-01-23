/******************************************************************************************************************
* @copyright  : Boost Software License
* @file       : weighted_vertex.h
* @brief      : a vertex associated with weighted edges
* @version    : 1.0
* @date       : 2021/08/03
******************************************************************************************************************/

#pragma once
#include "weighted_edge.h"

namespace scnu{
    class weighted_vertex {
    public:
        explicit weighted_vertex(uint32_t vertex_id);

        explicit weighted_vertex(const shared_ptr<weighted_vertex>& other_vertex);

        weighted_vertex(uint32_t vertex_id,
        const shared_ptr<unordered_map<uint32_t, shared_ptr<weighted_edge>>>& edge_map);

        bool exist_edge(uint32_t neighbor_vertex_id);

        shared_ptr<unordered_map<uint32_t , shared_ptr<weighted_edge>>> get_neighbor_map();

        uint32_t get_degree();

        shared_ptr<weighted_edge> get_edge(uint32_t neighbor_vertex_id);

        [[nodiscard]] uint32_t get_vertex_id() const;

        void insert_edge(uint32_t neighbor_vertex_id, const shared_ptr<weighted_edge> &edge);

        void remove_edge(uint32_t neighbor_vertex_id);

    private:
        uint32_t vertex_id;
        shared_ptr<unordered_map<uint32_t, shared_ptr<weighted_edge>>> neighbor_vertex_map;
    };
}

