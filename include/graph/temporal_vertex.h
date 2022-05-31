/******************************************************************************************************************
* @copyright  : Boost Software License
* @file       : temporal_graph.h
* @brief      : A graph  whose edge containing time information
* @version    : 1.0
* @date       : 2021/08/03
******************************************************************************************************************/

#pragma once
#include "temporal_edge.h"

namespace scnu{
    class temporal_vertex {
    public:
        explicit temporal_vertex(uint32_t vertex_id);

        explicit temporal_vertex(const shared_ptr<temporal_vertex>& other_vertex);

        temporal_vertex(uint32_t vertex_id,
                       const shared_ptr<unordered_map<uint32_t , shared_ptr<unordered_set<shared_ptr<temporal_edge>>>>>& edge_map);

        bool exist_edge(uint32_t neighbor_vertex_id, const shared_ptr<temporal_edge> &edge);

        shared_ptr<unordered_map<uint32_t , shared_ptr<unordered_set<shared_ptr<temporal_edge>>>>> get_neighbor_map();

        uint32_t get_neighbor_size();

        shared_ptr<unordered_set<shared_ptr<temporal_edge>>> get_temporal_edge_set(uint32_t neighbor_vertex_id);

        uint32_t get_temporal_edge_size(uint32_t neighbor_vertex_id);

        uint32_t get_vertex_id() const;

        void insert_edge(uint32_t neighbor_vertex_id, const shared_ptr<temporal_edge> &edge);

        void insert_edge(uint32_t neighbor_vertex_id, const shared_ptr<unordered_set<shared_ptr<temporal_edge>>> &temporal_edge_set);

        void insert_neighbor_vertex(uint32_t neighbor_vertex_id);

        void insert_neighbor_vertex(uint32_t neighbor_vertex_id, const shared_ptr<unordered_set<shared_ptr<temporal_edge>>> &temporal_edge_set);

        void remove_edge(uint32_t neighbor_vertex_id, const shared_ptr<temporal_edge>& edge);

        void remove_edge(uint32_t neighbor_vertex_id, const shared_ptr<unordered_set<shared_ptr<temporal_edge>>> &sub_temporal_edge_set);

        void remove_neighbor_vertex(uint32_t neighbor_vertex_id);

    private:
        uint32_t vertex_id;
        shared_ptr<unordered_map<uint32_t, shared_ptr<unordered_set<shared_ptr<temporal_edge>>>>> neighbor_vertex_map;
    };
}

