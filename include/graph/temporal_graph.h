/******************************************************************************************************************
* @copyright  : Boost Software License
* @file       : temporal_graph.h
* @brief      : A graph  whose edge containing time information
* @version    : 1.0
* @date       : 2021/08/03
******************************************************************************************************************/

#pragma  once
#include "temporal_vertex.h"
#include "graph_utility.h"

namespace scnu{
    class temporal_graph {
    public:
        temporal_graph();

        explicit temporal_graph(const shared_ptr<unordered_map<uint32_t , shared_ptr<temporal_vertex>>>& vertex_map);

        explicit temporal_graph(const shared_ptr<temporal_graph>& G);

        ~temporal_graph();

        bool exist_edge(const shared_ptr<temporal_edge>& edge);

        bool exist_edge(uint32_t source_vertex_id, uint32_t destination_vertex_id);

        bool exist_vertex(uint32_t vertex_id);

        uint32_t get_edge_size();

        uint32_t get_maximal_neighbor_vertex_size();

        uint32_t get_maximal_parallel_edge_size();

        uint32_t get_vertex_size();

        shared_ptr<unordered_set<shared_ptr<temporal_edge>>> get_edge_set();

        shared_ptr<unordered_set<shared_ptr<temporal_edge>>> get_edge_set(uint32_t source_vertex_id, uint32_t destination_vertex_id);

        shared_ptr<temporal_vertex> get_vertex(uint32_t vertex_id);

        shared_ptr<unordered_map<uint32_t, shared_ptr<temporal_vertex>>> get_vertex_map();

        void insert_edge(const shared_ptr<temporal_edge>& edge);

        void insert_edge(const shared_ptr<unordered_set<shared_ptr<temporal_edge>>> & edge_set);

        void insert_vertex(uint32_t vertex_id);

        void insert_vertex(const shared_ptr<temporal_vertex> &vertex);

        void remove_edge(const shared_ptr<temporal_edge>& edge);

        void remove_edges(uint32_t source_vertex_id, uint32_t destination_vertex_id);

        void remove_vertex(uint32_t vertex_id);

    private:
        shared_ptr<unordered_map<uint32_t , shared_ptr<temporal_vertex>>> vertex_map;
    };
}





