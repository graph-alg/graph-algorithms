/******************************************************************************************************************
* @copyright  : Boost Software License
* @file       : weighted_graph.h
* @brief      : A graph whose edge containing weight information
* @version    : 1.0
* @date       : 2022/6/15
******************************************************************************************************************/
#pragma once
#include "weighted_vertex.h"

namespace scnu{
    class weighted_graph {
    public:
        weighted_graph();

        explicit weighted_graph(const shared_ptr<unordered_map<uint32_t , shared_ptr<weighted_vertex>>>& vertex_map);

        explicit weighted_graph(const shared_ptr<weighted_graph>& G);

        ~weighted_graph();

        bool empty();

        bool exist_edge(const shared_ptr<weighted_edge>& edge);

        bool exist_edge(uint32_t source_vertex_id, uint32_t destination_vertex_id);

        bool exist_vertex(uint32_t vertex_id);

        uint32_t get_edge_size();

        uint32_t get_vertex_size();

        shared_ptr<weighted_vertex> get_vertex(uint32_t vertex_id);

        shared_ptr<weighted_edge> get_edge(uint32_t source_vertex_id, uint32_t destination_vertex_id);

        shared_ptr<unordered_set<shared_ptr<weighted_edge>>> get_edge_set();

        shared_ptr<unordered_map<uint32_t, shared_ptr<weighted_vertex>>> get_vertex_map();

        void insert_edge(const shared_ptr<weighted_edge>& edge);

        void insert_edge(const shared_ptr<unordered_set<shared_ptr<weighted_edge>>> & edge_set);

        void insert_vertex(uint32_t vertex_id);

        void insert_vertex(const shared_ptr<weighted_vertex> &vertex);

        void remove_edge(const shared_ptr<weighted_edge>& edge);

        void remove_edges(uint32_t source_vertex_id, uint32_t destination_vertex_id);

        void remove_vertex(uint32_t vertex_id);

    private:
        shared_ptr<unordered_map<uint32_t , shared_ptr<weighted_vertex>>> vertex_map;
    };
}


