/******************************************************************************************************************
* @copyright  : Boost Software License
* @file       : abstract_vertex.cpp
* @brief      : An abstract vertex structure
* @version    : 1.1
* @date       : 2020/10/15
******************************************************************************************************************/
#pragma once
#include "graph/abstract_edge.h"

namespace scnu
{
    /**
     * @details An abstract vertex class
     */
    class abstract_vertex {
    public:
        explicit abstract_vertex(uint32_t other_vertex_id);

        explicit abstract_vertex(const shared_ptr<abstract_vertex> &other_vertex);

        virtual ~abstract_vertex() = default;

        [[nodiscard]] uint32_t get_vertex_id() const;

        shared_ptr<abstract_edge> get_edge(uint32_t adjacent_vertex_id);

        shared_ptr<unordered_map<uint32_t ,shared_ptr<abstract_edge>>> get_edge_map();

        uint32_t get_degree();

        uint32_t get_degree(const shared_ptr<unordered_set<shared_ptr<abstract_edge>>>& removed_edge_set);

        void insert_edge(uint32_t adjacent_vertex_id, const shared_ptr<abstract_edge>& adjacent_edge);

        void remove_edge(uint32_t adjacent_vertex_id);

    private:
        uint32_t vertex_id;
        shared_ptr<unordered_map<uint32_t,shared_ptr<abstract_edge>>> edge_map;
    };
}



