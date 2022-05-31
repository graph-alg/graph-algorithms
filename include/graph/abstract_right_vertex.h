/******************************************************************************************************************
* @copyright  : Boost Software License
* @file       : AbstractVertex.cpp
* @brief      : An abstract right vertex structure
* @version    : 1.1
* @date       : 2020/10/15
******************************************************************************************************************/
#pragma once
#include "graph/abstract_bipartite_edge.h"

namespace scnu
{
    /**
     * @brief An abstract right vertex class for a abstract bipartite graph
     */
    class abstract_right_vertex
    {
    public:
        explicit abstract_right_vertex(uint32_t other_right_vertex_id);

        explicit abstract_right_vertex(const shared_ptr<abstract_right_vertex>& other_right_vertex);

        virtual ~abstract_right_vertex() = default;

        shared_ptr<abstract_bipartite_edge> get_edge(uint32_t other_right_vertex_id);

        shared_ptr<unordered_map<uint32_t,shared_ptr<abstract_bipartite_edge>>> get_edge_map();

        uint32_t get_degree();

        [[nodiscard]] uint32_t get_right_vertex_id() const;

        void insert_edge(uint32_t other_left_vertex_id, const shared_ptr<abstract_bipartite_edge>& edge);

        void remove_edge(uint32_t other_left_vertex_id);

    private:
        uint32_t right_vertex_id;
        shared_ptr<unordered_map<uint32_t,shared_ptr<abstract_bipartite_edge>>> edge_map;
    };
}


