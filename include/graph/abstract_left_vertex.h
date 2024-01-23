/******************************************************************************************************************
* @copyright  : Boost Software License
* @file       : AbstractVertex.cpp
* @brief      : An abstract left vertex structure
* @version    : 1.1
* @date       : 2020/10/15
******************************************************************************************************************/

#pragma once
#include "graph/abstract_bipartite_edge.h"

namespace  scnu
{
    /**
     * @details A left vertex class for a bipartite graph
     */
    class abstract_left_vertex
    {
    public:
        explicit abstract_left_vertex(uint32_t other_left_vertex_id);

        explicit abstract_left_vertex(const shared_ptr<abstract_left_vertex>& left_vertex);

        virtual ~abstract_left_vertex() = default;

        [[nodiscard]] uint32_t get_left_vertex_id() const;

        shared_ptr<abstract_bipartite_edge> get_edge(uint32_t right_vertex_id);

        shared_ptr<unordered_map<uint32_t,shared_ptr<abstract_bipartite_edge>>> get_edge_map();

        uint32_t get_degree();

        void insert_edge(uint32_t right_vertex_id, const shared_ptr<abstract_bipartite_edge>& edge);

        void remove_edge(uint32_t right_vertex_id);

    private:
        uint32_t left_vertex_id;
        shared_ptr<unordered_map<uint32_t,shared_ptr<abstract_bipartite_edge>>> edge_map;
    };
}




