/******************************************************************************************************************
* @copyright  : Boost Software License
* @file       : abstract_vertex.cpp
* @brief      : An abstract graph structure
* @version    : 1.1
* @date       : 2020/10/15
******************************************************************************************************************/

#pragma once
#include "graph/abstract_vertex.h"

namespace scnu
{
    /**
     * @details an simple graph with a simple vertex map
     */
    class abstract_graph {
    public:
        explicit abstract_graph();

        explicit abstract_graph(const shared_ptr<unordered_map<uint32_t,shared_ptr<abstract_vertex>>>& other_vertex_map);

        explicit abstract_graph(const shared_ptr<abstract_graph>& other_graph);

        ~abstract_graph() = default;

        bool empty();

        shared_ptr<abstract_edge> get_edge(uint32_t source_vertex_id, uint32_t destination_vertex_id);

        shared_ptr<unordered_set<shared_ptr<abstract_edge>>> get_edge_set();

        uint32_t get_maximal_degree();

        uint32_t get_edge_number();

        uint32_t get_vertex_number();

        shared_ptr<abstract_vertex> get_vertex(uint32_t vertexId);

        shared_ptr<unordered_map<uint32_t,shared_ptr<abstract_vertex>>> get_vertex_map();

        void insert_edge(const shared_ptr<abstract_edge>& edge);

        /**
         * @details insert a collection of edges into this graph
         * @details container type must be sequence container, such as vector,list, set, unordered_set, etc.
         * @tparam container_type 
         * @param edge_container 
         */
        template<typename container_type>
        void insert_edge_collection(const container_type& edge_container)
        {
            for(const auto& e:*edge_container)
            {
                insert_edge(e);
            }
        }

        void insert_vertex(uint32_t vertex_id);

        shared_ptr<unordered_set<uint32_t>> remove_edge(const shared_ptr<abstract_edge>& edge);

        shared_ptr<unordered_set<uint32_t>> remove_edge(const shared_ptr<abstract_edge> &edge,
                                                        const shared_ptr<unordered_set<uint32_t>> &isolated_vertex_set);

        /**
         * @details remove a collection of edges into this graph
         * @details container type must be sequence container, such as vector,list, set, unordered_set, etc.
         * @tparam container_type 
         * @param edge_container 
         */
        template<typename container_type>
        shared_ptr<unordered_set<uint32_t>> remove_edge_collection(const container_type& edge_container)
        {
            auto isolated_vertex_set = make_shared<unordered_set<uint32_t>>();
            for(const auto& e:*edge_container)
            {
                auto source_vertex_id = e->get_source_vertex_id();
                auto destination_vertex_id = e->get_destination_vertex_id();

                if (get_edge(source_vertex_id, destination_vertex_id)) {
                    auto source_vertex = get_vertex(source_vertex_id);
                    source_vertex->remove_edge(destination_vertex_id);

                    if (source_vertex->get_degree() == 0) {
                        remove_vertex(source_vertex_id);
                        isolated_vertex_set->insert(source_vertex_id);
                    }

                    auto destination_vertex = get_vertex(destination_vertex_id);
                    destination_vertex->remove_edge(source_vertex_id);

                    if (destination_vertex->get_degree() == 0) {
                        remove_vertex(destination_vertex_id);
                        isolated_vertex_set->insert(destination_vertex_id);
                    }
                }
            }
            return isolated_vertex_set;
        }

        void remove_vertex(uint32_t vertex_id);

    private:
        shared_ptr<unordered_map<uint32_t,shared_ptr<abstract_vertex>>> vertex_map;
    };
}


