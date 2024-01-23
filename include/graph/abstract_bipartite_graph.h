/******************************************************************************************************************
* @copyright  : Boost Software License
* @file       : abstract_bipartite_graph.h
* @brief      : An abstract bipartite graph structure
* @version    : 1.1
* @date       : 2020/10/15
******************************************************************************************************************/

#pragma once
#include "graph/abstract_left_vertex.h"
#include "graph/abstract_right_vertex.h"

namespace scnu {
    /**
     * @details an abstract bipartite graph with a left vertex map and a right vertex map
     */
    class abstract_bipartite_graph {
    public:
        abstract_bipartite_graph();

        explicit abstract_bipartite_graph(const shared_ptr<abstract_bipartite_graph> &graph);

        abstract_bipartite_graph(
                const shared_ptr<unordered_map<uint32_t, shared_ptr<abstract_left_vertex>>> &other_left_vertex_map,
                const shared_ptr<unordered_map<uint32_t, shared_ptr<abstract_right_vertex>>> &other_right_vertex_map);

        ~abstract_bipartite_graph() = default;

        double get_density();

        shared_ptr<abstract_bipartite_edge> get_edge(uint32_t left_vertex_id, uint32_t right_vertex_id);

        shared_ptr<unordered_set<shared_ptr<abstract_bipartite_edge>>> get_edge_set();

        shared_ptr<abstract_left_vertex> get_left_vertex(uint32_t left_vertex_id);

        shared_ptr<unordered_map<uint32_t, shared_ptr<abstract_left_vertex>>> get_left_vertex_map();

        double get_average_left_vertex_degree();

        double get_average_right_vertex_degree();

        double get_average_vertex_degree();

        uint32_t get_maximal_left_vertex_degree();

        uint32_t get_maximal_right_vertex_degree();



        uint32_t get_edge_number();

        uint32_t get_left_vertex_number();

        uint32_t get_right_vertex_number();

        shared_ptr<abstract_right_vertex> get_right_vertex(uint32_t right_vertex_id);

        shared_ptr<unordered_map<uint32_t, shared_ptr<abstract_right_vertex>>> get_right_vertex_map();

        void insert_edge(const shared_ptr<abstract_bipartite_edge> &e);

        /**
        * @details insert a collection of edges
        * @remarks container type must be sequence container, such as vector,list, set, unordered_set, etc.
        * @tparam container_type
        * @param edge_container
        */
        template<typename container_type>
        void insert_edge_collection(const shared_ptr<container_type>& edge_container)
        {
            for(const auto &e:*edge_container)
            {
                insert_edge(e);
            }
        }

        void insert_left_vertex(uint32_t left_vertex_id);

        void insert_right_vertex(uint32_t right_vertex_id);

        bool empty();

        void remove_edge(const shared_ptr<abstract_bipartite_edge> &edge);

        void remove_edge(const shared_ptr<abstract_bipartite_edge> &edge,
                         const shared_ptr<unordered_set<uint32_t>> &isolated_left_vertex_set,
                         const shared_ptr<unordered_set<uint32_t>> &isolated_right_vertex_set);

        /**
        * @details remove a collection of edges
        * @remarks container type must be sequence container, such as vector,list, set, unordered_set, etc.
        * @tparam container_type
        * @param edge_container
        */
        template<typename container_type>
        void remove_edge_collection(const shared_ptr<container_type>& edge_container)
        {
            for(const auto &e:*edge_container)
            {
                remove_edge(e);
            }
        }

        /**
        * @details remove a collection of edges
        * @remarks container type must be sequence container, such as vector,list, set, unordered_set, etc.
        * @tparam container_type
        * @param edge_container
        */
        template<typename container_type>
        void remove_edge_collection(const shared_ptr<container_type>& edge_container, const shared_ptr<unordered_set<uint32_t>> &isolated_left_vertex_set,
                                    const shared_ptr<unordered_set<uint32_t>> &isolated_right_vertex_set)
        {
            for(const auto &e:*edge_container)
            {
                remove_edge(e, isolated_left_vertex_set, isolated_right_vertex_set);
            }
        }

        void remove_right_vertex(uint32_t right_vertex_id);

        void remove_left_vertex(uint32_t left_vertex_id);

    private:
        shared_ptr<unordered_map<uint32_t, shared_ptr<abstract_left_vertex>>> left_vertex_map;
        shared_ptr<unordered_map<uint32_t, shared_ptr<abstract_right_vertex>>> right_vertex_map;
    };
}


