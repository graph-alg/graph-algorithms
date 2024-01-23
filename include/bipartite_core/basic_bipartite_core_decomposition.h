/******************************************************************************************************************
* @copyright  : Boost Software License
* @file       : bipartite_core_utility.hpp
* @details    : common head files for bipartite core
* @version    : 1.0
* @date       : 2021/02/01
******************************************************************************************************************/

#pragma once
#include "bipartite_core/bipartite_core.h"
#include "bipartite_core/bipartite_core_left_store_index.h"
#include "bipartite_core/bipartite_core_right_store_index.h"

namespace scnu{
    class basic_bipartite_core_decomposition{
    public:
        static shared_ptr<unordered_set<pair<uint32_t,uint32_t>,hash_pair,equal_pair>>
        decompose(const shared_ptr<abstract_bipartite_graph>& G,
                  const shared_ptr<unordered_map<uint32_t, shared_ptr<bipartite_core_left_store_index>>>& left_index_map,
                  const shared_ptr<unordered_map<uint32_t, shared_ptr<bipartite_core_right_store_index>>>& right_index_map);

        static shared_ptr<unordered_set<pair<uint32_t,uint32_t>,hash_pair,equal_pair>>
        decompose(const shared_ptr<abstract_bipartite_graph>& G,
                  const shared_ptr<unordered_map<uint32_t, shared_ptr<bipartite_core_left_store_index>>>& left_index_map,
                  const shared_ptr<unordered_map<uint32_t, shared_ptr<bipartite_core_right_store_index>>>& right_index_map,
                  uint32_t thread_number);


        static void find_left_core(const shared_ptr<abstract_bipartite_graph> &G,
                                   const shared_ptr<unordered_map<uint32_t, uint32_t>> &left_vertex_degree_map,
                                   const shared_ptr<unordered_map<uint32_t, uint32_t>> &right_vertex_degree_map,
                                   uint32_t i,
                                   uint32_t j);

        static void find_right_core(const shared_ptr<abstract_bipartite_graph> &G,
                                    const shared_ptr<unordered_map<uint32_t,uint32_t>>& left_vertex_degree_map,
                                    const shared_ptr<unordered_map<uint32_t,uint32_t>>& right_vertex_degree_map,
                                    uint32_t i,
                                    uint32_t j);


        static void remove_unsatisfied_vertices(const shared_ptr<abstract_bipartite_graph>& G,
                                                const shared_ptr<unordered_map<uint32_t,uint32_t>>& left_vertex_degree_map,
                                                const shared_ptr<unordered_map<uint32_t,uint32_t>>& right_vertex_degree_map,
                                                const shared_ptr<unordered_set<uint32_t>>& l_set,
                                                const shared_ptr<unordered_set<uint32_t>>& r_set,
                                                uint32_t i,
                                                uint32_t j);
    };
}
