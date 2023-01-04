/******************************************************************************************************************
* @copyright  : Boost Software License
* @file       : share_bipartite_core_decomposition.h
* @details    : An short core decomposition algorithm
* @version    : 1.0
* @date       : 2021/02/01
******************************************************************************************************************/

#pragma once
#include "bipartite_core/left_vertex_index.h"
#include "bipartite_core/right_vertex_index.h"

namespace scnu{
    class share_bipartite_core_decomposition {
    public:
        static uint32_t decompose(const shared_ptr<abstract_bipartite_graph> &B,
                                  const shared_ptr<unordered_map<uint32_t, shared_ptr<left_vertex_index>>> &left_index_map,
                                  const shared_ptr<unordered_map<uint32_t, shared_ptr<right_vertex_index>>> &right_index_map);

        static uint32_t decompose(const shared_ptr<abstract_bipartite_graph> &B,
                                  const shared_ptr<unordered_map<uint32_t, shared_ptr<left_vertex_index>>> &left_index_map,
                                  const shared_ptr<unordered_map<uint32_t, shared_ptr<right_vertex_index>>> &right_index_map,
                                  const shared_ptr<thread_pool>& pool);

    private:
        static void compute_alpha_core(const shared_ptr<abstract_bipartite_graph> &B,
                                       const shared_ptr<unordered_map<uint32_t, shared_ptr<left_vertex_index>>> &left_index_map,
                                       const shared_ptr<unordered_map<uint32_t, shared_ptr<right_vertex_index>>> &right_index_map,
                                       uint32_t alpha);

        static void compute_alpha_core(const shared_ptr<abstract_bipartite_graph>& B,
                                       const shared_ptr<unordered_map<uint32_t, shared_ptr<mutex>>>& left_mutex_map,
                                       const shared_ptr<unordered_map<uint32_t, shared_ptr<mutex>>>& right_mutex_map,
                                       const shared_ptr<unordered_map<uint32_t, shared_ptr<left_vertex_index>>> &left_index_map,
                                       const shared_ptr<unordered_map<uint32_t , shared_ptr<right_vertex_index>>> &right_index_map,
                                       uint32_t alpha);

        static void compute_beta_core(const shared_ptr<abstract_bipartite_graph> &B,
                                      const shared_ptr<unordered_map<uint32_t, shared_ptr<left_vertex_index>>> &left_index_map,
                                      const shared_ptr<unordered_map<uint32_t, shared_ptr<right_vertex_index>>> &right_index_map,
                                      uint32_t beta);

        static void compute_beta_core(const shared_ptr<abstract_bipartite_graph>& B,
                                      const shared_ptr<unordered_map<uint32_t, shared_ptr<mutex>>>& left_mutex_map,
                                      const shared_ptr<unordered_map<uint32_t, shared_ptr<mutex>>>& right_mutex_map,
                                       const shared_ptr<unordered_map<uint32_t , shared_ptr<left_vertex_index>>>& left_index_map,
                                       const shared_ptr<unordered_map<uint32_t , shared_ptr<right_vertex_index>>>& right_index_map,
                                      uint32_t beta);

        static uint32_t find_core(const shared_ptr<abstract_bipartite_graph> &B);

    };
}
