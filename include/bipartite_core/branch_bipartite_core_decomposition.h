/******************************************************************************************************************
* @copyright  : Boost Software License
* @file       : branch_bipartite_core_decomposition.h
* @details    : an branch bipartite core decomposition
* @version    : 1.0
* @date       : 2022/01/02
******************************************************************************************************************/

#pragma once
#include "bipartite_core/left_vertex_index.h"
#include "bipartite_core/right_vertex_index.h"
#include "bipartite_core/bipartite_core_degree_index.h"

namespace scnu{
    class branch_bipartite_core_decomposition {
    public:

        static uint32_t decompose(const shared_ptr<abstract_bipartite_graph> &B,
                                  const shared_ptr<unordered_map<uint32_t, shared_ptr<left_vertex_index>>> &left_index_map,
                                  const shared_ptr<unordered_map<uint32_t, shared_ptr<right_vertex_index>>> &right_index_map);


        static uint32_t decompose(const shared_ptr<abstract_bipartite_graph> &B,
                                  const shared_ptr<unordered_map<uint32_t, shared_ptr<left_vertex_index>>> &left_index_map,
                                  const shared_ptr<unordered_map<uint32_t, shared_ptr<right_vertex_index>>> &right_index_map,
                                  const shared_ptr<thread_pool>& pool);


        static void find_core(const shared_ptr<abstract_bipartite_graph> &B,
                              const shared_ptr<unordered_map<uint32_t, uint32_t>> &left_degree_map,
                              const shared_ptr<unordered_map<uint32_t, uint32_t>> &right_degree_map,
                              uint32_t k);


        static uint32_t find_left_core(const shared_ptr<abstract_bipartite_graph> &B,
                                       const shared_ptr<unordered_map<uint32_t, uint32_t>> &left_map,
                                       const shared_ptr<unordered_map<uint32_t, uint32_t>> &right_map,
                                       const shared_ptr<unordered_map<uint32_t, shared_ptr<left_vertex_index>>> &left_index_map,
                                       const shared_ptr<unordered_map<uint32_t, shared_ptr<right_vertex_index>>> &right_index_map,
                                       uint32_t k);

        static uint32_t find_right_core(const shared_ptr<abstract_bipartite_graph> &B,
                                        const shared_ptr<unordered_map<uint32_t, uint32_t>> &left_map,
                                        const shared_ptr<unordered_map<uint32_t, uint32_t>> &right_map,
                                        const shared_ptr<unordered_map<uint32_t, shared_ptr<left_vertex_index>>> &left_index_map,
                                        const shared_ptr<unordered_map<uint32_t, shared_ptr<right_vertex_index>>> &right_index_map,
                                        uint32_t k);



        static uint32_t find_left_core(const shared_ptr<abstract_bipartite_graph> &B,
                                       const shared_ptr<unordered_map<uint32_t, shared_ptr<mutex>>> &left_mutex_map,
                                       const shared_ptr<unordered_map<uint32_t, shared_ptr<mutex>>> &right_mutex_map,
                                    const shared_ptr<unordered_map<uint32_t, uint32_t>> &left_map,
                                    const shared_ptr<unordered_map<uint32_t, uint32_t>> &right_map,
                                    const shared_ptr<unordered_map<uint32_t, shared_ptr<left_vertex_index>>> &left_index_map,
                                    const shared_ptr<unordered_map<uint32_t, shared_ptr<right_vertex_index>>> &right_index_map,
                                    uint32_t k);




        static uint32_t find_right_core(const shared_ptr<abstract_bipartite_graph> &B,
                                        const shared_ptr<unordered_map<uint32_t, shared_ptr<mutex>>> &left_mutex_map,
                                        const shared_ptr<unordered_map<uint32_t, shared_ptr<mutex>>> &right_mutex_map,
                                        const shared_ptr<unordered_map<uint32_t, uint32_t>> &left_map,
                                        const shared_ptr<unordered_map<uint32_t, uint32_t>> &right_map,
                                        const shared_ptr<unordered_map<uint32_t, shared_ptr<left_vertex_index>>> &left_index_map,
                                        const shared_ptr<unordered_map<uint32_t, shared_ptr<right_vertex_index>>> &right_index_map,
                                       uint32_t k);


    };
}
