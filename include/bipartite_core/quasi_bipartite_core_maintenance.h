/******************************************************************************************************************
* @copyright  : Boost Software License
* @file       : share_bipartite_core_decomposition.h
* @details    : An short core decomposition algorithm
* @version    : 1.0
* @date       : 2021/02/01
******************************************************************************************************************/

#pragma once
#include "bipartite_core/bipartite_core_utility.h"
#include "bipartite_core/bipartite_core.h"
#include "bipartite_core/bipartite_core_compare.h"

namespace scnu{
    class quasi_bipartite_core_maintenance {
    public:

        static void init(const shared_ptr<abstract_bipartite_graph> &B,
                         const shared_ptr<unordered_map<uint32_t, shared_ptr<bipartite_core_left_store_index>>> &new_left_index_map,
                         const shared_ptr<unordered_map<uint32_t, shared_ptr<bipartite_core_right_store_index>>> &new_right_index_map);

        static void init(const shared_ptr<abstract_bipartite_graph> &B,
                         const shared_ptr<unordered_map<uint32_t, shared_ptr<bipartite_core_left_store_index>>> &new_left_index_map,
                         const shared_ptr<unordered_map<uint32_t, shared_ptr<bipartite_core_right_store_index>>> &new_right_index_map,
                         const shared_ptr<thread_pool>& pool);


        static void insert(const shared_ptr<abstract_bipartite_graph> &B,
                           const shared_ptr<unordered_set<shared_ptr<abstract_bipartite_edge>>> &edge_set,
                           const shared_ptr<unordered_map<uint32_t,shared_ptr<bipartite_core_left_store_index>>> &left_index_map,
                           const shared_ptr<unordered_map<uint32_t,shared_ptr<bipartite_core_right_store_index>>> &right_index_map,
                           const shared_ptr<unordered_map<uint32_t,shared_ptr<bipartite_core_left_store_index>>> &new_left_index_map,
                           const shared_ptr<unordered_map<uint32_t,shared_ptr<bipartite_core_right_store_index>>> &new_right_index_map);

        static void insert(const shared_ptr<abstract_bipartite_graph> &B,
                           const shared_ptr<unordered_set<shared_ptr<abstract_bipartite_edge>>> &edge_set,
                           const shared_ptr<unordered_map<uint32_t,shared_ptr<bipartite_core_left_store_index>>> &left_index_map,
                           const shared_ptr<unordered_map<uint32_t,shared_ptr<bipartite_core_right_store_index>>> &right_index_map,
                           const shared_ptr<unordered_map<uint32_t,shared_ptr<bipartite_core_left_store_index>>> &new_left_index_map,
                           const shared_ptr<unordered_map<uint32_t,shared_ptr<bipartite_core_right_store_index>>> &new_right_index_map,
                           const shared_ptr<thread_pool>& pool);


        static void remove(const shared_ptr<abstract_bipartite_graph> &B,
                           const shared_ptr<unordered_set<shared_ptr<abstract_bipartite_edge>>> &edge_set,
                           const shared_ptr<unordered_map<uint32_t, shared_ptr<bipartite_core_left_store_index>>> &left_index_map,
                           const shared_ptr<unordered_map<uint32_t, shared_ptr<bipartite_core_right_store_index>>> &right_index_map,
                           const shared_ptr<unordered_map<uint32_t,shared_ptr<bipartite_core_left_store_index>>>& new_left_index_map,
                           const shared_ptr<unordered_map<uint32_t,shared_ptr<bipartite_core_right_store_index>>>& new_right_index_map);

        static void remove(const shared_ptr<abstract_bipartite_graph> &B,
                           const shared_ptr<unordered_set<shared_ptr<abstract_bipartite_edge>>> &edge_set,
                           const shared_ptr<unordered_map<uint32_t, shared_ptr<bipartite_core_left_store_index>>> &left_index_map,
                           const shared_ptr<unordered_map<uint32_t, shared_ptr<bipartite_core_right_store_index>>> &right_index_map,
                           const shared_ptr<unordered_map<uint32_t,shared_ptr<bipartite_core_left_store_index>>>& new_left_index_map,
                           const shared_ptr<unordered_map<uint32_t,shared_ptr<bipartite_core_right_store_index>>>& new_right_index_map,
                           const shared_ptr<thread_pool>& pool);

    private:
        static uint32_t compute_left_vertex_core_degree(const shared_ptr<abstract_bipartite_graph> &G,
                const shared_ptr<unordered_map<uint32_t,shared_ptr<bipartite_core_right_store_index>>>& right_index_map,
                uint32_t l,
                uint32_t i,
                uint32_t j);


        static uint32_t compute_right_vertex_core_degree(
                const shared_ptr<abstract_bipartite_graph> &G,
                const shared_ptr<unordered_map<uint32_t,shared_ptr<bipartite_core_left_store_index>>> & left_index_map,
                uint32_t r,
                uint32_t i,
                uint32_t j);

        static bool left_candidate_graph(const shared_ptr<abstract_bipartite_graph> &B,
                                    const shared_ptr<unordered_set<shared_ptr<abstract_bipartite_edge>>> &edge_set,
                                    const shared_ptr<unordered_map<uint32_t,shared_ptr<bipartite_core_left_store_index>>>& left_index_map,
                                    const shared_ptr<unordered_map<uint32_t,shared_ptr<bipartite_core_right_store_index>>>& right_index_map,
                                    uint32_t i,
                                    uint32_t j,
                                    const shared_ptr<unordered_map<uint32_t, uint32_t>> & candidate_l_map,
                                    const shared_ptr<unordered_map<uint32_t, uint32_t>> & candidate_r_map);

        static bool right_candidate_graph(const shared_ptr<abstract_bipartite_graph> &B,
                                         const shared_ptr<unordered_set<shared_ptr<abstract_bipartite_edge>>> &edge_set,
                                         const shared_ptr<unordered_map<uint32_t,shared_ptr<bipartite_core_left_store_index>>>& left_index_map,
                                         const shared_ptr<unordered_map<uint32_t,shared_ptr<bipartite_core_right_store_index>>>& right_index_map,
                                         uint32_t i,
                                         uint32_t j,
                                         const shared_ptr<unordered_map<uint32_t, uint32_t>> & candidate_l_map,
                                         const shared_ptr<unordered_map<uint32_t, uint32_t>> & candidate_r_map);

        static void find_quasi_core(const shared_ptr<abstract_bipartite_graph> &G,
                                    const shared_ptr<unordered_set<shared_ptr<abstract_bipartite_edge>>>& edge_set,
                                    const shared_ptr<unordered_map<uint32_t, uint32_t>> &l_map,
                                    const shared_ptr<unordered_map<uint32_t, uint32_t>> &r_map,
                                    uint32_t i,
                                    uint32_t j);


        static void get_quasi_cores(const shared_ptr<abstract_bipartite_graph> &B,
                                    const shared_ptr<unordered_set<shared_ptr<abstract_bipartite_edge>>>& edge_set,
                                    const shared_ptr<unordered_map<pair<uint32_t, uint32_t>, shared_ptr<unordered_map<uint32_t, uint32_t>>, hash_pair, equal_pair>> &quasi_core_left_vertex_degree_map,
                                    const shared_ptr<unordered_map<pair<uint32_t, uint32_t>, shared_ptr<unordered_map<uint32_t, uint32_t>>, hash_pair, equal_pair>> &quasi_core_right_vertex_degree_map,
                                    uint64_t i,
                                    uint64_t j,
                                    const shared_ptr<unordered_map<pair<uint32_t, uint32_t>, shared_ptr<unordered_set<shared_ptr<abstract_bipartite_edge>>>, hash_pair, equal_pair>> &next_quasi_edge_map);


        static void remove_left_vertex(const shared_ptr<abstract_bipartite_graph> &B,
                                        const shared_ptr<unordered_map<uint32_t, uint32_t>> &candidate_l_map,
                                        const shared_ptr<unordered_map<uint32_t, uint32_t>> &candidate_r_map,
                                        const shared_ptr<unordered_set<uint32_t>> &evicted_l_set,
                                        const shared_ptr<unordered_set<uint32_t>> &evicted_r_set,
                                        uint32_t l,
                                        uint32_t i,
                                        uint32_t j);

        static void remove_right_vertex(const shared_ptr<abstract_bipartite_graph> &B,
                                        const shared_ptr<unordered_map<uint32_t, uint32_t>> &candidate_l_map,
                                        const shared_ptr<unordered_map<uint32_t, uint32_t>> &candidate_r_map,
                                        const shared_ptr<unordered_set<uint32_t>> &evicted_l_set,
                                        const shared_ptr<unordered_set<uint32_t>> &evicted_r_set,
                                        uint32_t r,
                                        uint32_t i,
                                        uint32_t j);

        static void removal_unsatisfied_vertices(const shared_ptr<abstract_bipartite_graph> &B,
                                                 const shared_ptr<unordered_set<shared_ptr<abstract_bipartite_edge>>>& edge_set,
                                                 const shared_ptr<unordered_map<uint32_t, uint32_t>> &l_map,
                                                 const shared_ptr<unordered_map<uint32_t, uint32_t>> &r_map,
                                                 const shared_ptr<unordered_set<uint32_t>> &l_set,
                                                 const shared_ptr<unordered_set<uint32_t>> &r_set,
                                                 uint32_t i,
                                                 uint32_t j);

        static uint32_t update_single_core(const shared_ptr<abstract_bipartite_graph> &B,
                                       const shared_ptr<unordered_set<shared_ptr<abstract_bipartite_edge>>> &edge_set,
                                       const shared_ptr<unordered_map<uint32_t, shared_ptr<bipartite_core_left_store_index>>> &left_index_map,
                                       const shared_ptr<unordered_map<uint32_t, shared_ptr<bipartite_core_right_store_index>>> &right_index_map,
                                       const shared_ptr<unordered_set<uint32_t>> &removed_l_set,
                                       const shared_ptr<unordered_set<uint32_t>> &removed_r_set,
                                       uint32_t i,
                                       uint32_t j);

        static void update_trivial_bipartite_cores(const shared_ptr<abstract_bipartite_graph> &B,
                                                   const shared_ptr<unordered_set<shared_ptr<abstract_bipartite_edge>>> &edge_set,
                                                   const shared_ptr<unordered_map<uint32_t, shared_ptr<bipartite_core_left_store_index>>> &left_index_map,
                                                   const shared_ptr<unordered_map<uint32_t, shared_ptr<bipartite_core_right_store_index>>> &right_index_map);

        static void update_trivial_bipartite_cores(const shared_ptr<abstract_bipartite_graph> &B,
                                                   const shared_ptr<unordered_set<shared_ptr<abstract_bipartite_edge>>> &edge_set,
                                                   const shared_ptr<unordered_map<uint32_t, shared_ptr<bipartite_core_left_store_index>>> &left_index_map,
                                                   const shared_ptr<unordered_map<uint32_t, shared_ptr<bipartite_core_right_store_index>>> &right_index_map,
                                                   const shared_ptr<thread_pool>& pool);
    };
}


