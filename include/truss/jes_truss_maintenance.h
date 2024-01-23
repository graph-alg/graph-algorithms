/******************************************************************************************************************
* @copyright  : Boost Software License
* @file       : jes_truss_maintenance.h
* @brief      : a parallel truss maintenance
* @version    : 1.0
* @date       : 2022/06/18
******************************************************************************************************************/

#pragma once
#include "truss/truss_utility.h"

namespace scnu
{
    class jes_truss_maintenance
    {
    public:

        static void insert(const shared_ptr<abstract_graph> &G,
                           const shared_ptr<unordered_set<shared_ptr<abstract_edge>>> &EI,
                           const shared_ptr<unordered_map<shared_ptr<abstract_edge>,uint32_t>> &edge_truss_map,
                           const shared_ptr<thread_pool>& pool);

        static void remove(const shared_ptr<abstract_graph> &G,
                           const shared_ptr<unordered_set<shared_ptr<abstract_edge>>> &ED,
                           const shared_ptr<unordered_map<shared_ptr<abstract_edge>,uint32_t>> &edge_truss_map,
                           const shared_ptr<thread_pool>& pool);

        static void compute_insert_edge_set(const shared_ptr<abstract_graph> &G,
                                            const shared_ptr<unordered_set<shared_ptr<abstract_edge>>> &EI,
                                            const shared_ptr<unordered_map<shared_ptr<abstract_edge>, uint32_t>> &edge_truss_map,
                                            const shared_ptr<unordered_set<shared_ptr<abstract_edge>>> &Ec,
                                            const shared_ptr<unordered_map<uint32_t, shared_ptr<unordered_set<shared_ptr<abstract_edge>>>>> &EU);

        static uint32_t compute_truss_number(const shared_ptr<abstract_graph> &G,
                                             const shared_ptr<unordered_map<shared_ptr<abstract_edge>, uint32_t>> &edge_truss_map,
                                             const shared_ptr<abstract_edge>& e);

        static void compute_delete_edge_set(const shared_ptr<abstract_graph> &G,
                                            const shared_ptr<unordered_set<shared_ptr<abstract_edge>>> &ED,
                                            const shared_ptr<unordered_map<shared_ptr<abstract_edge>, uint32_t>> &edge_truss_map,
                                            const shared_ptr<unordered_set<shared_ptr<abstract_edge>>> &Ec,
                                            const shared_ptr<unordered_map<uint32_t, shared_ptr<unordered_set<shared_ptr<abstract_edge>>>>> &EU);



        static uint32_t compute_pre_truss_number(const shared_ptr<abstract_graph>& G,
                                                 const shared_ptr<unordered_map<shared_ptr<abstract_edge>,uint32_t>> &edge_truss_map,
                                                 const shared_ptr<abstract_edge>& e);

        static void compute_maximal_three_hop_independent_set(const shared_ptr<abstract_graph> &G,
                                                              const shared_ptr<unordered_set<shared_ptr<abstract_edge>>> &superior_edge_set,
                                                              const shared_ptr<unordered_map<shared_ptr<abstract_edge>,
                                                                      shared_ptr<unordered_set<shared_ptr<abstract_edge>>>>> &superior_edge_triangle_map,
                                                              const shared_ptr<unordered_set<shared_ptr<abstract_edge>>> &Ec);

        static uint32_t
        compute_truss_support(const shared_ptr<abstract_graph> &G,
                              const shared_ptr<unordered_map<shared_ptr<abstract_edge>,uint32_t>> &edge_truss_map,
                              const shared_ptr<abstract_edge>& e,
                              uint32_t k);

        static void
        compute_insert_superior_edge_set(const shared_ptr<abstract_graph>&G,
                                         const shared_ptr<unordered_set<shared_ptr<abstract_edge>>> &edge_set,
                                         const shared_ptr<unordered_map<shared_ptr<abstract_edge>, uint32_t>> &edge_truss_map,
                                         const shared_ptr<unordered_set<shared_ptr<abstract_edge>>>& superior_edge_set,
                                         const shared_ptr<unordered_map<shared_ptr<abstract_edge>,
                                                 shared_ptr<unordered_set<shared_ptr<abstract_edge>>>>>& superior_edge_triangle_map);

        static void compute_delete_superior_edge_set(const shared_ptr<abstract_graph> &G,
                                                     const shared_ptr<unordered_set<shared_ptr<abstract_edge>>> &edge_set,
                                                     const shared_ptr<unordered_map<shared_ptr<abstract_edge>, uint32_t>> &edge_truss_map,
                                                     const shared_ptr<unordered_set<shared_ptr<abstract_edge>>> &superior_edge_set,
                                                     const shared_ptr<unordered_map<shared_ptr<abstract_edge>,
                                                             shared_ptr<unordered_set<shared_ptr<abstract_edge>>>>> &superior_edge_triangle_map);

        static void k_joint_delete(const shared_ptr<abstract_graph> &G,
                                   const shared_ptr<unordered_set<shared_ptr<abstract_edge>>> &Ek,
                                   const shared_ptr<unordered_map<shared_ptr<abstract_edge>, uint32_t>> &edge_truss_map,
                                   const shared_ptr<unordered_set<shared_ptr<abstract_edge>>> &candidate_edge_set,
                                   uint32_t k);

        static void k_joint_insert(const shared_ptr<abstract_graph> &G,
                                   const shared_ptr<unordered_set<shared_ptr<abstract_edge>>> &current_edge_set,
                                   const shared_ptr<unordered_map<shared_ptr<abstract_edge>, uint32_t>> &edge_truss_map,
                                   const shared_ptr<unordered_set<shared_ptr<abstract_edge>>> &candidate_edge_set,
                                   const shared_ptr<unordered_map<shared_ptr<abstract_edge>, uint32_t>> &candidate_edge_support_map,
                                   uint32_t k);

        static void remove_unsatisfied_edges(const shared_ptr<abstract_graph> &G,
                                             const shared_ptr<unordered_map<shared_ptr<abstract_edge>, uint32_t>> &edge_truss_map,
                                             const shared_ptr<abstract_edge>& e,
                                             const shared_ptr<unordered_set<shared_ptr<abstract_edge>>>& candidate_edge_set,
                                             const shared_ptr<unordered_map<shared_ptr<abstract_edge>, uint32_t>> &candidate_edge_support_map,
                                             const shared_ptr<unordered_set<shared_ptr<abstract_edge>>>& evicted_set,
                                             uint32_t k);
    };
}





