/******************************************************************************************************************
* @copyright  : Boost Software License
* @file       : index_jes_truss_maintenance.h
* @brief      : a parallel truss maintenance with index
* @version    : 1.0
* @date       : 2023/04/04
******************************************************************************************************************/

#pragma once
#include "truss/truss_utility.h"

namespace scnu
{
    class jes_order_truss_maintenance {
    public:

        static void init(const shared_ptr<abstract_graph> &G,
                         const std::shared_ptr<unordered_map<shared_ptr<abstract_edge>, uint32_t>> &edge_truss_map,
                         const shared_ptr<unordered_map<uint32_t, shared_ptr<map<uint32_t, shared_ptr<unordered_set<uint32_t>>>>>> &VE_map,
                         const shared_ptr<thread_pool> &pool);

        static void insert(const shared_ptr<abstract_graph> &G,
                           const shared_ptr<unordered_set<shared_ptr<abstract_edge>>> &EI,
                           const shared_ptr<unordered_map<shared_ptr<abstract_edge>, uint32_t>> &edge_truss_map,
                           const shared_ptr<unordered_map<uint32_t, shared_ptr<extend_list<int, shared_ptr<abstract_edge>>>>> &truss_order_map,
                           const shared_ptr<thread_pool> &pool);

        static uint32_t insert(const shared_ptr<abstract_graph> &G,
                               const shared_ptr<unordered_set<shared_ptr<abstract_edge>>> &EI,
                               const shared_ptr<unordered_map<shared_ptr<abstract_edge>, uint32_t>> &edge_truss_map,
                               const shared_ptr<unordered_map<uint32_t, shared_ptr<extend_list<int, shared_ptr<abstract_edge>>>>> &truss_order_map,
                               const shared_ptr<unordered_map<shared_ptr<abstract_edge>, uint32_t>> &TS,
                               const shared_ptr<thread_pool> &pool);

        static uint32_t insert(const shared_ptr<abstract_graph> &G,
                               const shared_ptr<unordered_set<shared_ptr<abstract_edge>>> &EI,
                               const shared_ptr<unordered_map<shared_ptr<abstract_edge>, uint32_t>> &edge_truss_map,
                               const shared_ptr<unordered_map<uint32_t, shared_ptr<extend_list<int, shared_ptr<abstract_edge>>>>> &truss_order_map,
                               const shared_ptr<unordered_map<shared_ptr<abstract_edge>, uint32_t>> &rem,
                               const shared_ptr<unordered_map<shared_ptr<abstract_edge>, uint32_t>> &TS,
                               const shared_ptr<thread_pool> &pool);

        static void insert(const shared_ptr<abstract_graph> &G,
                           const shared_ptr<unordered_set<shared_ptr<abstract_edge>>> &EI,
                           const shared_ptr<unordered_map<shared_ptr<abstract_edge>, uint32_t>> &edge_truss_map,
                           const shared_ptr<unordered_map<uint32_t, shared_ptr<extend_list<int, shared_ptr<abstract_edge>>>>> &truss_order_map,
                           const shared_ptr<unordered_map<shared_ptr<abstract_edge>, uint32_t>> &rem,
                           const shared_ptr<unordered_map<shared_ptr<abstract_edge>, uint32_t>> &TS,
                           const shared_ptr<unordered_map<uint32_t, shared_ptr<map<uint32_t, shared_ptr<unordered_set<uint32_t>>>>>> &VE_map,
                           const shared_ptr<thread_pool> &pool);

        static void remove(const shared_ptr<abstract_graph> &G,
                           const shared_ptr<unordered_set<shared_ptr<abstract_edge>>> &ED,
                           const shared_ptr<unordered_map<shared_ptr<abstract_edge>, uint32_t>> &edge_truss_map,
                           const shared_ptr<unordered_map<uint32_t, shared_ptr<extend_list<int, shared_ptr<abstract_edge>>>>> &truss_order_map,
                           const shared_ptr<thread_pool> &pool);

        static uint32_t remove(const shared_ptr<abstract_graph> &G,
                               const shared_ptr<unordered_set<shared_ptr<abstract_edge>>> &ED,
                               const shared_ptr<unordered_map<shared_ptr<abstract_edge>, uint32_t>> &edge_truss_map,
                               const shared_ptr<unordered_map<uint32_t, shared_ptr<extend_list<int, shared_ptr<abstract_edge>>>>> &truss_order_map,
                               const shared_ptr<unordered_map<shared_ptr<abstract_edge>, uint32_t>> &TS,
                               const shared_ptr<thread_pool> &pool);

        static uint32_t remove(const shared_ptr<abstract_graph> &G,
                               const shared_ptr<unordered_set<shared_ptr<abstract_edge>>> &ED,
                               const shared_ptr<unordered_map<shared_ptr<abstract_edge>, uint32_t>> &edge_truss_map,
                               const shared_ptr<unordered_map<uint32_t, shared_ptr<extend_list<int, shared_ptr<abstract_edge>>>>> &truss_order_map,
                               const shared_ptr<unordered_map<shared_ptr<abstract_edge>, uint32_t>> &rem,
                               const shared_ptr<unordered_map<shared_ptr<abstract_edge>, uint32_t>> &TS,
                               const shared_ptr<thread_pool> &pool);

        static void remove(const shared_ptr<abstract_graph> &G,
                           const shared_ptr<unordered_set<shared_ptr<abstract_edge>>> &ED,
                           const shared_ptr<unordered_map<shared_ptr<abstract_edge>, uint32_t>> &edge_truss_map,
                           const shared_ptr<unordered_map<uint32_t, shared_ptr<extend_list<int, shared_ptr<abstract_edge>>>>> &truss_order_map,
                           const shared_ptr<unordered_map<shared_ptr<abstract_edge>, uint32_t>> &rem,
                           const shared_ptr<unordered_map<shared_ptr<abstract_edge>, uint32_t>> &TS,
                           const shared_ptr<unordered_map<uint32_t, shared_ptr<map<uint32_t, shared_ptr<unordered_set<uint32_t>>>>>> &VE_map,
                           const shared_ptr<thread_pool> &pool);

    private:

        static uint32_t compute_truss_number(const shared_ptr<abstract_graph> &G,
                                             const shared_ptr<unordered_map<shared_ptr<abstract_edge>, uint32_t>> &edge_truss_map,
                                             const shared_ptr<abstract_edge>& e);

        static void compute_delete_edge_set(const shared_ptr<abstract_graph> &G,
                                            const shared_ptr<unordered_set<shared_ptr<abstract_edge>>> &ED,
                                            const shared_ptr<unordered_map<shared_ptr<abstract_edge>, uint32_t>> &edge_truss_map,
                                            const shared_ptr<unordered_map<uint32_t, shared_ptr<extend_list<int, shared_ptr<abstract_edge>>>>>& truss_order_map,
                                            const shared_ptr<unordered_set<shared_ptr<abstract_edge>>> &Ec,
                                            const shared_ptr<unordered_map<uint32_t, shared_ptr<unordered_set<shared_ptr<abstract_edge>>>>> &EU,
                                            const shared_ptr<thread_pool>& pool);

        static void compute_delete_edge_set(const shared_ptr<abstract_graph> &G,
                                            const shared_ptr<unordered_set<shared_ptr<abstract_edge>>> &ED,
                                            const shared_ptr<unordered_map<shared_ptr<abstract_edge>, uint32_t>> &edge_truss_map,
                                            const shared_ptr<unordered_map<uint32_t, shared_ptr<extend_list<int, shared_ptr<abstract_edge>>>>>& truss_order_map,
                                            const shared_ptr<unordered_map<shared_ptr<abstract_edge>, uint32_t>> &TS,
                                            const shared_ptr<unordered_set<shared_ptr<abstract_edge>>> &Ec,
                                            const shared_ptr<unordered_map<uint32_t, shared_ptr<unordered_set<shared_ptr<abstract_edge>>>>> &EU,
                                            const shared_ptr<thread_pool>& pool);

        static void compute_delete_edge_set(const shared_ptr<abstract_graph> &G,
                                            const shared_ptr<unordered_set<shared_ptr<abstract_edge>>> &ED,
                                            const shared_ptr<unordered_map<shared_ptr<abstract_edge>, uint32_t>> &edge_truss_map,
                                            const shared_ptr<unordered_map<uint32_t, shared_ptr<extend_list<int, shared_ptr<abstract_edge>>>>>& truss_order_map,
                                            const shared_ptr<unordered_map<shared_ptr<abstract_edge>, uint32_t>> &rem,
                                            const shared_ptr<unordered_map<shared_ptr<abstract_edge>, uint32_t>> &TS,
                                            const shared_ptr<unordered_set<shared_ptr<abstract_edge>>> &Ec,
                                            const shared_ptr<unordered_map<uint32_t, shared_ptr<unordered_set<shared_ptr<abstract_edge>>>>> &EU,
                                            const shared_ptr<thread_pool>& pool);

        static void compute_delete_edge_set(const shared_ptr<abstract_graph> &G,
                                            const shared_ptr<unordered_set<shared_ptr<abstract_edge>>> &ED,
                                            const shared_ptr<unordered_map<shared_ptr<abstract_edge>, uint32_t>> &edge_truss_map,
                                            const shared_ptr<unordered_map<uint32_t, shared_ptr<extend_list<int, shared_ptr<abstract_edge>>>>> &truss_order_map,
                                            const shared_ptr<unordered_map<shared_ptr<abstract_edge>, uint32_t>> &rem,
                                            const shared_ptr<unordered_map<shared_ptr<abstract_edge>, uint32_t>> &TS,
                                            const shared_ptr<unordered_map<uint32_t, shared_ptr<map<uint32_t, shared_ptr<unordered_set<uint32_t>>>>>> &VE_map,
                                            const shared_ptr<unordered_set<shared_ptr<abstract_edge>>> &Ec,
                                            const shared_ptr<unordered_map<uint32_t, shared_ptr<unordered_set<shared_ptr<abstract_edge>>>>> &EU,
                                            const shared_ptr<thread_pool> &pool);

        static void compute_insert_edge_set(const shared_ptr<abstract_graph> &G,
                                            const shared_ptr<unordered_set<shared_ptr<abstract_edge>>> &EI,
                                            const shared_ptr<unordered_map<shared_ptr<abstract_edge>, uint32_t>> &edge_truss_map,
                                            const shared_ptr<unordered_map<uint32_t, shared_ptr<extend_list<int, shared_ptr<abstract_edge>>>>>& truss_order_map,
                                            const shared_ptr<unordered_set<shared_ptr<abstract_edge>>> &Ec,
                                            const shared_ptr<unordered_map<uint32_t, shared_ptr<unordered_set<shared_ptr<abstract_edge>>>>> &EU,
                                            const shared_ptr<thread_pool>& pool);

        static void compute_insert_edge_set(const shared_ptr<abstract_graph> &G,
                                            const shared_ptr<unordered_set<shared_ptr<abstract_edge>>> &EI,
                                            const shared_ptr<unordered_map<shared_ptr<abstract_edge>, uint32_t>> &edge_truss_map,
                                            const shared_ptr<unordered_map<uint32_t, shared_ptr<extend_list<int, shared_ptr<abstract_edge>>>>>& truss_order_map,
                                            const shared_ptr<unordered_map<shared_ptr<abstract_edge>, uint32_t>> &TS,
                                            const shared_ptr<unordered_set<shared_ptr<abstract_edge>>> &Ec,
                                            const shared_ptr<unordered_map<uint32_t, shared_ptr<unordered_set<shared_ptr<abstract_edge>>>>> &EU,
                                            const shared_ptr<thread_pool>& pool);

        static void compute_insert_edge_set(const shared_ptr<abstract_graph> &G,
                                            const shared_ptr<unordered_set<shared_ptr<abstract_edge>>> &EI,
                                            const shared_ptr<unordered_map<shared_ptr<abstract_edge>, uint32_t>> &edge_truss_map,
                                            const shared_ptr<unordered_map<uint32_t, shared_ptr<extend_list<int, shared_ptr<abstract_edge>>>>> &truss_order_map,
                                            const shared_ptr<unordered_map<shared_ptr<abstract_edge>, uint32_t>> &rem,
                                            const shared_ptr<unordered_map<shared_ptr<abstract_edge>, uint32_t>> &TS,
                                            const shared_ptr<unordered_set<shared_ptr<abstract_edge>>> &Ec,
                                            const shared_ptr<unordered_map<uint32_t, shared_ptr<unordered_set<shared_ptr<abstract_edge>>>>> &EU,
                                            const shared_ptr<thread_pool> &pool);

        static void compute_insert_edge_set(const shared_ptr<abstract_graph> &G,
                                            const shared_ptr<unordered_set<shared_ptr<abstract_edge>>> &EI,
                                            const shared_ptr<unordered_map<shared_ptr<abstract_edge>, uint32_t>> &edge_truss_map,
                                            const shared_ptr<unordered_map<uint32_t, shared_ptr<extend_list<int, shared_ptr<abstract_edge>>>>> &truss_order_map,
                                            const shared_ptr<unordered_map<shared_ptr<abstract_edge>, uint32_t>> &rem,
                                            const shared_ptr<unordered_map<shared_ptr<abstract_edge>, uint32_t>> &TS,
                                            const shared_ptr<unordered_map<uint32_t, shared_ptr<map<uint32_t, shared_ptr<unordered_set<uint32_t>>>>>> &VE_map,
                                            const shared_ptr<unordered_set<shared_ptr<abstract_edge>>> &Ec,
                                            const shared_ptr<unordered_map<uint32_t, shared_ptr<unordered_set<shared_ptr<abstract_edge>>>>> &EU,
                                            const shared_ptr<thread_pool> &pool);

        static uint32_t compute_pre_truss_number(const shared_ptr<abstract_graph> &G,
                                                 const shared_ptr<unordered_map<shared_ptr<abstract_edge>, uint32_t>> &edge_truss_map,
                                                 const shared_ptr<abstract_edge> &e);


        static void compute_insertion_maximal_three_hop_independent_set(const shared_ptr<abstract_graph> &G,
                                                                        const shared_ptr<unordered_map<shared_ptr<abstract_edge>, uint32_t>> &edge_truss_map,
                                                                        const shared_ptr<unordered_set<shared_ptr<abstract_edge>>> &superior_edge_set,
                                                                        const shared_ptr<unordered_set<shared_ptr<abstract_edge>>> &Ec);

        static void compute_removal_maximal_three_hop_independent_set(const shared_ptr<abstract_graph> &G,
                                                                      const shared_ptr<unordered_map<shared_ptr<abstract_edge>, uint32_t>> &edge_truss_map,
                                                                      const shared_ptr<unordered_set<shared_ptr<abstract_edge>>> &superior_edge_set,
                                                                      const shared_ptr<unordered_set<shared_ptr<abstract_edge>>> &Ec);

        static uint32_t get_rem_set(const shared_ptr<abstract_graph> &G,
                                    const shared_ptr<unordered_map<shared_ptr<abstract_edge>, uint32_t>> &edge_truss_map,
                                    const shared_ptr<unordered_map<uint32_t, shared_ptr<extend_list<int, shared_ptr<abstract_edge>>>>> &truss_order_map,
                                    const shared_ptr<abstract_edge> &e,
                                const shared_ptr<unordered_set<shared_ptr<abstract_edge>>> &rem_set,
                                uint32_t k);

        static shared_ptr<unordered_set<shared_ptr<abstract_edge>>> get_rem_set(const shared_ptr<abstract_graph> &G,
                                                                                const shared_ptr<unordered_map<shared_ptr<abstract_edge>, uint32_t>> &edge_truss_map,
                                                                                const shared_ptr<extend_list<int, shared_ptr<abstract_edge>>> &order_list,
                                                                                const shared_ptr<abstract_edge> &e,
                                                                                uint32_t k);

        static shared_ptr<unordered_set<shared_ptr<abstract_edge>>> get_rem_set(const shared_ptr<abstract_graph> &G,
                                                                                const shared_ptr<unordered_map<shared_ptr<abstract_edge>, uint32_t>> &edge_truss_map,
                                                                                const shared_ptr<extend_list<int, shared_ptr<abstract_edge>>> &order_list,
                                                                                const shared_ptr<unordered_map<uint32_t, shared_ptr<map<uint32_t, shared_ptr<unordered_set<uint32_t>>>>>> &VE_map,
                                                                                const shared_ptr<abstract_edge> &e,
                                                                                uint32_t k);

        static uint32_t
        compute_truss_support(const shared_ptr<abstract_graph> &G,
                              const shared_ptr<unordered_map<shared_ptr<abstract_edge>,uint32_t>> &edge_truss_map,
                              const shared_ptr<abstract_edge>& e,
                              uint32_t k);


        static void compute_insert_superior_edge_set(const shared_ptr<abstract_graph> &G,
                                                     const shared_ptr<unordered_set<shared_ptr<abstract_edge>>> &edge_set,
                                                     const shared_ptr<unordered_map<shared_ptr<abstract_edge>, uint32_t>> &edge_truss_map,
                                                     const shared_ptr<unordered_map<uint32_t, shared_ptr<extend_list<int, shared_ptr<abstract_edge>>>>>& truss_order_map,
                                                     const shared_ptr<unordered_set<shared_ptr<abstract_edge>>> &superior_edge_set,
                                                     const shared_ptr<thread_pool> &pool);

        static void compute_insert_superior_edge_set(const shared_ptr<abstract_graph> &G,
                                                     const shared_ptr<unordered_set<shared_ptr<abstract_edge>>> &edge_set,
                                                     const shared_ptr<unordered_map<shared_ptr<abstract_edge>, uint32_t>> &edge_truss_map,
                                                     const shared_ptr<unordered_map<uint32_t, shared_ptr<extend_list<int, shared_ptr<abstract_edge>>>>>& truss_order_map,
                                                     const shared_ptr<unordered_map<shared_ptr<abstract_edge>, uint32_t>> &TS,
                                                     const shared_ptr<unordered_set<shared_ptr<abstract_edge>>> &superior_edge_set,
                                                     const shared_ptr<unordered_map<uint32_t, shared_ptr<unordered_set<shared_ptr<abstract_edge>>>>>& EU,
                                                     const shared_ptr<thread_pool> &pool);

        static void compute_insert_superior_edge_set(const shared_ptr<abstract_graph> &G,
                                                     const shared_ptr<unordered_set<shared_ptr<abstract_edge>>> &edge_set,
                                                     const shared_ptr<unordered_map<uint32_t, shared_ptr<extend_list<int, shared_ptr<abstract_edge>>>>> &truss_order_map,
                                                     const shared_ptr<unordered_set<shared_ptr<abstract_edge>>> &superior_edge_set,
                                                     const shared_ptr<thread_pool> &pool);

        static void compute_insert_superior_edge_set(const shared_ptr<abstract_graph> &G,
                                                     const shared_ptr<unordered_set<shared_ptr<abstract_edge>>> &edge_set,
                                                     const shared_ptr<unordered_map<shared_ptr<abstract_edge>, uint32_t>> &edge_truss_map,
                                                     const shared_ptr<unordered_map<uint32_t, shared_ptr<extend_list<int, shared_ptr<abstract_edge>>>>> &truss_order_map,
                                                     const shared_ptr<unordered_map<shared_ptr<abstract_edge>, uint32_t>> &rem,
                                                     const shared_ptr<unordered_map<shared_ptr<abstract_edge>, uint32_t>> &TS,
                                                     const shared_ptr<unordered_set<shared_ptr<abstract_edge>>> &superior_edge_set,
                                                     const shared_ptr<thread_pool> &pool);

        static void compute_insert_superior_edge_set(const shared_ptr<abstract_graph> &G,
                                                     const shared_ptr<unordered_set<shared_ptr<abstract_edge>>> &edge_set,
                                                     const shared_ptr<unordered_map<shared_ptr<abstract_edge>, uint32_t>> &edge_truss_map,
                                                     const shared_ptr<unordered_map<uint32_t, shared_ptr<extend_list<int, shared_ptr<abstract_edge>>>>> &truss_order_map,
                                                     const shared_ptr<unordered_map<shared_ptr<abstract_edge>, uint32_t>> &rem,
                                                     const shared_ptr<unordered_map<shared_ptr<abstract_edge>, uint32_t>> &TS,
                                                     const shared_ptr<unordered_map<uint32_t, shared_ptr<map<uint32_t, shared_ptr<unordered_set<uint32_t>>>>>> &VE_map,
                                                     const shared_ptr<unordered_set<shared_ptr<abstract_edge>>> &superior_edge_set,
                                                     const shared_ptr<thread_pool> &pool);


        static void compute_delete_superior_edge_set(const shared_ptr<abstract_graph> &G,
                                                     const shared_ptr<unordered_set<shared_ptr<abstract_edge>>> &edge_set,
                                                     const shared_ptr<unordered_map<shared_ptr<abstract_edge>, uint32_t>> &edge_truss_map,
                                                     const shared_ptr<unordered_map<uint32_t, shared_ptr<extend_list<int, shared_ptr<abstract_edge>>>>> &truss_order_map,
                                                     const shared_ptr<unordered_set<shared_ptr<abstract_edge>>> &superior_edge_set,
                                                     const shared_ptr<unordered_map<shared_ptr<abstract_edge>, shared_ptr<unordered_set<shared_ptr<abstract_edge>>>>> &affected_edge_map,
                                                     const shared_ptr<thread_pool> &pool);

        static void compute_delete_superior_edge_set(const shared_ptr<abstract_graph> &G,
                                                     const shared_ptr<unordered_set<shared_ptr<abstract_edge>>> &edge_set,
                                                     const shared_ptr<unordered_map<shared_ptr<abstract_edge>, uint32_t>> &edge_truss_map,
                                                     const shared_ptr<unordered_map<uint32_t, shared_ptr<extend_list<int, shared_ptr<abstract_edge>>>>> &truss_order_map,
                                                     const shared_ptr<unordered_map<shared_ptr<abstract_edge>, uint32_t>> &TS,
                                                     const shared_ptr<unordered_set<shared_ptr<abstract_edge>>> &superior_edge_set,
                                                     const shared_ptr<thread_pool> &pool);

        static void edge_support_computation(const shared_ptr<abstract_graph> &G,
                                             const shared_ptr<unordered_set<shared_ptr<abstract_edge>>> &candidate_edge_set,
                                             const shared_ptr<unordered_map<shared_ptr<abstract_edge>, uint32_t>> &edge_truss_map,
                                             const shared_ptr<unordered_map<shared_ptr<abstract_edge>, uint32_t>> &candidate_edge_support_map,
                                             const shared_ptr<unordered_set<shared_ptr<abstract_edge>>> &evicted_edge_set,
                                             uint32_t k);

        static void k_joint_delete(const shared_ptr<abstract_graph> &G,
                                   const shared_ptr<unordered_set<shared_ptr<abstract_edge>>> &Ek,
                                   const shared_ptr<unordered_map<shared_ptr<abstract_edge>, uint32_t>> &edge_truss_map,
                                   const shared_ptr<vector<shared_ptr<abstract_edge>>> &candidate_edge_vector,
                                   uint32_t k);

        static void k_joint_delete(const shared_ptr<abstract_graph> &G,
                                   const shared_ptr<unordered_set<shared_ptr<abstract_edge>>> &Ek,
                                   const shared_ptr<unordered_map<shared_ptr<abstract_edge>, uint32_t>> &edge_truss_map,
                                   const shared_ptr<unordered_map<shared_ptr<abstract_edge>, uint32_t>> &TS,
                                   const shared_ptr<vector<shared_ptr<abstract_edge>>> &candidate_edge_vector,
                                   uint32_t k);

        static void k_joint_delete(const shared_ptr<abstract_graph> &G,
                                   const shared_ptr<unordered_set<shared_ptr<abstract_edge>>> &Ek,
                                   const shared_ptr<unordered_map<shared_ptr<abstract_edge>, uint32_t>> &edge_truss_map,
                                   const shared_ptr<unordered_map<uint32_t,shared_ptr<extend_list<int, shared_ptr<abstract_edge>>>>> &truss_order_map,
                                   const shared_ptr<unordered_map<shared_ptr<abstract_edge>, uint32_t>> &rem,
                                   const shared_ptr<unordered_map<shared_ptr<abstract_edge>, uint32_t>> &TS,
                                   const shared_ptr<vector<shared_ptr<abstract_edge>>> &candidate_edge_vector,
                                   uint32_t k);

        static void k_joint_delete(const shared_ptr<abstract_graph> &G,
                                   const shared_ptr<unordered_set<shared_ptr<abstract_edge>>> &Ek,
                                   const shared_ptr<unordered_map<shared_ptr<abstract_edge>, uint32_t>> &edge_truss_map,
                                   const shared_ptr<unordered_map<uint32_t, shared_ptr<extend_list<int, shared_ptr<abstract_edge>>>>> &truss_order_map,
                                   const shared_ptr<unordered_map<shared_ptr<abstract_edge>, uint32_t>> &rem,
                                   const shared_ptr<unordered_map<shared_ptr<abstract_edge>, uint32_t>> &TS,
                                   const shared_ptr<unordered_map<uint32_t, shared_ptr<map<uint32_t, shared_ptr<unordered_set<uint32_t>>>>>> &VE_map,
                                   const shared_ptr<vector<shared_ptr<abstract_edge>>> &candidate_edge_vector,
                                   uint32_t k);


        static void k_joint_insert(const shared_ptr<abstract_graph> &G,
                                   const shared_ptr<unordered_set<shared_ptr<abstract_edge>>> &candidate_edge_set,
                                   const shared_ptr<unordered_map<shared_ptr<abstract_edge>, uint32_t>> &edge_truss_map,
                                   const shared_ptr<unordered_map<shared_ptr<abstract_edge>, uint32_t>> &candidate_edge_support_map,
                                   const shared_ptr<unordered_map<uint32_t, shared_ptr<extend_list<int, shared_ptr<abstract_edge>>>>>& truss_order_map,
                                   uint32_t k);

        static void k_joint_insert(const shared_ptr<abstract_graph> &G,
                                   const shared_ptr<unordered_set<shared_ptr<abstract_edge>>> &candidate_edge_set,
                                   const shared_ptr<unordered_map<shared_ptr<abstract_edge>, uint32_t>> &edge_truss_map,
                                   const shared_ptr<unordered_map<shared_ptr<abstract_edge>, uint32_t>> &candidate_edge_support_map,
                                   const shared_ptr<unordered_map<uint32_t, shared_ptr<extend_list<int, shared_ptr<abstract_edge>>>>>& truss_order_map,
                                   const shared_ptr<unordered_map<shared_ptr<abstract_edge>, uint32_t>> &TS,
                                   uint32_t k);

        static void k_joint_insert(const shared_ptr<abstract_graph> &G,
                                   const shared_ptr<unordered_set<shared_ptr<abstract_edge>>> &candidate_edge_set,
                                   const shared_ptr<unordered_map<shared_ptr<abstract_edge>, uint32_t>> &edge_truss_map,
                                   const shared_ptr<unordered_map<shared_ptr<abstract_edge>, uint32_t>> &candidate_edge_support_map,
                                   const shared_ptr<unordered_map<uint32_t, shared_ptr<extend_list<int, shared_ptr<abstract_edge>>>>> &truss_order_map,
                                   const shared_ptr<unordered_map<shared_ptr<abstract_edge>, uint32_t>> &rem,
                                   const shared_ptr<unordered_map<shared_ptr<abstract_edge>, uint32_t>> &TS,
                                   uint32_t k);

        static void k_joint_insert(const shared_ptr<abstract_graph> &G,
                                   const shared_ptr<unordered_set<shared_ptr<abstract_edge>>> &candidate_edge_set,
                                   const shared_ptr<unordered_map<shared_ptr<abstract_edge>, uint32_t>> &edge_truss_map,
                                   const shared_ptr<unordered_map<shared_ptr<abstract_edge>, uint32_t>> &candidate_edge_support_map,
                                   const shared_ptr<unordered_map<uint32_t, shared_ptr<extend_list<int, shared_ptr<abstract_edge>>>>> &truss_order_map,
                                   const shared_ptr<unordered_map<shared_ptr<abstract_edge>, uint32_t>> &rem,
                                   const shared_ptr<unordered_map<shared_ptr<abstract_edge>, uint32_t>> &TS,
                                   const shared_ptr<unordered_map<uint32_t, shared_ptr<map<uint32_t, shared_ptr<unordered_set<uint32_t>>>>>> &VE_map,
                                   uint32_t k);


        static void remove_unsatisfied_edges(const shared_ptr<abstract_graph> &G,
                                             const shared_ptr<unordered_set<shared_ptr<abstract_edge>>> &current_edge_set,
                                             const shared_ptr<unordered_set<shared_ptr<abstract_edge>>> &candidate_edge_set,
                                             const shared_ptr<unordered_map<shared_ptr<abstract_edge>, uint32_t>> &edge_truss_map,
                                             const shared_ptr<unordered_map<uint32_t, shared_ptr<extend_list<int, shared_ptr<abstract_edge>>>>> &truss_order_map,
                                             const shared_ptr<unordered_map<shared_ptr<abstract_edge>, uint32_t>> &candidate_edge_support_map,
                                             const shared_ptr<unordered_map<shared_ptr<abstract_edge>, uint32_t>> &ext,
                                             const shared_ptr<map<int, shared_ptr<extend_node<int, shared_ptr<abstract_edge>>>>> &current_edge_map,
                                             const shared_ptr<abstract_edge> &e_pivot,
                                             uint32_t k);

        static void remove_unsatisfied_edges(const shared_ptr<abstract_graph> &G,
                                             const shared_ptr<unordered_set<shared_ptr<abstract_edge>>> &current_edge_set,
                                             const shared_ptr<unordered_set<shared_ptr<abstract_edge>>> &candidate_edge_set,
                                             const shared_ptr<unordered_map<shared_ptr<abstract_edge>, uint32_t>> &edge_truss_map,
                                             const shared_ptr<unordered_map<uint32_t, shared_ptr<extend_list<int, shared_ptr<abstract_edge>>>>> &truss_order_map,
                                             const shared_ptr<unordered_map<shared_ptr<abstract_edge>, uint32_t>> &rem,
                                             const shared_ptr<unordered_map<shared_ptr<abstract_edge>, uint32_t>> &candidate_edge_support_map,
                                             const shared_ptr<unordered_map<shared_ptr<abstract_edge>, uint32_t>> &ext,
                                             const shared_ptr<map<int, shared_ptr<extend_node<int, shared_ptr<abstract_edge>>>>> &current_edge_map,
                                             const shared_ptr<abstract_edge> &e_pivot,
                                             uint32_t k);

        static void remove_unsatisfied_edges(const shared_ptr<abstract_graph> &G,
                                             const shared_ptr<unordered_set<shared_ptr<abstract_edge>>> &current_edge_set,
                                             const shared_ptr<unordered_set<shared_ptr<abstract_edge>>> &candidate_edge_set,
                                             const shared_ptr<unordered_map<shared_ptr<abstract_edge>, uint32_t>> &edge_truss_map,
                                             const shared_ptr<unordered_map<uint32_t, shared_ptr<extend_list<int, shared_ptr<abstract_edge>>>>> &truss_order_map,
                                             const shared_ptr<unordered_map<shared_ptr<abstract_edge>, uint32_t>> &rem,
                                             const shared_ptr<unordered_map<uint32_t, shared_ptr<map<uint32_t, shared_ptr<unordered_set<uint32_t>>>>>> &VE_map,
                                             const shared_ptr<unordered_map<shared_ptr<abstract_edge>, uint32_t>> &candidate_edge_support_map,
                                             const shared_ptr<unordered_map<shared_ptr<abstract_edge>, uint32_t>> &ext,
                                             const shared_ptr<map<int, shared_ptr<extend_node<int, shared_ptr<abstract_edge>>>>> &current_edge_map,
                                             const shared_ptr<abstract_edge> &e_pivot,
                                             uint32_t k);

        static bool test_order(const shared_ptr<unordered_map<shared_ptr<abstract_edge>, uint32_t>> &edge_truss_map,
                               const shared_ptr<unordered_map<uint32_t, shared_ptr<extend_list<int, shared_ptr<abstract_edge>>>>>&truss_order_map,
                               const shared_ptr<abstract_edge>& e1,
                               const shared_ptr<abstract_edge>& e2);

        static bool test_order(const shared_ptr<extend_list<int, shared_ptr<abstract_edge>>>& order_list,
                               const shared_ptr<abstract_edge>& e1,
                               const shared_ptr<abstract_edge>& e2);

        static void update_edge_truss_support(const shared_ptr<abstract_graph> &G,
                                              const shared_ptr<unordered_set<shared_ptr<abstract_edge>>> &edge_set,
                                              const shared_ptr<unordered_map<shared_ptr<abstract_edge>, uint32_t>> &edge_truss_map,
                                              const shared_ptr<unordered_map<shared_ptr<abstract_edge>, uint32_t>> &TS,
                                              uint32_t k);

        static void update_edge_truss_support(const shared_ptr<abstract_graph> &G,
                                              const shared_ptr<unordered_set<shared_ptr<abstract_edge>>> &edge_set,
                                              const shared_ptr<unordered_map<shared_ptr<abstract_edge>, uint32_t>> &edge_truss_map,
                                              const shared_ptr<unordered_map<shared_ptr<abstract_edge>, uint32_t>> &TS,
                                              const shared_ptr<unordered_map<uint32_t, shared_ptr<map<uint32_t, shared_ptr<unordered_set<uint32_t>>>>>> &VE_map,
                                              uint32_t k);

        static void update_edge_support(const shared_ptr<abstract_graph> &G,
                                        const shared_ptr<unordered_set<shared_ptr<abstract_edge>>> &candidate_edge_set,
                                        const shared_ptr<unordered_map<shared_ptr<abstract_edge>, uint32_t>> &edge_truss_map,
                                        const shared_ptr<unordered_map<uint32_t, shared_ptr<extend_list<int, shared_ptr<abstract_edge>>>>> &truss_order_map,
                                        const shared_ptr<unordered_map<shared_ptr<abstract_edge>, uint32_t>> &ext,
                                        const shared_ptr<map<int, shared_ptr<extend_node<int, shared_ptr<abstract_edge>>>>> &current_edge_map,
                                        uint32_t k);

        static void update_edge_support(const shared_ptr<abstract_graph> &G,
                                        const shared_ptr<unordered_set<shared_ptr<abstract_edge>>> &candidate_edge_set,
                                        const shared_ptr<unordered_map<shared_ptr<abstract_edge>, uint32_t>> &edge_truss_map,
                                        const shared_ptr<unordered_map<uint32_t, shared_ptr<extend_list<int, shared_ptr<abstract_edge>>>>> &truss_order_map,
                                        const shared_ptr<unordered_map<uint32_t, shared_ptr<map<uint32_t, shared_ptr<unordered_set<uint32_t>>>>>> &VE_map,
                                        const shared_ptr<unordered_map<shared_ptr<abstract_edge>, uint32_t>> &ext,
                                        const shared_ptr<map<int, shared_ptr<extend_node<int, shared_ptr<abstract_edge>>>>> &current_edge_map,
                                        uint32_t k);

        static void update_order_list(const shared_ptr<abstract_graph> &G,
                                      const shared_ptr<unordered_set<shared_ptr<abstract_edge>>> &current_edge_set,
                                      const shared_ptr<unordered_set<shared_ptr<abstract_edge>>> &candidate_edge_set,
                                      const shared_ptr<unordered_map<shared_ptr<abstract_edge>, uint32_t>> &candidate_edge_support_map,
                                      const shared_ptr<unordered_map<shared_ptr<abstract_edge>, uint32_t>> &edge_truss_map,
                                      const shared_ptr<extend_list<int, shared_ptr<abstract_edge>>> &current_order_list,
                                      const shared_ptr<vector<shared_ptr<scnu::abstract_edge>>> &e_edge_vector,
                                      uint32_t k);

        static void update_order_list(const shared_ptr<abstract_graph> &G,
                                      const shared_ptr<unordered_set<shared_ptr<abstract_edge>>> &current_edge_set,
                                      const shared_ptr<unordered_set<shared_ptr<abstract_edge>>> &candidate_edge_set,
                                      const shared_ptr<unordered_map<shared_ptr<abstract_edge>, uint32_t>> &candidate_edge_support_map,
                                      const shared_ptr<unordered_map<shared_ptr<abstract_edge>, uint32_t>> &edge_truss_map,
                                      const shared_ptr<extend_list<int, shared_ptr<abstract_edge>>> &current_order_list,
                                      const shared_ptr<unordered_map<shared_ptr<abstract_edge>, uint32_t>> &TS,
                                      const shared_ptr<vector<shared_ptr<scnu::abstract_edge>>> &e_edge_vector,
                                      uint32_t k);

        static void update_order_list(const shared_ptr<abstract_graph> &G,
                                      const shared_ptr<unordered_set<shared_ptr<abstract_edge>>> &current_edge_set,
                                      const shared_ptr<unordered_set<shared_ptr<abstract_edge>>> &candidate_edge_set,
                                      const shared_ptr<unordered_map<shared_ptr<abstract_edge>, uint32_t>> &candidate_edge_support_map,
                                      const shared_ptr<unordered_map<shared_ptr<abstract_edge>, uint32_t>> &edge_truss_map,
                                      const shared_ptr<extend_list<int, shared_ptr<abstract_edge>>> &current_order_list,
                                      const shared_ptr<unordered_map<shared_ptr<abstract_edge>, uint32_t>> &rem,
                                      const shared_ptr<unordered_map<shared_ptr<abstract_edge>, uint32_t>> &TS,
                                      const shared_ptr<vector<shared_ptr<scnu::abstract_edge>>> &e_edge_vector,
                                      uint32_t k);

        static void update_order_list(const shared_ptr<abstract_graph> &G,
                                      const shared_ptr<unordered_set<shared_ptr<abstract_edge>>> &current_edge_set,
                                      const shared_ptr<unordered_set<shared_ptr<abstract_edge>>> &candidate_edge_set,
                                      const shared_ptr<unordered_map<shared_ptr<abstract_edge>, uint32_t>> &candidate_edge_support_map,
                                      const shared_ptr<unordered_map<shared_ptr<abstract_edge>, uint32_t>> &edge_truss_map,
                                      const shared_ptr<extend_list<int, shared_ptr<abstract_edge>>> &current_order_list,
                                      const shared_ptr<unordered_map<shared_ptr<abstract_edge>, uint32_t>> &rem,
                                      const shared_ptr<unordered_map<shared_ptr<abstract_edge>, uint32_t>> &TS,
                                      const shared_ptr<unordered_map<uint32_t, shared_ptr<map<uint32_t, shared_ptr<unordered_set<uint32_t>>>>>> &VE_map,
                                      const shared_ptr<vector<shared_ptr<scnu::abstract_edge>>> &e_edge_vector,
                                      uint32_t k);

        template<class container_type>
        static void update_VE_map_for_insertion(const container_type &edge_container,
                                                const shared_ptr<unordered_map<shared_ptr<abstract_edge>, uint32_t>>& edge_truss_map,
                                                const shared_ptr<unordered_map<uint32_t, shared_ptr<map<uint32_t, shared_ptr<unordered_set<uint32_t>>>>>> &VE_map){
            for (const auto &e: *edge_container) {
                auto u = e->get_source_vertex_id();
                auto v = e->get_destination_vertex_id();

                auto k = edge_truss_map->at(e);

                if(!VE_map->count(u)){
                    VE_map->insert({u, make_shared<map<uint32_t, shared_ptr<unordered_set<uint32_t>>>>()});
                }
                auto u_map = VE_map->at(u);
                if (!u_map->count(k)) {
                    u_map->insert({k, make_shared<unordered_set<uint32_t>>()});
                }
                u_map->at(k)->insert(v);

                if(!VE_map->count(v)){
                    VE_map->insert({v, make_shared<map<uint32_t, shared_ptr<unordered_set<uint32_t>>>>()});
                }
                auto v_map = VE_map->at(v);
                if (!v_map->count(k)) {
                    v_map->insert({k, make_shared<unordered_set<uint32_t>>()});
                }
                v_map->at(k)->insert(u);
            }
        }

        template<class container_type>
        static void update_VE_map_for_removal(const container_type &edge_container,
                                              const shared_ptr<unordered_map<shared_ptr<abstract_edge>, uint32_t>>& edge_truss_map,
                                              const shared_ptr<unordered_map<uint32_t, shared_ptr<map<uint32_t, shared_ptr<unordered_set<uint32_t>>>>>> &VE_map){
            for (const auto &e: *edge_container) {
                auto u = e->get_source_vertex_id();
                auto v = e->get_destination_vertex_id();

                auto k = edge_truss_map->at(e);

                auto u_map = VE_map->at(u);
                if (u_map->count(k)) {
                    u_map->at(k)->erase(v);
                    if(u_map->at(k)->empty()){
                        u_map->erase(k);
                        if(u_map->empty()){
                            VE_map->erase(u);
                        }
                    }
                }

                auto v_map = VE_map->at(v);
                if (v_map->count(k)) {
                    v_map->at(k)->erase(u);
                    if(v_map->at(k)->empty()){
                        v_map->erase(k);
                        if(v_map->empty()){
                            VE_map->erase(v);
                        }
                    }
                }
            }
        }
    };
}
