/******************************************************************************************************************
* @copyright  : Boost Software License
* @file       : jes_truss_maintenance.h
* @brief      : a parallel truss maintenance
* @version    : 1.0
* @date       : 2023/03/27
******************************************************************************************************************/

#pragma once
#include "truss/truss_utility.h"

namespace scnu {

    class mixed_structure_truss_maintenance {

    public:

        static void insert(const shared_ptr<abstract_graph> &G,
                           const shared_ptr<unordered_set<shared_ptr<abstract_edge>>> &EI,
                           const shared_ptr<unordered_map<shared_ptr<abstract_edge>, uint32_t>> &edge_truss_map);

        static uint32_t insert(const shared_ptr<abstract_graph> &G,
                               const shared_ptr<unordered_set<shared_ptr<abstract_edge>>> &EI,
                               const shared_ptr<unordered_map<shared_ptr<abstract_edge>, uint32_t>> &edge_truss_map,
                               const shared_ptr<thread_pool> &pool);

        static void remove(const shared_ptr<abstract_graph> &G,
                           const shared_ptr<unordered_set<shared_ptr<abstract_edge>>> &ED,
                           const shared_ptr<unordered_map<shared_ptr<abstract_edge>, uint32_t>> &edge_truss_map);

        static uint32_t remove(const shared_ptr<abstract_graph> &G,
                               const shared_ptr<unordered_set<shared_ptr<abstract_edge>>> &ED,
                               const shared_ptr<unordered_map<shared_ptr<abstract_edge>, uint32_t>> &edge_truss_map,
                               const shared_ptr<thread_pool> &pool);

    private:

        static uint32_t compute_pre_truss_number(const shared_ptr<abstract_graph> &G,
                                                 const shared_ptr<unordered_map<shared_ptr<abstract_edge>, uint32_t>> &edge_truss_map,
                                                 const shared_ptr<abstract_edge> &e);

        static uint32_t compute_truss_support(const shared_ptr<abstract_graph> &G,
                                              const shared_ptr<unordered_map<shared_ptr<abstract_edge>, uint32_t>> &edge_truss_map,
                                              const shared_ptr<::abstract_edge> &e,
                                              uint32_t k);

        static void decremental_traversal(const shared_ptr<abstract_graph> &G,
                                          const shared_ptr<unordered_map<shared_ptr<abstract_edge>, uint32_t>> &edge_truss_map,
                                          const shared_ptr<map<uint32_t, shared_ptr<unordered_set<shared_ptr<abstract_edge>>>>> &M);

        static void decremental_traversal(const shared_ptr<abstract_graph> &G,
                                          const shared_ptr<unordered_map<shared_ptr<abstract_edge>, uint32_t>> &edge_truss_map,
                                          const shared_ptr<map<uint32_t, shared_ptr<unordered_set<shared_ptr<abstract_edge>>>>> &M,
                                          const shared_ptr<thread_pool> &pool);


        static void incremental_traversal(const shared_ptr<abstract_graph> &G,
                                          const shared_ptr<unordered_map<shared_ptr<abstract_edge>, uint32_t>> &edge_truss_map,
                                          const shared_ptr<map<uint32_t, shared_ptr<unordered_set<shared_ptr<abstract_edge>>>>> &M);

        static void incremental_traversal(const shared_ptr<abstract_graph> &G,
                                          const shared_ptr<unordered_map<shared_ptr<abstract_edge>, uint32_t>> &edge_truss_map,
                                          const shared_ptr<map<uint32_t, shared_ptr<unordered_set<shared_ptr<abstract_edge>>>>> &M,
                                          const shared_ptr<thread_pool>& pool);
    };

}
