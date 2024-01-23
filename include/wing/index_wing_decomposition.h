/******************************************************************************************************************
* @copyright  : Boost Software License
* @file       : index_wing_decomposition.h
* @version    : 1.0
* @date       : 2021/4/6
******************************************************************************************************************/

#pragma once
#include "wing_utility.h"
#include "BE_index.h"

namespace scnu{

    /**
     * @brief a index decomposition method
     * @cite Wang K,  Lin X,  Qin L, et al. Towards efficient solutions of bitruss decomposition for
     * large-scale bipartite graphs[J]. The VLDB Journal, 2021(1):1-24.
     */
    class index_wing_decomposition{
    public:

        static void init(const shared_ptr<abstract_bipartite_graph> &G,
                         const shared_ptr<unordered_map<shared_ptr<abstract_bipartite_edge>, shared_ptr<mutex>>> &edge_mutex_map,
                         const shared_ptr<unordered_map<shared_ptr<abstract_bipartite_edge>, uint32_t>> &edge_rank_map,
                         const shared_ptr<unordered_map<shared_ptr<abstract_bipartite_edge>, uint32_t>> &edge_support_map,
                         const shared_ptr<unordered_map<shared_ptr<abstract_bipartite_edge>, uint32_t>> &edge_wing_map,
                         const shared_ptr<thread_pool> &pool);

        static uint32_t WS_init(const shared_ptr<abstract_bipartite_graph> &G,
                                const shared_ptr<unordered_map<shared_ptr<abstract_bipartite_edge>, shared_ptr<mutex>>> &edge_mutex_map,
                                const shared_ptr<unordered_map<shared_ptr<abstract_bipartite_edge>, uint32_t>> &edge_rank_map,
                                const shared_ptr<unordered_map<shared_ptr<abstract_bipartite_edge>, uint32_t>> &edge_wing_map,
                                const shared_ptr<unordered_map<shared_ptr<abstract_bipartite_edge>, uint32_t>> &WS,
                                const shared_ptr<thread_pool> &pool);

        static uint32_t WS_WL_init(const shared_ptr<abstract_bipartite_graph> &G,
                                   const shared_ptr<unordered_map<shared_ptr<abstract_bipartite_edge>,shared_ptr<mutex>>>& edge_mutex_map,
                                   const shared_ptr<unordered_map<shared_ptr<abstract_bipartite_edge>, uint32_t>>& edge_rank_map,
                                   const shared_ptr<unordered_map<shared_ptr<abstract_bipartite_edge>, uint32_t>> &edge_wing_map,
                                   const shared_ptr<unordered_map<shared_ptr<abstract_bipartite_edge>, uint32_t>>& WS,
                                   const shared_ptr<unordered_map<uint32_t,shared_ptr<map<uint32_t,shared_ptr<unordered_set<uint32_t>>>>>> &WL,
                                   const shared_ptr<thread_pool>& pool);

        static void decompose(const shared_ptr<abstract_bipartite_graph>& G,
                              const shared_ptr<unordered_map<shared_ptr<abstract_bipartite_edge>, shared_ptr<mutex>>> &edge_mutex_map,
                              const shared_ptr<unordered_map<shared_ptr<abstract_bipartite_edge>, uint32_t>> &edge_rank_map,
                              const shared_ptr<unordered_map<shared_ptr<abstract_bipartite_edge>, uint32_t>> &edge_support_map,
                              const shared_ptr<unordered_map<shared_ptr<abstract_bipartite_edge>, uint32_t>> &edge_wing_map,
                              const shared_ptr<thread_pool> &pool);

        static void decompose(const shared_ptr<abstract_bipartite_graph>& G,
                              const shared_ptr<unordered_map<shared_ptr<abstract_bipartite_edge>, shared_ptr<mutex>>> &edge_mutex_map,
                              const shared_ptr<unordered_map<shared_ptr<abstract_bipartite_edge>, uint32_t>> &edge_rank_map,
                              const shared_ptr<unordered_map<shared_ptr<abstract_bipartite_edge>, uint32_t>> &edge_support_map,
                              const shared_ptr<unordered_map<shared_ptr<abstract_bipartite_edge>, uint32_t>> &edge_wing_map,
                              double tau,
                              const shared_ptr<thread_pool> &pool);


    private:
        static shared_ptr<BE_index> compressed_index_construction(const shared_ptr<abstract_bipartite_graph> &G,
                                                                  const shared_ptr<unordered_set<shared_ptr<abstract_bipartite_edge>>> &edge_set,
                                                                  const shared_ptr<unordered_map<shared_ptr<abstract_bipartite_edge>,uint32_t>>& edge_support_map,
                                                                  const shared_ptr<unordered_map<shared_ptr<abstract_bipartite_edge>, uint32_t>> &edge_wing_map,
                                                                  const shared_ptr<unordered_map<shared_ptr<abstract_bipartite_edge>, shared_ptr<mutex>>> &edge_mutex_map,
                                                                  const shared_ptr<thread_pool>& pool);

        static void edge_support_computation(const shared_ptr<abstract_bipartite_graph> &G,
                                             const shared_ptr<unordered_set<shared_ptr<abstract_bipartite_edge>>> &edge_set,
                                             const shared_ptr<unordered_map<shared_ptr<abstract_bipartite_edge>, shared_ptr<mutex>>> &edge_mutex_map,
                                             const shared_ptr<unordered_map<shared_ptr<abstract_bipartite_edge>, uint32_t>> &edge_support_map,
                                             const shared_ptr<unordered_map<uint32_t, uint32_t>> &vertex_priority_map,
                                             const shared_ptr<thread_pool> &pool);

        static uint32_t estimate_maximal_wing(const shared_ptr<unordered_map<shared_ptr<abstract_bipartite_edge>,uint32_t>>& edge_wing_map);


        static shared_ptr<BE_index> index_construction(const shared_ptr<abstract_bipartite_graph> &G,
                                                       const shared_ptr<unordered_set<shared_ptr<abstract_bipartite_edge>>>& edge_set,
                                                       const shared_ptr<unordered_map<shared_ptr<abstract_bipartite_edge>, uint32_t>> &edge_support_map,
                                                       const shared_ptr<unordered_map<shared_ptr<abstract_bipartite_edge>,uint32_t>>& edge_wing_map,
                                                       const shared_ptr<unordered_map<shared_ptr<abstract_bipartite_edge>,shared_ptr<mutex>>>& edge_mutex_map,
                                                       const shared_ptr<thread_pool>& pool);

        static void left_index_construction(const shared_ptr<abstract_bipartite_graph> &G,
                                            const shared_ptr<unordered_map<uint32_t, uint32_t>> &vertex_priority,
                                            const shared_ptr<unordered_map<shared_ptr<abstract_bipartite_edge>, uint32_t>> &edge_support_map,
                                            const shared_ptr<unordered_map<shared_ptr<abstract_bipartite_edge>,shared_ptr<mutex>>>& edge_mutex_map,
                                            const shared_ptr<BE_index> &bloom_index,
                                            const shared_ptr<thread_pool>& pool);

        static void left_index_construction(const shared_ptr<abstract_bipartite_graph> &G,
                                            const shared_ptr<unordered_map<uint32_t,uint32_t>>& vertex_priority,
                                            const shared_ptr<unordered_map<shared_ptr<abstract_bipartite_edge>, uint32_t>> &edge_support_map,
                                            const shared_ptr<unordered_map<shared_ptr<abstract_bipartite_edge>, uint32_t>> &edge_wing_map,
                                            const shared_ptr<unordered_map<shared_ptr<abstract_bipartite_edge>,shared_ptr<mutex>>>& edge_mutex_map,
                                            const shared_ptr<BE_index> &bloom_index,
                                            const shared_ptr<thread_pool>& pool);


        static uint32_t minimal_butterfly_support_finding(const shared_ptr<unordered_map<shared_ptr<abstract_bipartite_edge>, uint32_t>> &edge_support_map,
                                                          const shared_ptr<unordered_map<shared_ptr<abstract_bipartite_edge>, uint32_t>> &edge_wing_map);


        static uint32_t minimal_butterfly_support_finding(const shared_ptr<unordered_map<shared_ptr<abstract_bipartite_edge>, uint32_t>> &edge_support_map,
                                                          const shared_ptr<unordered_map<shared_ptr<abstract_bipartite_edge>, uint32_t>> &edge_wing_map,
                                                          const shared_ptr<thread_pool>& pool);

        static void remove_unsatisfied_edge(const shared_ptr<abstract_bipartite_graph> &G,
                                            const shared_ptr<unordered_set<shared_ptr<abstract_bipartite_edge>>>& edge_set,
                                            const shared_ptr<unordered_map<shared_ptr<abstract_bipartite_edge>, uint32_t>> &edge_support_map,
                                            uint32_t epsilon);


        static void right_index_construction(const shared_ptr<abstract_bipartite_graph> &G,
                                             const shared_ptr<unordered_map<uint32_t, uint32_t>> &vertex_priority_map,
                                             const shared_ptr<unordered_map<shared_ptr<abstract_bipartite_edge>, uint32_t>> &edge_support_map,
                                             const shared_ptr<unordered_map<shared_ptr<abstract_bipartite_edge>, shared_ptr<mutex>>> &edge_mutex_map,
                                             const shared_ptr<BE_index> &bloom_index,
                                             const shared_ptr<thread_pool> &pool);

        static void right_index_construction(const shared_ptr<abstract_bipartite_graph> &G,
                                             const shared_ptr<unordered_map<uint32_t,uint32_t>>& vertex_priority,
                                             const shared_ptr<unordered_map<shared_ptr<abstract_bipartite_edge>, uint32_t>> &edge_support_map,
                                             const shared_ptr<unordered_map<shared_ptr<abstract_bipartite_edge>, uint32_t>> &edge_wing_map,
                                             const shared_ptr<unordered_map<shared_ptr<abstract_bipartite_edge>,shared_ptr<mutex>>>& edge_mutex_map,
                                             const shared_ptr<BE_index> &bloom_index,
                                             const shared_ptr<thread_pool>& pool);

        static void scan(const shared_ptr<unordered_set<shared_ptr<abstract_bipartite_edge>>> &edge_set,
                         const shared_ptr<unordered_map<shared_ptr<abstract_bipartite_edge>, uint32_t>> &edge_support_map,
                         const shared_ptr<unordered_map<shared_ptr<abstract_bipartite_edge>, uint32_t>> &edge_wing_map,
                         uint32_t MBS,
                         const shared_ptr<unordered_set<shared_ptr<abstract_bipartite_edge>>> &S_cur);


        static void scan(const shared_ptr<unordered_set<shared_ptr<abstract_bipartite_edge>>> &edge_set,
                         const shared_ptr<unordered_map<shared_ptr<abstract_bipartite_edge>, uint32_t>> &edge_support_map,
                         const shared_ptr<unordered_map<shared_ptr<abstract_bipartite_edge>, uint32_t>> &edge_wing_map,
                         uint32_t MBS,
                         const shared_ptr<unordered_set<shared_ptr<abstract_bipartite_edge>>> &S_cur,
                         const shared_ptr<thread_pool> &pool);

        static void set_edge_rank(const shared_ptr<unordered_set<shared_ptr<scnu::abstract_bipartite_edge>>> &edge_set,
                                  const shared_ptr<unordered_map<shared_ptr<scnu::abstract_bipartite_edge>, uint32_t>> &edge_rank_map,
                                  const shared_ptr<uint32_t>& rank_id);

        static void vertex_priority_computation(const shared_ptr<abstract_bipartite_graph>& G,
                                                const shared_ptr<unordered_map<uint32_t,uint32_t>>& vertex_priority_map,
                                                const shared_ptr<thread_pool> &pool);

    };
}