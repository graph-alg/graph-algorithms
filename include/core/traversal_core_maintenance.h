/******************************************************************************************************************
* @copyright  : Boost Software License
* @file       : traversal_core_maintenance.h
* @brief      : traversal-based core maintenance methods
* @version    : 1.0
* @date       : 2020/9/7
******************************************************************************************************************/
#pragma once
#include "core/core_utility.h"
#include "core/basic_core_decomposition.h"

namespace scnu
{
    /**
     * @class traversal_core_maintenance
     * @brief traversal-based core maintenance methods
     * @cite incremental k-core decomposition: algorithms and evaluation
     */
    class traversal_core_maintenance
    {
    public:
        static void insert(const std::shared_ptr<abstract_graph> &G,
                           const std::shared_ptr<unordered_map<uint32_t, uint32_t>> &core,
                           const std::shared_ptr<unordered_map<uint32_t, shared_ptr<unordered_map<uint32_t, uint32_t>>>> &rcd,
                           uint32_t n,
                           const std::shared_ptr<abstract_edge> &e);


        static void remove(const std::shared_ptr<abstract_graph> &G,
                           const std::shared_ptr<unordered_map<uint32_t, uint32_t>> &core,
                           const std::shared_ptr<unordered_map<uint32_t, shared_ptr<unordered_map<uint32_t, uint32_t>>>> &rcd,
                           uint32_t n,
                           const std::shared_ptr<abstract_edge> &e);

        static void init_rcd(const std::shared_ptr<abstract_graph> &G,
                             const std::shared_ptr<unordered_map<uint32_t,uint32_t>> &core,
                             const std::shared_ptr<unordered_map<uint32_t,shared_ptr<unordered_map<uint32_t,uint32_t>>>> &rcd,
                             uint32_t n);
    private:
        static uint32_t compute_core_rcd(const std::shared_ptr<abstract_graph> &G, uint32_t k, uint32_t u,
                                    const std::shared_ptr<unordered_map<uint32_t,uint32_t>> &core,
                                    const std::shared_ptr<unordered_map<uint32_t,shared_ptr<unordered_map<uint32_t,uint32_t>>>> &rcd, uint32_t h);



        static int compute_rcd(const std::shared_ptr<abstract_graph> &G,
                               uint32_t u,
                               const std::shared_ptr<unordered_map<uint32_t,uint32_t>> &core,
                               const std::shared_ptr<unordered_map<uint32_t,shared_ptr<unordered_map<uint32_t,uint32_t>>>> &rcd,
                               uint32_t h);


        static void multi_hop_prepare_rcd_insertion(const std::shared_ptr<abstract_graph> &G,
                                                    const std::shared_ptr<unordered_map<uint32_t,uint32_t>> &core,
                                                    const std::shared_ptr<unordered_map<uint32_t,shared_ptr<unordered_map<uint32_t,uint32_t>>>> &rcd,
                                                    uint32_t n,
                                                    const std::shared_ptr<abstract_edge> &e);

        static void multi_hop_recompute_rcd(const std::shared_ptr<abstract_graph> &G,
                                            const std::shared_ptr<unordered_map<uint32_t,uint32_t>> &core,
                                            const std::shared_ptr<unordered_map<uint32_t,shared_ptr<unordered_map<uint32_t,uint32_t>>>> &rcd,
                                            uint32_t n,
                                            const std::shared_ptr<list<uint32_t>> &changed,
                                            bool state);


        static void multi_hop_prepare_rcd_removal(const std::shared_ptr<abstract_graph> &G,
                                                  const std::shared_ptr<unordered_map<uint32_t,uint32_t>> &core,
                                                  const std::shared_ptr<unordered_map<uint32_t,shared_ptr<unordered_map<uint32_t,uint32_t>>>> &rcd,
                                                  uint32_t n,
                                                  const std::shared_ptr<abstract_edge> &e);

        static void propagate_eviction(const std::shared_ptr<abstract_graph> &G,
                                       const std::shared_ptr<unordered_map<uint32_t,uint32_t>> &core,
                                       const std::shared_ptr<unordered_map<uint32_t,long long>> &cd,
                                       const std::shared_ptr<unordered_set<uint32_t>> &evicted,
                                       uint32_t k,
                                       uint32_t v);

        static void propagate_dismissal(const std::shared_ptr<abstract_graph> &G,
                                        const std::shared_ptr<unordered_map<uint32_t,uint32_t>> &core,
                                        const std::shared_ptr<unordered_map<uint32_t,shared_ptr<unordered_map<uint32_t,uint32_t>>>> &rcd,
                                        const std::shared_ptr<unordered_map<uint32_t,long long>> &cd,
                                        const std::shared_ptr<unordered_set<uint32_t>> &dismissed,
                                        const std::shared_ptr<unordered_set<uint32_t>> &visited,
                                        uint32_t k,
                                        uint32_t v,
                                        uint32_t n,
                                        const std::shared_ptr<list<uint32_t>>& changed);
    };
}




