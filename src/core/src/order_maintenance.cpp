
#include "core/order_maintenance.h"

namespace scnu {

    /**
     * @details batch insert maintenance, it updates the core number of vertices and maintain core-order map
     * remarks: some codes are optimized
     * @param G the given graph
     * @param EI a set of inserted edges
     * @param core a map of vertices and their core numbers
     * @param k_order a map of core numbers and corresponding core-order list 
     * @param rem a map of vertices and remaining degrees 
     * @param ext a map of vertices and candidate degrees
     */
    void order_maintenance::batch_insert(const shared_ptr<abstract_graph> &G,
                                         const shared_ptr<unordered_set<shared_ptr<abstract_edge>>> &EI,
                                         const shared_ptr<unordered_map<uint32_t, uint32_t>> &core,
                                         const shared_ptr<unordered_map<uint32_t, shared_ptr<extend_list<double, uint32_t>>>> &k_order,
                                         const shared_ptr<unordered_map<uint32_t, uint32_t>> &rem,
                                         const shared_ptr<unordered_map<uint32_t, uint32_t>> &ext) {
        G->insert_edge_collection(EI);
        unordered_map<uint32_t, shared_ptr<unordered_set<uint32_t>>> N;
        unordered_set<uint32_t> Ck;
        /**
         * @brief sort affected vertices
         */
        map<uint32_t, shared_ptr<extend_node<double, uint32_t>>> skip_map;
        uint32_t k = 1;
        N.insert({k, make_shared<unordered_set<uint32_t>>()});
        unordered_set<shared_ptr<abstract_edge>> removed_edge_set;

        while (!EI->empty()||!Ck.empty()) {
            for (const auto &e : *EI) {
                auto u = e->get_source_vertex_id();
                auto v = e->get_destination_vertex_id();
                if (!core->count(u)) {
                    core->insert({u, 1});
                    rem->insert({u,0});
                    ext->insert({u, 0});
                    N.at(k)->insert(u);
                }
                if (!core->count(v)) {
                    core->insert({v, 1});
                    rem->insert({v,0});
                    ext->insert({v, 0});
                    N.at(k)->insert(v);
                }
                if(core->at(u)<k||core->at(v)<k)
                {
                    removed_edge_set.insert(e);
                    continue;
                }
                if(k_order->at(k)->count_value(u)&&core->at(u)<=core->at(v))
                {
                    ++ext->at(u);
                    auto u_node = k_order->at(k)->find(u);
                    skip_map.insert({u_node->get_key(), u_node});
                }
                if(k_order->at(k)->count_value(v)&&core->at(u)>core->at(v))
                {
                    ++ext->at(v);
                    auto v_node = k_order->at(k)->find(v);
                    skip_map.insert({v_node->get_key(), v_node});
                }
            }
            for(const auto &e:removed_edge_set)
            {
                EI->erase(e);
            }


            if (!k_order->count(k)) {
                k_order->insert({k, make_shared<extend_list<double, uint32_t>>()});
            }
            auto v0_node = k_order->at(k)->get_head();
            unordered_map<uint32_t, uint32_t> s;
            copy(N.at(k)->begin(), N.at(k)->end(), std::inserter(Ck, Ck.end()));

            /**
             * @remarks Q must be ordered
             */
            queue<uint32_t> Q;
            for (const auto &u:Ck) {
                s.insert({u, 0});
                auto u_vertex = G->get_vertex(u);
                for (const auto&[v, e]:*u_vertex->get_edge_map()) {
                    /**
                     * @brief optimization: v must be located in Ck or k-core
                     */
                    if (Ck.count(v) || core->at(v) >= k) {
                        ++s.at(u);
                    }
                }
                if (s.at(u) <= k) {
                    Q.push(u);
                }
            }

            while (!Q.empty()) {
                auto u = Q.front();
                Q.pop();

                rem->insert({u, s.at(u)});

                auto u_vertex = G->get_vertex(u);
                for (const auto &[v, e]:*u_vertex->get_edge_map()) {
                    if (Ck.count(v)) {
                        --s.at(v);
                        if (s.at(v) <= k) {
                            Q.push(v);
                        }
                    }
                }
                Ck.erase(u);

                auto u_node =  make_shared<extend_node<double, uint32_t>>(0, u);
                k_order->at(k)->insert_before(u_node, v0_node);
                core->at(u) = k;
            }

            for (const auto &u:Ck) {
                auto u_vertex = G->get_vertex(u);
                for (const auto&[v, e] :*u_vertex->get_edge_map()) {
                    if (k_order->at(k)->count_value(v)) {
                        ++ext->at(v);
                        auto v_node = k_order->at(k)->find(v);
                        skip_map.insert({v_node->get_key(), v_node});
                    }
                }
            }


            /**
             * @remarks v0 will be reset
             */
            while (v0_node) {
                auto v0 = v0_node->get_value();
                auto v0_next_node = v0_node->get_next();

                skip_map.erase(v0_node->get_key());
                if (ext->at(v0) == 0) {
                    if (!skip_map.empty()) {
                        v0_next_node = skip_map.begin()->second;
                    } else {
                        v0_next_node.reset();
                    }
                } else if (ext->at(v0) + rem->at(v0) > k) {
                    Ck.insert(v0);
                    s.insert({v0, ext->at(v0) + rem->at(v0)});

                    ext->at(v0) = 0;
                    auto v0_vertex = G->get_vertex(v0);
                    for (const auto &[w, e]:*v0_vertex->get_edge_map()) {
                        if (k_order->at(k)->count_value(w) && test_order(core, k_order, v0, w)) {
                            ++ext->at(w);
                            auto w_node = k_order->at(k)->find(w);
                            skip_map.insert({w_node->get_key(), w_node});
                        }
                    }
                    k_order->at(k)->remove(v0);
                } else {
                    rem->at(v0) = ext->at(v0) + rem->at(v0);
                    ext->at(v0) = 0;
                    auto v0_vertex = G->get_vertex(v0);
                    for (const auto &[w, e]:*v0_vertex->get_edge_map()) {
                        if (Ck.count(w)) {
                            --s.at(w);
                            if (s.at(w) <= k) {
                                Q.push(w);
                            }
                        }
                    }


                    while (!Q.empty()) {
                        auto u = Q.front();
                        Q.pop();

                        Ck.erase(u);
                        auto u_node = make_shared<extend_node<double, uint32_t>>(0, u);
                        k_order->at(k)->insert_before(u_node, v0_next_node);
                        core->at(u) = k;

                        rem->at(u) = s.at(u);

                        auto u_vertex = G->get_vertex(u);
                        for (const auto &[v, e]:*u_vertex->get_edge_map()) {
                            if (Ck.count(v)) {
                                --s.at(v);
                                if (s.at(v) <= k) {
                                    Q.push(v);
                                }
                            }
                            else if (core->at(v)==k&&test_order(core,k_order,u,v)) {
                                --ext->at(v);
                                if(ext->at(v)==0)
                                {
                                    skip_map.erase(k_order->at(k)->find_key(v).value());
                                }
                            }
                        }
                    }
                }
                v0_node = v0_next_node;
            }
            if (!N.count(k + 1)) {
                N.insert({k + 1, make_shared<unordered_set<uint32_t>>()});
            }
            copy(Ck.begin(), Ck.end(), inserter(*N.at(k + 1), N.at(k + 1)->end()));
            for (const auto &v:*N.at(k)) {
                if (!N.at(k + 1)->count(v)) {
                    core->at(v) = k;
                }
            }
            k = k + 1;
        }
    }

    void order_maintenance::batch_insert(const shared_ptr<abstract_graph> &G,
                                         const shared_ptr<unordered_set<shared_ptr<abstract_edge>>> &EI,
                                         const shared_ptr<unordered_map<uint32_t, uint32_t>> &core,
                                         const shared_ptr<unordered_map<uint32_t, shared_ptr<extend_list<double, uint32_t>>>> &k_order,
                                         const shared_ptr<unordered_map<uint32_t ,shared_ptr<extend_node<double,uint32_t>>>>& node_map,
                                         const shared_ptr<unordered_map<uint32_t,shared_ptr<unordered_map<uint32_t,uint32_t>>>> &CD,
                                         uint32_t previous_max_k) {
        /**
         * @brief update 1-core
         */
         if(!CD->count(1))
         {
             CD->insert({1,make_shared<unordered_map<uint32_t ,uint32_t>>()});
             k_order->insert({1,make_shared<extend_list<double,uint32_t>>()});
         }
        auto Ck = make_shared<unordered_set<uint32_t>>();

        for (const auto &e:*EI) {
            G->insert_edge(e);
            auto u = e->get_source_vertex_id();
            auto v = e->get_destination_vertex_id();

            if (!core->count(u)) {
                core->insert({u, 1});
                Ck->insert(u);
                CD->at(1)->insert({u,0});
                auto u_node = make_shared<extend_node<double,uint32_t>>(0,u);
                node_map->insert({u,u_node});
            }
            ++CD->at(1)->at(u);

            if (!core->count(v)) {
                core->insert({v, 1});
                Ck->insert(v);
                CD->at(1)->insert({v,0});
                auto v_node = make_shared<extend_node<double,uint32_t>>(0,v);
                node_map->insert({v,v_node});
            }
            ++CD->at(1)->at(v);
        }
        /**
         * @brief update 1-core
         */
        auto evicted_set = make_shared<unordered_set<uint32_t>>();
        partial_core(G, Ck, core, k_order, node_map, CD,evicted_set, 1);
        int k = 2;
        while (!EI->empty()) {
            if (k > previous_max_k) {
                search_influenced_vertices(G,EI , core,CD, k_order,node_map, Ck, evicted_set,k);
                evicted_set->clear();
                partial_core_decomposition(G, Ck, evicted_set, core,CD, k_order, node_map, k);
                break;
            } else {
                search_influenced_vertices(G, EI, core,CD, k_order,node_map,  Ck, evicted_set,k);
                //evictedSet->clear();
                partial_core(G, Ck, core, k_order, node_map, CD,evicted_set, k);
            }
            ++k;
        }
    }

    /**
     * @details core-order based remove maintenance, it updates the core number of vertices and maintain core-order map
     * @param G the given graph
     * @param ER a set of remove edges
     * @param core a map of vertices and core numbers
     * @param k_order a map of core numbers and corresponding core-order lists
     * @param rem a map of vertices and remaining degrees
     * @param ts a map of vertices and maximal core degrees
     */
    void order_maintenance::batch_remove(const shared_ptr<abstract_graph> &G,
                                         const shared_ptr<unordered_set<shared_ptr<abstract_edge>>> &ER,
                                         const shared_ptr<unordered_map<uint32_t, uint32_t>> &core,
                                         const shared_ptr<unordered_map<uint32_t, shared_ptr<extend_list<double, uint32_t>>>> &k_order,
                                         const shared_ptr<unordered_map<uint32_t, uint32_t>> &rem,
                                         const shared_ptr<unordered_map<uint32_t, uint32_t>> &ts) {
        queue<uint32_t> Q;
        unordered_set<uint32_t> visited;
        list<uint32_t> removed_vertex_list;
        for (const auto &e:*ER) {
            auto u = e->get_source_vertex_id();
            auto v = e->get_destination_vertex_id();
            auto k = min(core->at(u), core->at(v));
            auto u_vertex = G->get_vertex(u);
            for (const auto&[w, e]:*u_vertex->get_edge_map()) {
                if (core->at(w) >= k) {
                    --ts->at(w);
                    if (ts->at(w) < k && !visited.count(w)) {
                        Q.push(w);
                        visited.insert(w);
                    }
                }
            }
            for (const auto&[w, e]:*G->get_vertex(v)->get_edge_map()) {
                if (core->at(w) >= k) {
                    --ts->at(w);
                    if (ts->at(w) < k && !visited.count(w)) {
                        Q.push(w);
                        visited.insert(w);
                    }
                }
            }
            G->remove_edge(e);
        }
        while (!Q.empty()) {
            auto u = Q.front();
            Q.pop();


            auto k_new = recompute_core_number(G, core, u);
            if (k_new == 0) {
                core->erase(u);
                rem->erase(u);
            }
            auto u_vertex = G->get_vertex(u);
            for (const auto&[v, e]:*u_vertex->get_edge_map()) {
                if (core->at(v) <= core->at(u)) {
                    --ts->at(v);
                    if (ts->at(v) < core->at(v) && !visited.count(v)) {
                        Q.push(v);
                        visited.insert(v);
                    }
                }
            }
            auto u_node = k_order->at(core->at(u))->remove(u);
            k_order->at(k_new)->right_insert(u_node);

            core->at(u) = k_new;
            rem->at(u) = 0;
            ts->at(u) = 0;
            for (const auto&[v, e]:*u_vertex->get_edge_map()) {
                if (core->at(v) >= core->at(u)) {
                    ++ts->at(u);
                }
                if (test_order(core, k_order, u, v)) {
                    ++rem->at(u);
                } else {
                    --rem->at(v);
                }
            }
        }
    }

    void order_maintenance::partial_core_decomposition(const shared_ptr<abstract_graph>& G,
                                                       const shared_ptr<unordered_set<uint32_t>>& Vc,
                                                       const shared_ptr<unordered_set<uint32_t>>& evicted_set,
                                                       const shared_ptr<unordered_map<uint32_t,uint32_t>>& vertex_core_map,
                                                       const shared_ptr<unordered_map<uint32_t,shared_ptr<unordered_map<uint32_t,uint32_t>>>>& CD,
                                                       const shared_ptr<unordered_map<uint32_t,shared_ptr<extend_list<double,uint32_t>>>>& k_order,
                                                       const shared_ptr<unordered_map<uint32_t,shared_ptr<extend_node<double,uint32_t>>>>& node_map,
                                                       uint32_t k)
    {
        //remove_unsatisfied_vertices(G,Vc,evicted_set,k,CD);
        while(!Vc->empty())
        {
            auto CDk = CD->at(k);
            k_order->insert({k, make_shared<extend_list<double,uint32_t>>()});
            CD->insert({k + 1, make_shared<unordered_map<uint32_t,uint32_t>>()});
            auto CD_next_k = CD->at(k + 1);
            for(auto u:*Vc) {
                vertex_core_map->at(u) = k;

                if(CDk->at(u) < k + 1)
                {
                    evicted_set->insert(u);
                } else
                {
                    CD_next_k->insert({u, CDk->at(u)});
                }
            }
            ++k;
            remove_unsatisfied_vertices(G, Vc, evicted_set, CD, k_order, node_map, k);
        }
    }

    void order_maintenance::remove_unsatisfied_vertices(const shared_ptr<abstract_graph> &G,
                                                        const shared_ptr<unordered_set<uint32_t>> &Vc,
                                                        const shared_ptr<unordered_set<uint32_t>> &evictedSet,
                                                        const shared_ptr<unordered_map<uint32_t, shared_ptr<unordered_map<uint32_t, uint32_t>>>> &kMcd,
                                                        const shared_ptr<unordered_map<uint32_t, shared_ptr<extend_list<double, uint32_t>>>> &kOrder,
                                                        const shared_ptr<unordered_map<uint32_t, shared_ptr<extend_node<double, uint32_t>>>> &nodeMap,
                                                        uint32_t k) {

        auto mcd = kMcd->at(k);
        while (!evictedSet->empty()) {
            auto u = *evictedSet->begin();
            evictedSet->erase(u);
            Vc->erase(u);

            auto u_vertex = G->get_vertex(u);
            for (const auto& [v,e]:*u_vertex->get_edge_map()) {
                if(!mcd->count(v)) {
                    continue;
                }
                --kMcd->at(k)->at(v);
                if (mcd->at(v)< k) {
                    evictedSet->insert(v);
                }
            }
            mcd->erase(u);
            kOrder->at(k-1)->left_insert(nodeMap->at(u));
        }
    }

    void order_maintenance::search_influenced_vertices(const shared_ptr<abstract_graph> &G,
                                                       const shared_ptr<unordered_set<shared_ptr<abstract_edge>>> &edge_set,
                                                       const shared_ptr<unordered_map<uint32_t, uint32_t>> &core,
                                                       const shared_ptr<unordered_map<uint32_t, shared_ptr<extend_list<double, uint32_t>>>> &k_order,
                                                       const shared_ptr<unordered_map<uint32_t, shared_ptr<unordered_map<uint32_t, uint32_t>>>> &CD,
                                                       const shared_ptr<unordered_set<uint32_t>> &Ck,
                                                       const shared_ptr<unordered_set<uint32_t>> &evicted_set,
                                                       const shared_ptr<extend_node<double,uint32_t>> &old_head,
                                                       uint32_t k){
        auto skip_set = make_shared<set<shared_ptr<extend_node<double,uint32_t>>,extend_node_compare>>();
        auto vertex_set = make_shared<unordered_set<uint32_t>>();
        auto CD_previous_k = CD->at(k - 1);
        auto CDk = CD->at(k);
        auto order_list = k_order->at(k - 1);

        auto iter = edge_set->begin();
        while (iter != edge_set->end()) {
            auto e = *iter;
            ++iter;
            auto u = e->get_source_vertex_id();
            auto v = e->get_destination_vertex_id();

            if (evicted_set->count(u) || !CD_previous_k->count(u)||CD_previous_k->at(u)<k) {
                edge_set->erase(e);
                continue;
            }

            if (evicted_set->count(v) || !CD_previous_k->count(v)||CD_previous_k->at(v)<k) {
                edge_set->erase(e);
                continue;
            }


            if (core->at(u) >= k && core->at(v) >= k) {
                ++CDk->at(u);
                ++CDk->at(v);
                continue;
            }

            if (order_list->count_value(u)) {
                auto u_node = order_list->find(u);
                skip_set->insert(u_node);
            }
            if (order_list->count_value(v)) {
                auto v_node = order_list->find(v);
                skip_set->insert(v_node);
            }
        }

        for(auto u:*Ck)
        {
            auto u_vertex = G->get_vertex(u);
            for(const auto &[v,e]:*u_vertex->get_edge_map())
            {
                if(!CD_previous_k->count(v) || evicted_set->count(v) || Ck->count(v))
                {
                    continue;
                }
                if(core->at(v)>=k) {
                    ++CDk->at(v);
                    continue;
                }
                if(order_list->count_value(v))
                {
                    skip_set->insert(order_list->find(v));
                }
            }
        }

        unordered_set<uint32_t> affected_set;
        while(!skip_set->empty())
        {
            affected_set.clear();
            auto u_node = *skip_set->begin();
            skip_set->erase(u_node);
            auto u = u_node->get_value();
            if(Ck->count(u))
            {
                continue;
            }
            CDk->insert({u, 0});
            Ck->insert(u);
            auto u_vertex = G->get_vertex(u);
            for (const auto &[v,e]:*u_vertex->get_edge_map()) {
                if(!CD_previous_k->count(v) || evicted_set->count(v) || CD_previous_k->at(v) < k)
                {
                    continue;
                }

                if(core->at(v)>=k)
                {
                    ++CDk->at(u);
                    ++CDk->at(v);
                    continue;
                }
                if (Ck->count(v)) {
                    ++CDk->at(u);
                    continue;
                }

                if (test_order(core, k_order, u, v)) {
                    ++CDk->at(u);
                    affected_set.insert(v);
                }
            }
            if (CDk->at(u) >= k) {
                for (auto v:affected_set) {
                    skip_set->insert( order_list->find(v));
                }
            } else {
                remove_unsatisfied_vertices(G, Ck, u, k, CD, k_order,evicted_set);
            }
        }
    }

    void order_maintenance::search_influenced_vertices(const shared_ptr<abstract_graph>& G,
                                                       const shared_ptr<unordered_set<shared_ptr<abstract_edge>>>& edge_set,
                                                       const shared_ptr<unordered_map<uint32_t,uint32_t>>& core,
                                                       const shared_ptr<unordered_map<uint32_t, shared_ptr<unordered_map<uint32_t ,uint32_t>>>>& CD,
                                                       const shared_ptr<unordered_map<uint32_t, shared_ptr<extend_list<double,uint32_t>>>>& k_order,
                                                       const shared_ptr<unordered_map<uint32_t ,shared_ptr<extend_node<double,uint32_t>>>>& node_map,
                                                       const shared_ptr<unordered_set<uint32_t>>& Vc,
                                                       const shared_ptr<unordered_set<uint32_t>>& evicted_set,
                                                       uint32_t k){
        auto B = make_shared<set<shared_ptr<extend_node<double,uint32_t>>,extend_node_compare>>();
        auto CDk = CD->at(k - 1);
        auto CD_next_k = CD->at(k);
        auto order_list = k_order->at(k - 1);

        auto iter = edge_set->begin();
        while (iter != edge_set->end()) {
            auto e = *iter;
            ++iter;
            auto u = e->get_source_vertex_id();
            auto v = e->get_destination_vertex_id();

            if (evicted_set->count(u) || CDk->at(u) < k) {
                edge_set->erase(e);
                continue;
            }

            if (evicted_set->count(v) || CDk->at(v) < k) {
                edge_set->erase(e);
                continue;
            }


            if (core->at(u) >= k && core->at(v) >= k) {
                ++CD_next_k->at(u);
                ++CD_next_k->at(v);
                continue;
            }

            if (order_list->count_value(u)) {
                B->insert(order_list->find(u));
            }
            if (order_list->count_value(v)) {
                B->insert(order_list->find(v));
            }
        }

        for(auto u:*Vc)
        {
            auto u_vertex = G->get_vertex(u);
            for(const auto &[v,e]:*u_vertex->get_edge_map())
            {
                if(!CDk->count(v) || evicted_set->count(v) || Vc->count(v) || CDk->at(v) < k)
                {
                    continue;
                }
                if(core->at(v) >= k) {
                    ++CD_next_k->at(v);
                }
                else
                {
                    B->insert(order_list->find(v));
                }
            }
        }

        unordered_set<uint32_t> affected_set;
        while(!B->empty())
        {
            affected_set.clear();
            auto u_node = *B->begin();
            B->erase(u_node);
            auto u = u_node->get_value();
            CD_next_k->insert({u, 0});
            Vc->insert(u);
            auto u_vertex = G->get_vertex(u);
            for (const auto &[v,e]:*u_vertex->get_edge_map()) {
                if(!CDk->count(v) || evicted_set->count(v) || CDk->at(v) < k)
                {
                    continue;
                }

                if(core->at(v) >= k)
                {
                    ++CD_next_k->at(u);
                    ++CD_next_k->at(v);
                } else if (Vc->count(v)) {
                    ++CD_next_k->at(u);
                } else if (test_order(core, k_order, u, v)) {
                    ++CD_next_k->at(u);
                    affected_set.insert(v);
                }
            }
            if (CD_next_k->at(u) >= k) {
                for (auto v:affected_set) {
                    B->insert(order_list->find(v));
                }
            } else {
                remove_unsatisfied_vertices(G, Vc, CD, k_order, node_map,
                                          evicted_set, u, k);
            }
        }
    }

    uint32_t order_maintenance::get_core_degree(const shared_ptr<abstract_graph> &G,
                                    const shared_ptr<unordered_map<uint32_t, uint32_t>> &core,
                                    uint32_t k,
                                    uint32_t u)
    {
        auto core_degree = 0;
        auto u_vertex = G->get_vertex(u);
        for(const auto &[v,e]:*u_vertex->get_edge_map())
        {
            if(core->at(v)>=k)
            {
                ++core_degree;
            }
        }
        return core_degree;
    }

    /**
     * @details initialize a core-order map and remaining degree for each vertex
     * @param G the given graph
     * @param core a map of vertices and their core number
     * @param k_order a map of core numbers and corresponding lists
     * @param rem a map of vertices and its remaining degree (the number of neighbors whose core-order is greater than it)
     * @param ext a map of vertices and candidate degree (the number of neighbors whose core-order is less than it)
     */
    void order_maintenance::init_order(const shared_ptr<abstract_graph> &G,
                                       const shared_ptr<unordered_map<uint32_t, uint32_t>> &core,
                                       const shared_ptr<unordered_map<uint32_t, shared_ptr<extend_list<double, uint32_t>>>> &k_order,
                                       const shared_ptr<unordered_map<uint32_t, uint32_t>> &rem,
                                       const shared_ptr<unordered_map<uint32_t, uint32_t>> &ext) {

        for (auto[u, core_number]:*core) {
            auto u_vertex = G->get_vertex(u);
            auto remaining_degree = 0;
            for (const auto &[v, e]:*u_vertex->get_edge_map()) {
                if (test_order(core, k_order, u, v)) {
                    ++remaining_degree;
                }
            }
            rem->insert({u, remaining_degree});
            ext->insert({u, 0});
        }
    }

    /**
     * @details recompute the core number of a given vertex
     * @param G the given graph
     * @param core a map of vertices and core numbers
     * @param u a given vertex
     * @return
     */
    uint32_t order_maintenance::recompute_core_number(const std::shared_ptr<abstract_graph> &G,
                                                      const std::shared_ptr<unordered_map<uint32_t, uint32_t>> &core,
                                                      uint32_t u) {
        auto u_vertex = G->get_vertex(u);
        if (!u_vertex) {
            return 0;
        }

        vector<uint32_t> core_number_vector(core->at(u) + 1, 0);
        for (const auto&[v, e]:*u_vertex->get_edge_map()) {
            if (core->at(v) >= core->at(u)) {
                ++core_number_vector.at(core->at(u));
            } else {
                ++core_number_vector.at(core->at(v));
            }
        }

        auto k = 0;
        auto sum = 0;
        for (auto index = core_number_vector.size() - 1; index >= 0; --index) {
            auto count = core_number_vector.at(index);
            sum += count;
            if (sum >= index) {
                k = index;
                break;
            }
        }
        return k;
    }

    /**
     * @details test core-order of two given vertices
     * @param core a map of vertices and their core numbers
     * @param k_order a map of core numbers and corresponding core-order lists
     * @param u the given vertex
     * @param v the given vertex
     * @return boolean value
     */
    bool order_maintenance::test_order(const shared_ptr<unordered_map<uint32_t, uint32_t>> &core,
                                       const shared_ptr<unordered_map<uint32_t, shared_ptr<extend_list<double, uint32_t>>>> &k_order,
                                       uint32_t u,
                                       uint32_t v) {
        auto k1 = core->at(u);
        auto k2 = core->at(v);
        if (k1 < k2) {
            return true;
        } else if (k1 > k2) {
            return false;
        } else {
            /**
             * @remarks the type of find results is optional
             */
            return k_order->at(k1)->find_key(u).value() < k_order->at(k2)->find_key(v).value();
        }
    }

    shared_ptr<extend_node<double,uint32_t>> order_maintenance::partial_core(const shared_ptr<abstract_graph> &G,
                                                                             const shared_ptr<unordered_map<uint32_t, uint32_t>> &core,
                                                                             const shared_ptr<unordered_map<uint32_t, shared_ptr<extend_list<double, uint32_t>>>> &k_order,
                                                                             const shared_ptr<unordered_map<uint32_t, shared_ptr<unordered_map<uint32_t, uint32_t>>>> &CD,
                                                                             const shared_ptr<unordered_set<uint32_t>> &Ck,
                                                                             const shared_ptr<unordered_set<uint32_t>> &evicted_set,
                                                                             uint32_t k) {

        unordered_set<uint32_t> vertex_set;
        auto CDk = CD->at(k);
        if(!CD->count(k+1))
        {
            CD->insert({k + 1, make_shared<unordered_map<uint32_t,uint32_t>>()});
            k_order->insert({k+1,make_shared<extend_list<double,uint32_t>>()});
        }

        auto CD_next_k = CD->at(k + 1);
        for (auto u:*Ck) {
            core->at(u) = k;
            CD_next_k->insert({u, CDk->at(u)});
            if (CD_next_k->at(u) == k) {
                vertex_set.insert(u);
            }
        }

        auto old_head = k_order->at(k)->get_head();
        while (!vertex_set.empty()) {
            auto u = *vertex_set.begin();
            vertex_set.erase(u);
            evicted_set->insert(u);

            auto u_vertex = G->get_vertex(u);
            for (const auto& [v,e]:*u_vertex->get_edge_map()) {
                if (!Ck->count(v)) {
                    continue;
                }
                --CD_next_k->at(v);
                if (CD_next_k->at(v) <= k) {
                    vertex_set.insert(v);
                }
            }
            Ck->erase(u);
            CD_next_k->erase(u);
            auto u_node = make_shared<extend_node<double,uint32_t>>(0,u);
            k_order->at(k)->insert_before(u_node,old_head);
        }
        return old_head;
    }

    void order_maintenance::partial_core(const shared_ptr<abstract_graph>& G,
                                         const shared_ptr<unordered_set<uint32_t>>& Vc,
                                         const shared_ptr<unordered_map<uint32_t,uint32_t>>& core,
                                         const shared_ptr<unordered_map<uint32_t,shared_ptr<extend_list<double,uint32_t>>>>& k_order,
                                         const shared_ptr<unordered_map<uint32_t,shared_ptr<extend_node<double,uint32_t>>>>& node_map,
                                         const shared_ptr<unordered_map<uint32_t,shared_ptr<unordered_map<uint32_t,uint32_t>>>>& CD,
                                         const shared_ptr<unordered_set<uint32_t>>& evicted_set,
                                         uint32_t k) {

        unordered_set<uint32_t> vertex_set;
        auto CDk = CD->at(k);
        CD->insert({k + 1, make_shared<unordered_map<uint32_t, uint32_t>>()});
        auto CD_next_k = CD->at(k + 1);
        for (auto u:*Vc) {
            core->at(u) = k;
            CD_next_k->insert({u, CDk->at(u)});
            if (CD_next_k->at(u) == k) {
                vertex_set.insert(u);
            }
        }

        auto old_head = k_order->at(k)->get_head();
        while (!vertex_set.empty()) {
            auto u = *vertex_set.begin();
            vertex_set.erase(u);
            evicted_set->insert(u);

            auto u_vertex = G->get_vertex(u);
            for (const auto& [v,e]:*u_vertex->get_edge_map()) {
                if (!Vc->count(v)) {
                    continue;
                }
                --CD_next_k->at(v);
                if (CD_next_k->at(v) <= k) {
                    vertex_set.insert(v);
                }
            }
            Vc->erase(u);
            CD_next_k->erase(u);
            k_order->at(k)->insert_before(node_map->at(u), old_head);
        }
    }

    shared_ptr<extend_node<double,uint32_t>> order_maintenance::partial_core(const shared_ptr<abstract_graph> &G,
                                                                             const shared_ptr<unordered_map<uint32_t, uint32_t>> &core,
                                                                             const shared_ptr<unordered_map<uint32_t, shared_ptr<extend_list<double, uint32_t>>>> &k_order,
                                                                             const shared_ptr<unordered_set<uint32_t>> &Ck,
                                                                             const shared_ptr<unordered_map<uint32_t, uint32_t>> &s,
                                                                             uint32_t k)
    {
        if(!k_order->count(k))
        {
            k_order->insert({k,make_shared<extend_list<double,uint32_t>>()});
        }
        auto old_head_node = k_order->at(k)->get_head();
        queue<uint32_t> Q;
        for (const auto &u:*Ck) {
            core->at(u) = k;

            auto u_vertex = G->get_vertex(u);
            if (s->at(u) <= k) {
                Q.push(u);
            }
        }

        while (!Q.empty()) {
            auto u = Q.front();
            Q.pop();

            auto u_vertex = G->get_vertex(u);
            for (const auto &[v, e]:*u_vertex->get_edge_map()) {
                if (Ck->count(v)) {
                    --s->at(v);
                    if (s->at(v) <= k) {
                        Q.push(v);
                    }
                }
            }
            Ck->erase(u);
            s->erase(u);

            auto u_node = make_shared<extend_node<double, uint32_t>>(0, u);
            k_order->at(k)->insert_before(u_node, old_head_node);
        }
        return old_head_node;
    }


    void order_maintenance::remove_unsatisfied_vertices(const shared_ptr<abstract_graph>& G,
                                                        const shared_ptr<unordered_set<uint32_t>>& Vc,
                                                        const shared_ptr<unordered_map<uint32_t,shared_ptr<unordered_map<uint32_t,uint32_t>>>>& CD,
                                                        const shared_ptr<unordered_map<uint32_t,shared_ptr<extend_list<double,uint32_t>>>>& k_order,
                                                        const shared_ptr<unordered_map<uint32_t,shared_ptr<extend_node<double,uint32_t>>>>& node_map,
                                                        const shared_ptr<unordered_set<uint32_t>>& evicted_set,
                                                        uint32_t w,
                                                        uint32_t k) {
        unordered_set<uint32_t> vertex_set;
        vertex_set.insert(w);
        auto w_node_next= k_order->at(k - 1)->find(w)->get_next();
        auto CDk = CD->at(k);
        while (!vertex_set.empty()) {
            auto u = *vertex_set.begin();
            vertex_set.erase(u);
            evicted_set->insert(u);

            auto u_vertex = G->get_vertex(u);
            for (const auto& [v,e]:*u_vertex->get_edge_map()) {
                if(!CDk->count(v)) {
                    continue;
                }
                --CDk->at(v);
                if (CDk->at(v) < k) {
                    vertex_set.insert(v);
                }
            }
            Vc->erase(u);
            k_order->at(k-1)->insert_before(node_map->at(u), w_node_next);
            CDk->erase(u);
        }
    }

    void order_maintenance::remove_unsatisfied_vertices(const shared_ptr<abstract_graph>& G,
                                                        const shared_ptr<unordered_set<uint32_t>>& Ck,
                                                        uint32_t w,
                                                        uint32_t k,
                                                        const shared_ptr<unordered_map<uint32_t,shared_ptr<unordered_map<uint32_t ,uint32_t>>>>& CD,
                                                        const shared_ptr<unordered_map<uint32_t ,shared_ptr<extend_list<double,uint32_t>>>>& k_order,
                                                        const shared_ptr<unordered_set<uint32_t>>& evicted_set) {
        unordered_set<uint32_t> vertex_set;
        vertex_set.insert(w);
        auto w_next= k_order->at(k - 1)->find(w)->get_next();
        auto CDk = CD->at(k);
        while (!vertex_set.empty()) {
            auto u = *vertex_set.begin();
            vertex_set.erase(u);
            evicted_set->insert(u);

            auto u_vertex = G->get_vertex(u);
            for (const auto &[v,e]:*u_vertex->get_edge_map()) {
                if(!CDk->count(v)) {
                    continue;
                }
                --CDk->at(v);
                if (CDk->at(v) < k) {
                    vertex_set.insert(v);
                }
            }
            Ck->erase(u);
            CDk->erase(u);
            auto u_node = make_shared<extend_node<double,uint32_t>>(0,u);
            k_order->at(k - 1)->insert_before(u_node, w_next);
        }
    }

}