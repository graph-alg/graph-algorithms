
#include "wing/quasi_order_wing_maintenance.h"


namespace scnu {

    void quasi_order_wing_maintenance::candidate_graph_finding(const shared_ptr<abstract_bipartite_graph> &G,
                                                         const shared_ptr<unordered_set<shared_ptr<abstract_bipartite_edge>>> &edge_set,
                                                         const shared_ptr<unordered_map<shared_ptr<abstract_bipartite_edge>, uint32_t>> &edge_wing_map,
                                                         const shared_ptr<unordered_map<uint32_t, shared_ptr<extend_list<uint32_t, shared_ptr<abstract_bipartite_edge>>>>> &wing_order_map,
                                                         const shared_ptr<unordered_map<shared_ptr<abstract_bipartite_edge>, uint32_t>> &WS_previous_k,
                                                         const shared_ptr<unordered_map<shared_ptr<abstract_bipartite_edge>, uint32_t>> &WS_k,
                                                         const shared_ptr<unordered_set<shared_ptr<abstract_bipartite_edge>>> &F_k,
                                                         const shared_ptr<unordered_set<shared_ptr<abstract_bipartite_edge>>> &evicted_edge_set,
                                                         uint32_t k) {
        auto affected_edge_set = make_shared<unordered_set<shared_ptr<abstract_bipartite_edge>>>();
        /**
         * @brief find the k-insert graph
         */
        auto S_k = make_shared<unordered_set<shared_ptr<abstract_bipartite_edge>>>();
        for (auto iter = edge_set->begin(); iter != edge_set->end();) {
            const auto &e = *iter;
            ++iter;
            if (edge_wing_map->at(e) == k - 1 && WS_previous_k->at(e) >= k) {
                S_k->insert(e);
            } else {
                edge_set->erase(e);
            }
        }

        /**
         * @brief find a set of candidate edges
         */
        auto wing_order = wing_order_map->at(k - 1);
        while (!S_k->empty()) {
            auto e1 = *S_k->begin();
            S_k->erase(e1);

            auto e1_support = 0;
            affected_edge_set->clear();

            auto l1 = e1->get_left_vertex_id();
            auto r1 = e1->get_right_vertex_id();
            auto l1_vertex = G->get_left_vertex(l1);
            auto r1_vertex = G->get_right_vertex(r1);

            for (const auto&[r2, e2]:*l1_vertex->get_edge_map()) {
                if (r2 == r1 || !WS_previous_k->count(e2) || WS_previous_k->at(e1) < k
                    || wing_order->find_key(e2) < wing_order->find_key(e1)) {
                    continue;
                }
                for (const auto&[l2, e3]:*r1_vertex->get_edge_map()) {
                    if (l2 == l1 || !WS_previous_k->count(e3) || WS_previous_k->at(e3) < k
                        || wing_order->find_key(e3) < wing_order->find_key(e1)) {
                        continue;
                    }
                    auto e4 = G->get_edge(l2, r2);
                    if (!e4 || !WS_previous_k->count(e4) || WS_previous_k->at(e4) < k
                        || wing_order->find_key(e4) < wing_order->find_key(e1)) {
                        continue;
                    }

                    ++e1_support;
                    if (edge_wing_map->at(e2) < k) {
                        affected_edge_set->insert(e2);
                    } else {
                        ++WS_k->at(e2);
                    }

                    if (edge_wing_map->at(e3) < k) {
                        affected_edge_set->insert(e3);
                        ++WS_k->at(e3);
                    }

                    if (edge_wing_map->at(e4) < k) {
                        affected_edge_set->insert(e4);
                        ++WS_k->at(e4);
                    }
                }
            }

            if (e1_support >= k) {
                F_k->insert(e1);
                WS_k->insert({e1, e1_support});
                for (const auto &e:*affected_edge_set) {
                    S_k->insert(e);
                }
            } else {
                evicted_edge_set->insert(e1);
            }
        }

        remove_unsatisfied_edges(G,WS_k,F_k,evicted_edge_set,k);
    }

    void quasi_order_wing_maintenance::candidate_graph_finding(const shared_ptr<abstract_bipartite_graph> &G,
                                                         const shared_ptr<unordered_map<shared_ptr<abstract_bipartite_edge>,shared_ptr<mutex>>> &edge_mutex_map,
                                                         const shared_ptr<unordered_set<shared_ptr<abstract_bipartite_edge>>> &edge_set,
                                                         const shared_ptr<unordered_map<shared_ptr<abstract_bipartite_edge>, uint32_t>> &edge_wing_map,
                                                         const shared_ptr<unordered_map<uint32_t, shared_ptr<extend_list<uint32_t, shared_ptr<abstract_bipartite_edge>>>>> &wing_order_map,
                                                         const shared_ptr<unordered_map<shared_ptr<abstract_bipartite_edge>, uint32_t>> &WS_previous_k,
                                                         const shared_ptr<unordered_map<shared_ptr<abstract_bipartite_edge>, uint32_t>> &WS_k,
                                                         const shared_ptr<unordered_set<shared_ptr<abstract_bipartite_edge>>> &F_k,
                                                         const shared_ptr<mutex> &F_k_mutex,
                                                         const shared_ptr<unordered_set<shared_ptr<abstract_bipartite_edge>>> &evicted_edge_set,
                                                         const shared_ptr<mutex> &evicted_edge_set_mutex,
                                                         uint32_t k,
                                                         uint32_t thread_count) {
        /**
         * @brief find the k-insert graph
         */
        auto S_k = make_shared<unordered_set<shared_ptr<abstract_bipartite_edge>>>();
        for (auto iter = edge_set->begin(); iter != edge_set->end();) {
            const auto &e = *iter;
            ++iter;
            if (edge_wing_map->at(e) == k - 1 && WS_previous_k->at(e) >= k) {
                S_k->insert(e);
            } else {
                edge_set->erase(e);
            }
        }

        /**
         * @brief find a set of candidate edges
         */
        auto task_vector = make_shared<vector<shared_ptr<unordered_set<shared_ptr<abstract_bipartite_edge>>>>>();
        for(uint32_t i = 0;i<thread_count;++i){
            task_vector->emplace_back(make_shared<unordered_set<shared_ptr<abstract_bipartite_edge>>>());
        }
        thread_pool pool(thread_count);
        auto current_order_list = wing_order_map->at(k - 1);
        while (!S_k->empty()) {
            uint32_t count = 0;
            for(const auto&e:*S_k){
                uint32_t index = count % thread_count;
                ++count;
                task_vector->at(index)->insert(e);
            }
            S_k->clear();
            for(uint32_t i = 0;i<thread_count;++i){
                auto sub_edge_set = task_vector->at(i);
                auto task = [=]{
                    while (!sub_edge_set->empty()){
                        auto e1 = *sub_edge_set->begin();
                        sub_edge_set->erase(e1);

                        auto e1_support = 0;
                        auto neighbor_edge_set = make_shared<unordered_set<shared_ptr<abstract_bipartite_edge>>>();

                        auto l1 = e1->get_left_vertex_id();
                        auto r1 = e1->get_right_vertex_id();
                        auto l1_vertex = G->get_left_vertex(l1);
                        auto r1_vertex = G->get_right_vertex(r1);

                        for (const auto&[r2, e2]:*l1_vertex->get_edge_map()) {
                            if (r2 == r1 || !WS_previous_k->count(e2) || WS_previous_k->at(e2) < k
                                || current_order_list->find_key(e2) < current_order_list->find_key(e1)) {
                                continue;
                            }
                            for (const auto&[l2, e3]:*r1_vertex->get_edge_map()) {
                                if (l2 == l1 || !WS_previous_k->count(e3) || WS_previous_k->at(e3) < k
                                    || current_order_list->find_key(e3) < current_order_list->find_key(e1)) {
                                    continue;
                                }
                                auto e4 = G->get_edge(l2, r2);
                                if (!e4 || !WS_previous_k->count(e4) || WS_previous_k->at(e4) < k
                                    || current_order_list->find_key(e4) < current_order_list->find_key(e1)) {
                                    continue;
                                }

                                ++e1_support;

                                if (edge_wing_map->at(e2) < k) {
                                    neighbor_edge_set->insert(e2);
                                } else {
                                    edge_mutex_map->at(e2)->lock();
                                    ++WS_k->at(e2);
                                    edge_mutex_map->at(e2)->unlock();
                                }

                                if (edge_wing_map->at(e3) < k) {
                                    neighbor_edge_set->insert(e3);
                                }else{
                                    edge_mutex_map->at(e3)->lock();
                                    ++WS_k->at(e2);
                                    edge_mutex_map->at(e3)->unlock();
                                }

                                if (edge_wing_map->at(e4) < k) {
                                    neighbor_edge_set->insert(e4);
                                }else
                                {
                                    edge_mutex_map->at(e4)->lock();
                                    ++WS_k->at(e4);
                                    edge_mutex_map->at(e4)->unlock();
                                }
                            }
                        }

                        if (e1_support >= k) {
                            F_k_mutex->lock();
                            F_k->insert(e1);
                            WS_k->insert({e1, e1_support});
                            copy(neighbor_edge_set->begin(),neighbor_edge_set->end(), inserter(*S_k,S_k->end()));
                            F_k_mutex->unlock();
                        } else {
                            evicted_edge_set_mutex->lock();
                            evicted_edge_set->insert(e1);
                            evicted_edge_set_mutex->unlock();
                        }
                    }
                };
                pool.submit_task(task);
            }
        }
        remove_unsatisfied_edges(G,WS_k,F_k,evicted_edge_set,k);
    }

    void quasi_order_wing_maintenance::find_affected_edge_set(
            const shared_ptr<abstract_bipartite_graph> &G,
            const shared_ptr<unordered_set<shared_ptr<abstract_bipartite_edge>>> &edge_set,
            const shared_ptr<unordered_map<uint32_t, shared_ptr<unordered_map<shared_ptr<abstract_bipartite_edge>, uint32_t>>>> &WS,
            const shared_ptr<unordered_set<shared_ptr<abstract_bipartite_edge>>> &affected_edge_set,
            uint32_t k) {
        auto WS_k = WS->at(k);
        for(const auto &e1:*edge_set){
            WS_k->erase(e1);
            auto l1 = e1->get_left_vertex_id();
            auto r1 = e1->get_right_vertex_id();
            for (const auto &[r2, e2]:*G->get_left_vertex(l1)->get_edge_map()) {
                if (r2 == r1 ||!WS_k->count(e2)) {
                    continue;
                }
                for (const auto &[l2, e3]:*G->get_right_vertex(r1)->get_edge_map()) {
                    if (l2 == l1 ||!WS_k->count(e3)) {
                        continue;
                    }

                    auto e4 = G->get_edge(l2, r2);
                    if (!e4  ||!WS_k->count(e4)) {
                        continue;
                    }

                    --WS_k->at(e2);
                    if (WS_k->at(e2) < k) {
                        affected_edge_set->insert(e2);
                    }

                    --WS_k->at(e3);
                    if (WS_k->at(e3) < k) {
                        affected_edge_set->insert(e3);
                    }

                    --WS_k->at(e4);
                    if (WS_k->at(e4) < k) {
                        affected_edge_set->insert(e4);
                    }
                }
            }
        }
    }

    void quasi_order_wing_maintenance::find_affected_edge_set(const shared_ptr<abstract_bipartite_graph> &G,
                                                        const shared_ptr<unordered_set<shared_ptr<abstract_bipartite_edge>>> &edge_set,
                                                        const shared_ptr<unordered_map<uint32_t, shared_ptr<unordered_map<shared_ptr<abstract_bipartite_edge>, uint32_t>>>> &WS,
                                                        const shared_ptr<unordered_map<shared_ptr<abstract_bipartite_edge>,shared_ptr<mutex>>>& edge_mutex_map,
                                                        const shared_ptr<unordered_set<shared_ptr<abstract_bipartite_edge>>>& affected_edge_set,
                                                        const shared_ptr<mutex>& affected_edge_set_mutex,
                                                        uint32_t k,
                                                        uint32_t thread_count) {
        auto WS_k = WS->at(k);
        auto task_vector = make_shared<vector<shared_ptr<unordered_set<shared_ptr<abstract_bipartite_edge>>>>>();
        uint32_t count = 0;
        for(const auto &e:*edge_set){
            auto index = count % thread_count;
            ++count;
            task_vector->at(index)->insert(e);
        }
        thread_pool pool(thread_count);
        for(uint32_t i = 0;i<thread_count;++i){
            auto task = [=]{
                auto sub_edge_set = task_vector->at(i);
                auto sub_affected_edge_set = make_shared<unordered_set<shared_ptr<abstract_bipartite_edge>>>();
                while(!sub_edge_set->empty()){
                    auto e1 = *sub_edge_set->begin();
                    sub_edge_set->erase(e1);

                    edge_mutex_map->at(e1)->lock();
                    WS_k->at(e1) = 0;
                    edge_mutex_map->at(e1)->unlock();

                    auto l1 = e1->get_left_vertex_id();
                    auto r1 = e1->get_right_vertex_id();
                    for (const auto &[r2, e2]:*G->get_left_vertex(l1)->get_edge_map()) {
                        if (r2 == r1 ||!WS_k->count(e2)) {
                            continue;
                        }
                        edge_mutex_map->at(e2)->lock();
                        if(WS_k->at(e2) < k){
                            edge_mutex_map->at(e2)->unlock();
                            continue;
                        }else
                        {
                            edge_mutex_map->at(e2)->unlock();
                        }
                        for (const auto &[l2, e3]:*G->get_right_vertex(r1)->get_edge_map()) {
                            if (l2 == l1 ||!WS_k->count(e3)) {
                                continue;
                            }
                            edge_mutex_map->at(e3)->lock();
                            if(WS_k->at(e3) < k){
                                edge_mutex_map->at(e3)->unlock();
                                continue;
                            }else
                            {
                                edge_mutex_map->at(e3)->unlock();
                            }
                            auto e4 = G->get_edge(l2, r2);
                            if (!e4  ||!WS_k->count(e4)) {
                                continue;
                            }
                            edge_mutex_map->at(e4)->lock();
                            if(WS_k->at(e4) < k){
                                edge_mutex_map->at(e4)->unlock();
                                continue;
                            }else
                            {
                                edge_mutex_map->at(e4)->unlock();
                            }

                            edge_mutex_map->at(e2)->lock();
                            --WS_k->at(e2);
                            if (WS_k->at(e2) < k) {
                                sub_affected_edge_set->insert(e2);
                            }
                            edge_mutex_map->at(e2)->unlock();

                            edge_mutex_map->at(e3)->lock();
                            --WS_k->at(e3);
                            if (WS_k->at(e3) < k) {
                                sub_affected_edge_set->insert(e3);
                            }
                            edge_mutex_map->at(e3)->unlock();

                            edge_mutex_map->at(e4)->lock();
                            --WS_k->at(e4);
                            if (WS_k->at(e4) < k) {
                                sub_affected_edge_set->insert(e4);
                            }
                            edge_mutex_map->at(e4)->unlock();
                        }
                    }
                }
                affected_edge_set_mutex->lock();
                copy(sub_affected_edge_set->begin(),sub_affected_edge_set->end(), inserter(*affected_edge_set,affected_edge_set->end()));
                affected_edge_set_mutex->unlock();
            };
            pool.submit_task(task);
        }
        for(const auto&e:*edge_set){
            WS_k->erase(e);
        }
    }


    /**
     * @details maintain wing number of a given graph under the insert case
     * @param G: the previous graph
     * @param edge_set: a set of inserted edges
     * @param edge_wing_map: a map of pairs of edges and their wing number
     * @param wing_order_map: the wing order of each k-wing
     * @param WS: the support of edges in each k-wing
     * @param previous_k_max: the maximal wing number of the previous graph
     */
    void quasi_order_wing_maintenance::insert(const shared_ptr<abstract_bipartite_graph> &G,
                                        const shared_ptr<unordered_set<shared_ptr<abstract_bipartite_edge>>> &edge_set,
                                        const shared_ptr<unordered_map<shared_ptr<abstract_bipartite_edge>, uint32_t>> &edge_wing_map,
                                        const shared_ptr<unordered_map<uint32_t, shared_ptr<extend_list<uint32_t, shared_ptr<abstract_bipartite_edge>>>>> &wing_order_map,
                                        const shared_ptr<unordered_map<uint32_t, shared_ptr<unordered_map<shared_ptr<abstract_bipartite_edge>, uint32_t>>>> &WS,
                                        const shared_ptr<uint32_t>& previous_k_max) {
        /**
         * @brief deal with 0-wing and 1-wing
         */
        init_insertion(G, edge_set, edge_wing_map, wing_order_map,WS);

        /**
         * @brief update the remainder k-wings
         */
        uint32_t k = 1;
        auto F_k = make_shared<unordered_set<shared_ptr<abstract_bipartite_edge>>>();
        auto evicted_edge_set = make_shared<unordered_set<shared_ptr<abstract_bipartite_edge>>>();
        while(!edge_set->empty()){
            if (k < *previous_k_max) {
                candidate_graph_finding(G,edge_set,edge_wing_map,wing_order_map,WS->at(k-1),WS->at(k),F_k,evicted_edge_set,k);
                partial_wing(F_k,edge_wing_map,wing_order_map,k);
            } else {
                candidate_graph_finding(G,edge_set,edge_wing_map,wing_order_map,WS->at(k-1),WS->at(k),F_k,evicted_edge_set,k);
                *previous_k_max = partial_wing_decomposition(G,F_k,edge_wing_map,wing_order_map,WS,k);
                break;
            }
            /**
              * @brief continue the loop
              */
            ++k;
            F_k->clear();
            evicted_edge_set->clear();
        }
    }

    void quasi_order_wing_maintenance::insert(const shared_ptr<abstract_bipartite_graph> &G,
                                        const shared_ptr<unordered_set<shared_ptr<abstract_bipartite_edge>>> &edge_set,
                                        const shared_ptr<unordered_map<shared_ptr<abstract_bipartite_edge>, uint32_t>> &edge_wing_map,
                                        const shared_ptr<unordered_map<uint32_t, shared_ptr<extend_list<uint32_t, shared_ptr<abstract_bipartite_edge>>>>> &wing_order_map,
                                        const shared_ptr<unordered_map<uint32_t, shared_ptr<unordered_map<shared_ptr<abstract_bipartite_edge>, uint32_t>>>> &WS,
                                        const shared_ptr<uint32_t>& previous_k_max,
                                        uint32_t thread_count) {
        /**
         * @brief deal with 0-wing and 1-wing
         */
        init_insertion(G, edge_set, edge_wing_map, wing_order_map,WS);

        /**
         * @brief update the remainder k-wings
         */
        uint32_t k = 1;
        auto F_k = make_shared<unordered_set<shared_ptr<abstract_bipartite_edge>>>();
        auto evicted_edge_set = make_shared<unordered_set<shared_ptr<abstract_bipartite_edge>>>();
        while(!edge_set->empty()){
            if (k < *previous_k_max) {
                candidate_graph_finding(G,edge_set,edge_wing_map,wing_order_map,WS->at(k-1),WS->at(k),F_k,evicted_edge_set,k);
                partial_wing(F_k,edge_wing_map,wing_order_map,k);
            } else {
                candidate_graph_finding(G,edge_set,edge_wing_map,wing_order_map,WS->at(k-1),WS->at(k),F_k,evicted_edge_set,k);
                *previous_k_max = partial_wing_decomposition(G,F_k,edge_wing_map,wing_order_map,WS,k);
                break;
            }
            /**
              * @brief continue the loop
              */
            ++k;
            F_k->clear();
            evicted_edge_set->clear();
        }
    }

    void quasi_order_wing_maintenance::partial_wing(const shared_ptr<unordered_set<shared_ptr<abstract_bipartite_edge>>> &F_k,
                                              const shared_ptr<unordered_map<shared_ptr<abstract_bipartite_edge>, uint32_t>> &edge_wing_map,
                                              const shared_ptr<unordered_map<uint32_t, shared_ptr<extend_list<uint32_t, shared_ptr<abstract_bipartite_edge>>>>> &wing_order_map,
                                              uint32_t k) {
        /**
            * @brief update wing number of edges in the candidate graph
            */
        if (!wing_order_map->count(k)) {
            wing_order_map->insert({k, make_shared<extend_list<uint32_t , shared_ptr<abstract_bipartite_edge>>>()});
        }
        for (const auto &e:*F_k) {
            auto e_node = wing_order_map->at(k - 1)->remove(e);
            edge_wing_map->at(e) = k;
            wing_order_map->at(k)->left_insert(e_node);
        }
    }

    /**
     * @details decompose the remainder k-wings
     * @param G
     * @param F_k
     * @param edge_wing_map
     * @param WS
     * @param wing_order_map
     * @param k
     */
    uint32_t quasi_order_wing_maintenance::partial_wing_decomposition(const shared_ptr<abstract_bipartite_graph> &G,
                                                                const shared_ptr<unordered_set<shared_ptr<abstract_bipartite_edge>>> &F_k,
                                                                const shared_ptr<unordered_map<shared_ptr<abstract_bipartite_edge>, uint32_t>> &edge_wing_map,
                                                                const shared_ptr<unordered_map<uint32_t, shared_ptr<extend_list<uint32_t, shared_ptr<abstract_bipartite_edge>>>>> &wing_order_map,
                                                                const shared_ptr<unordered_map<uint32_t,shared_ptr<unordered_map<shared_ptr<abstract_bipartite_edge>, uint32_t>>>> &WS,
                                                                uint32_t k) {
        uint32_t k_max = 0;
        auto evicted_set = make_shared<unordered_set<shared_ptr<abstract_bipartite_edge>>>();
        while(!F_k->empty()){
            k_max = k;
            if(!WS->count(k+1)){
                WS->insert({k+1, make_shared<unordered_map<shared_ptr<abstract_bipartite_edge>,uint32_t>>()});
            }
            auto WS_k = WS->at(k);
            auto WS_next_k = WS->at(k+1);
            for(auto iter = F_k->begin();iter!=F_k->end();){
                auto e = *iter;
                ++iter;
                edge_wing_map->at(e) = k;
                auto e_node = wing_order_map->at(k-1)->remove(e);
                wing_order_map->at(k)->left_insert(e_node);

                if (WS_k->at(e) < k + 1) {
                    evicted_set->insert(e);
                } else {
                    WS_next_k->insert({e, WS_k->at(e)});
                }
            }

            remove_unsatisfied_edges(G, WS_k, F_k, evicted_set, k + 1);
            ++k;
        }
        return k_max;
    }

    uint32_t quasi_order_wing_maintenance::partial_wing_decomposition(const shared_ptr<abstract_bipartite_graph> &G,
                                                                const shared_ptr<unordered_set<shared_ptr<abstract_bipartite_edge>>> &F_k,
                                                                const shared_ptr<unordered_map<shared_ptr<abstract_bipartite_edge>, uint32_t>> &edge_wing_map,
                                                                const shared_ptr<unordered_map<uint32_t, shared_ptr<extend_list<uint32_t, shared_ptr<abstract_bipartite_edge>>>>> &wing_order_map,
                                                                const shared_ptr<unordered_map<uint32_t,shared_ptr<unordered_map<shared_ptr<abstract_bipartite_edge>, uint32_t>>>> &WS,
                                                                const shared_ptr<unordered_map<shared_ptr<abstract_bipartite_edge>,shared_ptr<mutex>>>& edge_mutex_map,
                                                                uint32_t k,
                                                                uint32_t thread_count) {
        uint32_t k_max = 0;
        auto evicted_edge_set = make_shared<unordered_set<shared_ptr<abstract_bipartite_edge>>>();
        auto evicted_edge_set_mutex = make_shared<mutex>();
        while(!F_k->empty()){
            k_max = k;
            if(!WS->count(k+1)){
                WS->insert({k+1, make_shared<unordered_map<shared_ptr<abstract_bipartite_edge>,uint32_t>>()});
            }
            auto WS_k = WS->at(k);
            auto WS_next_k = WS->at(k+1);
            for(auto iter = F_k->begin();iter!=F_k->end();){
                auto e = *iter;
                ++iter;
                edge_wing_map->at(e) = k;
                auto e_node = wing_order_map->at(k-1)->remove(e);
                wing_order_map->at(k)->left_insert(e_node);

                if (WS_k->at(e) < k + 1) {
                    evicted_edge_set->insert(e);
                } else {
                    WS_next_k->insert({e, WS_k->at(e)});
                }
            }

            remove_unsatisfied_edges(G, WS_k, F_k, evicted_edge_set, evicted_edge_set_mutex, edge_mutex_map, k + 1,thread_count);
            ++k;
        }
        return k_max;
    }

    void quasi_order_wing_maintenance::remove(const shared_ptr<abstract_bipartite_graph> &G,
                                        const shared_ptr<unordered_set<shared_ptr<abstract_bipartite_edge>>> &edge_set,
                                        const shared_ptr<unordered_map<shared_ptr<abstract_bipartite_edge>, uint32_t>> &edge_wing_map,
                                        const shared_ptr<unordered_map<uint32_t, shared_ptr<extend_list<uint32_t, shared_ptr<abstract_bipartite_edge>>>>> &wing_order_map,
                                        const shared_ptr<unordered_map<uint32_t, shared_ptr<unordered_map<shared_ptr<abstract_bipartite_edge>, uint32_t>>>> &WS,
                                        const shared_ptr<uint32_t> &k_max) {
        auto removed_edge_set = container_copy::to_unordered_set<shared_ptr<abstract_bipartite_edge>>(edge_set);

        uint32_t k = 1;
        auto S_r = make_shared<unordered_set<shared_ptr<abstract_bipartite_edge>>>();
        while(!edge_set->empty()){
            auto WS_k = WS->at(k);
            for(auto iter = edge_set->begin();iter!=edge_set->end();){
                auto e = *iter;
                ++iter;
                if(!WS_k->count(e)){
                    WS_k->erase(e);
                    edge_set->erase(e);
                }else
                {
                    S_r->insert(e);
                }
            }

            update_single_wing(G,edge_set,edge_wing_map,wing_order_map,WS,k);
            if(WS_k->empty())
            {
                *k_max = k > 1 ? k - 1:0;
            }

            S_r->clear();
            ++k;
        }

        for(const auto&e:*removed_edge_set){
            G->remove_edge(e);
        }
    }

    void quasi_order_wing_maintenance::remove(const shared_ptr<abstract_bipartite_graph> &G,
                                        const shared_ptr<unordered_set<shared_ptr<abstract_bipartite_edge>>> &edge_set,
                                        const shared_ptr<unordered_map<shared_ptr<abstract_bipartite_edge>, uint32_t>> &edge_wing_map,
                                        const shared_ptr<unordered_map<uint32_t, shared_ptr<extend_list<uint32_t, shared_ptr<abstract_bipartite_edge>>>>> &wing_order_map,
                                        const shared_ptr<unordered_map<uint32_t, shared_ptr<unordered_map<shared_ptr<abstract_bipartite_edge>, uint32_t>>>> &WS,
                                        const shared_ptr<uint32_t> &k_max,
                                        uint32_t thread_count) {
        /**
         * @brief edge_support_computation data structure
         */
        auto edge_mutex_map = make_shared<unordered_map<shared_ptr<abstract_bipartite_edge>,shared_ptr<mutex>>>();
        for(const auto&e:*G->get_edge_set()){
            edge_mutex_map->insert({e, make_shared<mutex>()});
        }

        auto removed_edge_set = make_shared<unordered_set<shared_ptr<abstract_bipartite_edge>>>();

        auto evicted_edge_set = make_shared<unordered_set<shared_ptr<abstract_bipartite_edge>>>();
        auto evicted_edge_set_mutex = make_shared<mutex>();
        /**
         * @brief indicate k_max_flag is updated or not
         */
        auto update_flag = false;
        uint32_t k = 1;
        while(!edge_set->empty()){
            auto WS_k = WS->at(k);
            for(auto iter = edge_set->begin();iter!=edge_set->end();){
                auto e = *iter;
                ++iter;
                if(!WS_k->count(e)){
                    WS_k->erase(e);
                    edge_set->erase(e);
                    removed_edge_set->insert(e);
                }
            }
            find_affected_edge_set(G, edge_set, WS, edge_mutex_map,evicted_edge_set,evicted_edge_set_mutex, k,thread_count);
            update_single_wing(G, evicted_edge_set,evicted_edge_set_mutex, edge_wing_map, edge_mutex_map, wing_order_map, WS, k, thread_count);
            if(WS_k->empty() && !update_flag)
            {
                *k_max = k > 1 ? k - 1:0;
                update_flag = true;
            }
            ++k;
        }

        for(const auto&e:*removed_edge_set){
            G->remove_edge(e);
        }
    }

    void quasi_order_wing_maintenance::remove_unsatisfied_edges(const shared_ptr<abstract_bipartite_graph> &G,
                                                          const shared_ptr<unordered_map<shared_ptr<abstract_bipartite_edge>, uint32_t>> &WS_k,
                                                          const shared_ptr<unordered_set<shared_ptr<abstract_bipartite_edge>>> &F_k,
                                                          const shared_ptr<unordered_set<shared_ptr<abstract_bipartite_edge>>> &evicted_edge_set,
                                                          uint32_t k) {
        /**
            * @brief remove unsatisfied edges
            */
        while (!evicted_edge_set->empty()) {
            auto e1 = *evicted_edge_set->begin();
            evicted_edge_set->erase(e1);
            WS_k->erase(e1);
            F_k->erase(e1);

            auto l1 = e1->get_left_vertex_id();
            auto r1 = e1->get_right_vertex_id();
            for (const auto &[r2, e2]:*G->get_left_vertex(l1)->get_edge_map()) {
                if (r2 == r1 || WS_k->count(e2)) {
                    continue;
                }
                for (const auto &[l2, e3]:*G->get_right_vertex(r1)->get_edge_map()) {
                    if (l2 == l1 || WS_k->count(e3)) {
                        continue;
                    }

                    auto e4 = G->get_edge(l2, r2);
                    if (!e4 || WS_k->count(e4)) {
                        continue;
                    }

                    --WS_k->at(e2);
                    if (WS_k->at(e2) < k) {
                        evicted_edge_set->insert(e2);
                    }

                    --WS_k->at(e3);
                    if (WS_k->at(e3) < k) {
                        evicted_edge_set->insert(e3);
                    }

                    --WS_k->at(e4);
                    if (WS_k->at(e4) < k) {
                        evicted_edge_set->insert(e4);
                    }
                }
            }
        }
    }

    void quasi_order_wing_maintenance::remove_unsatisfied_edges(const shared_ptr<abstract_bipartite_graph> &G,
                                                          const shared_ptr<unordered_map<shared_ptr<abstract_bipartite_edge>, uint32_t>> &WS_k,
                                                          const shared_ptr<unordered_set<shared_ptr<abstract_bipartite_edge>>> &F_k,
                                                          const shared_ptr<unordered_set<shared_ptr<abstract_bipartite_edge>>> &evicted_edge_set,
                                                          const shared_ptr<mutex>& evicted_edge_set_mutex,
                                                          const shared_ptr<unordered_map<shared_ptr<abstract_bipartite_edge>,shared_ptr<mutex>>>& edge_mutex_map,
                                                          uint32_t k,
                                                          uint32_t thread_count) {
        /**
         * @brief remove unsatisfied edges
         */
        auto task_vector = make_shared<vector<shared_ptr<unordered_set<shared_ptr<abstract_bipartite_edge>>>>>();
        for(uint32_t i = 0;i<thread_count;++i){
            task_vector->emplace_back(make_shared<unordered_set<shared_ptr<abstract_bipartite_edge>>>());
        }
        thread_pool pool(thread_count);
        while (!evicted_edge_set->empty()) {
            uint32_t count = 0;
            for(const auto&e:*evicted_edge_set){
                task_vector->at(count++)->insert(e);
                WS_k->erase(e);
                F_k->erase(e);
            }
            evicted_edge_set->clear();
            for(uint32_t i = 0;i<thread_count;++i){
                auto task = [=]{
                    auto sub_edge_set = task_vector->at(i);
                    auto affected_edge_set = make_shared<unordered_set<shared_ptr<abstract_bipartite_edge>>>();
                    while (!sub_edge_set->empty()){
                        auto e1 = *sub_edge_set->begin();
                        sub_edge_set->erase(e1);

                        auto l1 = e1->get_left_vertex_id();
                        auto r1 = e1->get_right_vertex_id();
                        for (const auto &[r2, e2]:*G->get_left_vertex(l1)->get_edge_map()) {
                            if (r2 == r1 || WS_k->count(e2)) {
                                continue;
                            }

                            for (const auto &[l2, e3]:*G->get_right_vertex(r1)->get_edge_map()) {
                                if (l2 == l1 || WS_k->count(e3)) {
                                    continue;
                                }

                                auto e4 = G->get_edge(l2, r2);
                                if (!e4 || WS_k->count(e4)) {
                                    continue;
                                }

                                edge_mutex_map->at(e2)->lock();
                                --WS_k->at(e2);
                                if (WS_k->at(e2) < k) {
                                    affected_edge_set->insert(e2);
                                }
                                edge_mutex_map->at(e2)->unlock();

                                edge_mutex_map->at(e3)->lock();
                                --WS_k->at(e3);
                                if (WS_k->at(e3) < k) {
                                    affected_edge_set->insert(e3);
                                }
                                edge_mutex_map->at(e3)->unlock();

                                edge_mutex_map->at(e4)->lock();
                                --WS_k->at(e4);
                                if (WS_k->at(e4) < k) {
                                    affected_edge_set->insert(e4);
                                }
                                edge_mutex_map->at(e4)->unlock();
                            }
                        }
                    }

                    evicted_edge_set_mutex->lock();
                    copy(affected_edge_set->begin(),affected_edge_set->end(), inserter(*evicted_edge_set,evicted_edge_set->end()));
                    evicted_edge_set_mutex->unlock();
                };
                pool.submit_task(task);
            }
            pool.barrier();
        }

        remove_unsatisfied_edges(G,WS_k,F_k,evicted_edge_set,evicted_edge_set_mutex,edge_mutex_map,k,thread_count);
    }

    void quasi_order_wing_maintenance::init_insertion(const shared_ptr<abstract_bipartite_graph> &G,
                                                const shared_ptr<unordered_set<shared_ptr<abstract_bipartite_edge>>> &edge_set,
                                                const shared_ptr<unordered_map<shared_ptr<abstract_bipartite_edge>,uint32_t>>& edge_wing_map,
                                                const shared_ptr<unordered_map<uint32_t,shared_ptr<extend_list<uint32_t,shared_ptr<abstract_bipartite_edge>>>>>& wing_order_map,
                                                const shared_ptr<unordered_map<uint32_t, shared_ptr<unordered_map<shared_ptr<abstract_bipartite_edge>, uint32_t>>>> &WS) {

        if(!WS->count(0))
        {
            WS->insert({0, make_shared<unordered_map<shared_ptr<abstract_bipartite_edge>,uint32_t>>()});
        }
        auto WS0 = WS->at(0);

        if(!WS->count(1))
        {
            WS->insert({1, make_shared<unordered_map<shared_ptr<abstract_bipartite_edge>,uint32_t>>()});
        }
        auto WS1 = WS->at(1);

        if (!wing_order_map->count(0)) {
            wing_order_map->insert(
                    {0, make_shared<extend_list<uint32_t, shared_ptr<abstract_bipartite_edge>>>()});
            WS->insert({0, make_shared<unordered_map<shared_ptr<abstract_bipartite_edge>, uint32_t>>()});
        }
        auto order_list = wing_order_map->at(1);

        for(const auto&e1:*edge_set){
            WS1->insert({e1,0});
            G->insert_edge(e1);
            edge_wing_map->insert({e1,0});
            order_list->left_insert(e1);
        }

        auto visited_edge_set = make_shared<unordered_set<shared_ptr<abstract_bipartite_edge>>>();
        for(const auto&e1:*edge_set){
            auto l1 = e1->get_left_vertex_id();
            auto r1 = e1->get_right_vertex_id();

            visited_edge_set->insert(e1);
            for (const auto &[r2, e2]:*G->get_left_vertex(l1)->get_edge_map()) {
                if (r2 == r1 || visited_edge_set->count(e2)) {
                    continue;
                }
                for (const auto &[l2, e3]:*G->get_right_vertex(r1)->get_edge_map()) {
                    if (l2 == l1 || visited_edge_set->count(e3)) {
                        continue;
                    }

                    auto e4 = G->get_edge(l2, r2);
                    if (!e4 || visited_edge_set->count(e4)) {
                        continue;
                    }

                    ++WS0->at(e2);
                    ++WS1->at(e2);

                    ++WS0->at(e2);
                    ++WS1->at(e2);

                    ++WS0->at(e3);
                    ++WS1->at(e3);

                    ++WS0->at(e4);
                    ++WS1->at(e4);
                }
            }

            if(WS1->at(e1) == 0){
                WS1->erase(e1);
            }else
            {
                edge_wing_map->at(e1) = 1;
            }
        }
    }

    void quasi_order_wing_maintenance::init_insertion(const shared_ptr<abstract_bipartite_graph> &G,
                                                const shared_ptr<unordered_set<shared_ptr<abstract_bipartite_edge>>> &edge_set,
                                                const shared_ptr<unordered_map<shared_ptr<abstract_bipartite_edge>,uint32_t>>& edge_wing_map,
                                                const shared_ptr<unordered_map<uint32_t,shared_ptr<extend_list<uint32_t,shared_ptr<abstract_bipartite_edge>>>>>& wing_order_map,
                                                const shared_ptr<unordered_map<uint32_t, shared_ptr<unordered_map<shared_ptr<abstract_bipartite_edge>, uint32_t>>>> &WS,
                                                const shared_ptr<unordered_map<shared_ptr<abstract_bipartite_edge>,shared_ptr<mutex>>>& edge_mutex_map) {

        if(!WS->count(0))
        {
            WS->insert({0, make_shared<unordered_map<shared_ptr<abstract_bipartite_edge>,uint32_t>>()});
        }
        auto WS0 = WS->at(0);

        if(!WS->count(1))
        {
            WS->insert({1, make_shared<unordered_map<shared_ptr<abstract_bipartite_edge>,uint32_t>>()});
        }
        auto WS1 = WS->at(1);

        if (!wing_order_map->count(0)) {
            wing_order_map->insert(
                    {0, make_shared<extend_list<uint32_t, shared_ptr<abstract_bipartite_edge>>>()});
            WS->insert({0, make_shared<unordered_map<shared_ptr<abstract_bipartite_edge>, uint32_t>>()});
        }
        auto order_list = wing_order_map->at(1);

        for(const auto&e:*edge_set){
            WS1->insert({e, 0});
            G->insert_edge(e);
            edge_wing_map->insert({e, 0});
            order_list->left_insert(e);

            edge_mutex_map->insert({e, make_shared<mutex>()});
        }

        auto visited_edge_set = make_shared<unordered_set<shared_ptr<abstract_bipartite_edge>>>();
        for(const auto&e1:*edge_set){
            auto l1 = e1->get_left_vertex_id();
            auto r1 = e1->get_right_vertex_id();

            visited_edge_set->insert(e1);
            for (const auto &[r2, e2]:*G->get_left_vertex(l1)->get_edge_map()) {
                if (r2 == r1 || visited_edge_set->count(e2)) {
                    continue;
                }
                for (const auto &[l2, e3]:*G->get_right_vertex(r1)->get_edge_map()) {
                    if (l2 == l1 || visited_edge_set->count(e3)) {
                        continue;
                    }

                    auto e4 = G->get_edge(l2, r2);
                    if (!e4 || visited_edge_set->count(e4)) {
                        continue;
                    }

                    ++WS0->at(e2);
                    ++WS1->at(e2);

                    ++WS0->at(e2);
                    ++WS1->at(e2);

                    ++WS0->at(e3);
                    ++WS1->at(e3);

                    ++WS0->at(e4);
                    ++WS1->at(e4);
                }
            }

            if(WS1->at(e1) == 0){
                WS1->erase(e1);
            }else
            {
                edge_wing_map->at(e1) = 1;
            }
        }
    }

    void quasi_order_wing_maintenance::update_single_wing(const shared_ptr<abstract_bipartite_graph> &G,
                                                    const shared_ptr<unordered_set<shared_ptr<abstract_bipartite_edge>>> &edge_set,
                                                    const shared_ptr<unordered_map<shared_ptr<abstract_bipartite_edge>,uint32_t>>& edge_wing_map,
                                                    const shared_ptr<unordered_map<uint32_t,shared_ptr<extend_list<uint32_t,shared_ptr<abstract_bipartite_edge>>>>>& wing_order_map,
                                                    const shared_ptr<unordered_map<uint32_t, shared_ptr<unordered_map<shared_ptr<abstract_bipartite_edge>, uint32_t>>>> &WS,
                                                    uint32_t k) {

        auto WS_k = WS->at(k);
        auto super_order_list = wing_order_map->at(k - 1);
        auto current_order_list = wing_order_map->at(k);

        auto evicted_edge_set = make_shared<unordered_set<shared_ptr<abstract_bipartite_edge>>>();
        for(const auto&e1:*edge_set){
            WS_k->erase(e1);

            auto l1 = e1->get_left_vertex_id();
            auto r1 = e1->get_right_vertex_id();
            for (const auto &[r2, e2]:*G->get_left_vertex(l1)->get_edge_map()) {
                if (r2 == r1 || !WS_k->count(e2)) {
                    continue;
                }
                for (const auto &[l2, e3]:*G->get_right_vertex(r1)->get_edge_map()) {
                    if (l2 == l1 || !WS_k->count(e3)) {
                        continue;
                    }

                    auto e4 = G->get_edge(l2, r2);
                    if (!e4 || !WS_k->count(e4)) {
                        continue;
                    }

                    --WS_k->at(e2);
                    if (WS_k->at(e2) < k) {
                        evicted_edge_set->insert(e2);
                    }

                    --WS_k->at(e3);
                    if (WS_k->at(e3) < k) {
                        evicted_edge_set->insert(e3);
                    }

                    --WS_k->at(e4);
                    if (WS_k->at(e4) < k) {
                        evicted_edge_set->insert(e4);
                    }
                }
            }
        }

        while(!evicted_edge_set->empty()){
            auto e1 = *edge_set->begin();
            evicted_edge_set->erase(e1);
            WS_k->erase(e1);
            if(edge_wing_map->at(e1) >= k){
                edge_wing_map->at(e1) = k - 1;
            }
            current_order_list->remove(e1);

            if(!edge_set->count(e1)){
                auto e1_node = current_order_list->remove(e1);
                super_order_list->left_insert(e1_node);
            }

            auto l1 = e1->get_left_vertex_id();
            auto r1 = e1->get_right_vertex_id();
            for (const auto &[r2, e2]:*G->get_left_vertex(l1)->get_edge_map()) {
                if (r2 == r1 || !WS_k->count(e2)) {
                    continue;
                }
                for (const auto &[l2, e3]:*G->get_right_vertex(r1)->get_edge_map()) {
                    if (l2 == l1 || !WS_k->count(e3)) {
                        continue;
                    }

                    auto e4 = G->get_edge(l2, r2);
                    if (!e4 || !WS_k->count(e4)) {
                        continue;
                    }

                    --WS_k->at(e2);
                    if(WS_k->at(e2) < k){
                        evicted_edge_set->insert(e2);
                    }

                    --WS_k->at(e3);
                    if(WS_k->at(e3) < k){
                        evicted_edge_set->insert(e3);
                    }

                    --WS_k->at(e4);
                    if(WS_k->at(e4) < k){
                        evicted_edge_set->insert(e4);
                    }
                }
            }
        }
    }

    void quasi_order_wing_maintenance::update_single_wing(const shared_ptr<abstract_bipartite_graph> &G,
                                                    const shared_ptr<unordered_set<shared_ptr<abstract_bipartite_edge>>> &evicted_edge_set,
                                                    const shared_ptr<mutex> &evicted_edge_set_mutex,
                                                    const shared_ptr<unordered_map<shared_ptr<abstract_bipartite_edge>,uint32_t>>& edge_wing_map,
                                                    const shared_ptr<unordered_map<shared_ptr<abstract_bipartite_edge>,shared_ptr<mutex>>>& edge_mutex_map,
                                                    const shared_ptr<unordered_map<uint32_t,shared_ptr<extend_list<uint32_t,shared_ptr<abstract_bipartite_edge>>>>>& wing_order_map,
                                                    const shared_ptr<unordered_map<uint32_t, shared_ptr<unordered_map<shared_ptr<abstract_bipartite_edge>, uint32_t>>>> &WS,
                                                    uint32_t k,
                                                    uint32_t thread_count) {
        auto WS_k = WS->at(k);


        auto removal_edge_set = make_shared<unordered_set<shared_ptr<abstract_bipartite_edge>>>();
        auto task_vector = make_shared<vector<shared_ptr<unordered_set<shared_ptr<abstract_bipartite_edge>>>>>(thread_count,
                                                                                                               make_shared<unordered_set<shared_ptr<abstract_bipartite_edge>>>());
        thread_pool pool(thread_count);
        while(!evicted_edge_set->empty()){
            /**
             * @brief assign_value tasks
             */
            uint32_t count = 0;
            for(const auto &e:*evicted_edge_set){
                uint32_t index = count % thread_count;
                ++count;
                task_vector->at(index)->insert(e);
                removal_edge_set->insert(e);
            }

            evicted_edge_set->clear();
            for(auto i = 0;i< thread_count;++i){
                auto task = [=]{
                    auto sub_edge_set = task_vector->at(i);
                    auto sub_affected_edge_set = make_shared<unordered_set<shared_ptr<abstract_bipartite_edge>>>();
                    while(!sub_edge_set->empty()) {
                        auto e1 = *sub_edge_set->begin();
                        sub_edge_set->erase(e1);

                        edge_mutex_map->at(e1)->lock();
                        WS_k->at(e1) = 0;
                        edge_mutex_map->at(e1)->unlock();

                        auto l1 = e1->get_left_vertex_id();
                        auto r1 = e1->get_right_vertex_id();
                        for (const auto &[r2, e2]:*G->get_left_vertex(l1)->get_edge_map()) {
                            if (r2 == r1 || !WS_k->count(e2)) {
                                continue;
                            }
                            edge_mutex_map->at(e2)->lock();
                            if(WS_k->at(e2) < k)
                            {
                                edge_mutex_map->at(e2)->unlock();
                                continue;
                            } else
                            {
                                edge_mutex_map->at(e2)->unlock();
                            }
                            for (const auto &[l2, e3]:*G->get_right_vertex(r1)->get_edge_map()) {
                                if (l2 == l1) {
                                    continue;
                                }

                                auto e4 = G->get_edge(l2, r2);
                                if (!e4) {
                                    continue;
                                }

                                edge_mutex_map->at(e2)->lock();
                                if (WS_k->at(e2) > k) {
                                    --WS_k->at(e2);
                                    if (WS_k->at(e2) < k) {
                                        sub_affected_edge_set->insert(e2);
                                    }
                                }
                                edge_mutex_map->at(e2)->unlock();

                                edge_mutex_map->at(e3)->lock();
                                if (WS_k->at(e3) > k) {
                                    --WS_k->at(e3);
                                    if (WS_k->at(e3) < k) {
                                        sub_affected_edge_set->insert(e3);
                                    }
                                }
                                edge_mutex_map->at(e3)->unlock();

                                edge_mutex_map->at(e4)->lock();
                                if (WS_k->at(e4) > k) {
                                    --WS_k->at(e4);
                                    if (WS_k->at(e4) < k) {
                                        sub_affected_edge_set->insert(e4);
                                    }
                                }
                                edge_mutex_map->at(e4)->unlock();
                            }
                        }
                    }
                    task_vector->at(i)->clear();

                    evicted_edge_set_mutex->lock();
                    copy(sub_affected_edge_set->begin(), sub_affected_edge_set->end(), inserter(*evicted_edge_set, evicted_edge_set->end()));
                    evicted_edge_set_mutex->unlock();
                };
                pool.submit_task(task);
            }
            pool.barrier();
        }

        auto super_order_list = wing_order_map->at(k - 1);
        auto current_order_list = wing_order_map->at(k);
        for(const auto&e:*removal_edge_set){
            WS_k->erase(e);

            if(current_order_list->count_value(e)){
                current_order_list->remove(e);
            }

            if(edge_wing_map->at(e) >= k)
            {
                super_order_list->left_insert(e);
                edge_wing_map->at(e) = k-1;
            }
        }
    }
}