
#include "core/jes_core_maintenance.h"

namespace scnu {
    uint32_t jes_core_maintenance::compute_core_number(const shared_ptr<abstract_graph> &G,
                                                      const shared_ptr<unordered_map<uint32_t, uint32_t>> &core,
                                                      uint32_t v) {
        map<uint32_t,uint32_t> core_number_count_map;
        auto v_vertex = G->get_vertex(v);
        if(!v_vertex)
        {
            return 0;
        }
        for (const auto &[u,e]:*v_vertex->get_edge_map()) {
           for(auto k = 1;k<=core->at(u);++k)
           {
               if (!core_number_count_map.count(k)) {
                   core_number_count_map.insert({k, 0});
               }
               ++core_number_count_map.at(k);
           }
        }

        auto k_max = 1;
        for (auto iter = core_number_count_map.rbegin(); iter != core_number_count_map.rend(); ++iter) {
            auto core_number = iter->first;
            auto count = iter->second;
            if (count >= core_number) {
                k_max = core_number;
                break;
            }
        }
        return k_max;
    }


    /**
     * @details divide deleted edges into several edge sets
     * @param G the given graph
     * @param ED a set of deleted edges
     * @param core a map of vertices and core numbers
     * @param thread_number the number of threads
     * @return
     */
    pair<shared_ptr<map<uint32_t, shared_ptr<unordered_set<shared_ptr<abstract_edge>>>>>,shared_ptr<unordered_set<uint32_t>>>
    jes_core_maintenance::compute_delete_edge_set(const shared_ptr<abstract_graph> &G,
                                                      const shared_ptr<unordered_set<shared_ptr<abstract_edge>>> &ED,
                                                      const shared_ptr<unordered_map<uint32_t, uint32_t>> &core,
                                                      const shared_ptr<thread_pool>& pool) {
        auto EU = make_shared<map<uint32_t, shared_ptr<unordered_set<shared_ptr<abstract_edge>>>>>();
        auto superior_vertex_set = get_superior_vertex_set(ED, core);
        auto Vc = get_maximal_three_hop_independent_set(G, superior_vertex_set);

        G->remove_edge_collection(ED);
        auto pre_core_map = make_shared<unordered_map<uint32_t,uint32_t>>();
        auto pre_core_mutex = make_shared<mutex>();
        for(const auto&u:*Vc)
        {
            pool->submit_task([=] {
                auto pre_core_value = compute_pre_core_number(G, core, u);
                pre_core_mutex->lock();
                pre_core_map->insert({u, pre_core_value});
                pre_core_mutex->unlock();
            });
        }
        pool->barrier();
        G->insert_edge_collection(ED);

        auto mark_mutex = make_shared<mutex>() ;
        auto mark = make_shared<unordered_set<uint32_t>>();

        auto marked_edge_set = make_shared<unordered_set<shared_ptr<abstract_edge>>>();
        /**
         * @brief parallel part
         */
        for (const auto &u:*Vc) {
            pool->submit_task([=]() {
                auto init_core = core->at(u);
                /**
                 * @note: thread security
                 */
                core->at(u) = pre_core_map->at(u);
                auto u_vertex = G->get_vertex(u);
                for (const auto &[v, e]:*u_vertex->get_edge_map()) {
                    if (pre_core_map->at(u) <= core->at(v) && core->at(v) <= init_core) {
                        mark_mutex->lock();
                        if (!EU->count(core->at(v))) {
                            EU->insert({core->at(v), make_shared<unordered_set<shared_ptr<abstract_edge>>>()});
                        }
                        EU->at(core->at(v))->insert(e);
                        marked_edge_set->insert(e);
                        e->swap(v, u);
                        mark->insert(v);
                        mark_mutex->unlock();
                    }
                }
            });
        }
        pool->barrier();
        /**
         * @brief parallel part
         */
        for (const auto&e:*ED) {
            if(marked_edge_set->count(e)){
                continue;
            }

            pool->submit_task([=]() {
                auto x = e->get_source_vertex_id();
                auto y = e->get_destination_vertex_id();

                if (core->at(x) < core->at(y)) {
                    swap(x, y);
                }

                if (Vc->count(y) || core->at(y) < core->at(x)) {
                    mark_mutex->lock();
                    if (!mark->count(y)) {
                        if (!EU->count(core->at(y))) {
                            EU->insert({core->at(y), make_shared<unordered_set<shared_ptr<abstract_edge>>>()});
                        }
                        EU->at(core->at(y))->insert(e);
                        e->swap(y, x);
                        mark->insert(y);
                    }
                    mark_mutex->unlock();
                }

                if (core->at(x) == core->at(y)) {
                    mark_mutex->lock();
                    if (!mark->count(x) && !mark->count(y)) {
                        if (!EU->count(core->at(x))) {
                            EU->insert({core->at(x), make_shared<unordered_set<shared_ptr<abstract_edge>>>()});
                        }
                        EU->at(core->at(x))->insert(e);
                        mark->insert(x);
                        mark->insert(y);
                    }
                    mark_mutex->unlock();
                }
            });
        }
        pool->barrier();
        return make_pair(EU, Vc);
    }

    /**
     * @details divide inserted edges into several edge sets
     * @param G the given graph
     * @param EI a set of inserted edges
     * @param core a map of vertices and core numbers
     * @param thread_number the number of running threads
     * @return 
     */
    pair<shared_ptr<unordered_map<uint32_t, shared_ptr<unordered_set<shared_ptr<abstract_edge>>>>>,shared_ptr<unordered_set<uint32_t>>>
    jes_core_maintenance::compute_insert_edge_set(const shared_ptr<abstract_graph> &G,
                                                      const shared_ptr<unordered_set<shared_ptr<abstract_edge>>> &EI,
                                                      const shared_ptr<unordered_map<uint32_t, uint32_t>> &core,
                                                      const shared_ptr<thread_pool>& pool) {
        auto EU = make_shared<unordered_map<uint32_t, shared_ptr<unordered_set<shared_ptr<abstract_edge>>>>>();
        auto EU_mutex = make_shared<mutex>();
        /**
          * @brief get G+EI
          * @remarks  next loop will operate on current graph
        */
        auto superior_vertex_set = get_superior_vertex_set(EI, core);
        auto marked_edge_set = make_shared<unordered_set<shared_ptr<abstract_edge>>>();
        for(const auto&e:*EI){
            G->insert_edge(e);
        }
        auto Vc = get_maximal_three_hop_independent_set(G, superior_vertex_set);
        /**
         * @brief parallel classify associated edges in Vc
         * @remarks revise core number of vertices in Vc do not affect its neighbors
         */
        auto mark_mutex = make_shared<mutex>() ;
        auto mark = make_shared<unordered_set<uint32_t>>();
        for (const auto &u:*Vc) {
            pool->submit_task([=] {
                auto init_core_number = core->at(u);
                auto pre_core_value = compute_pre_core_number(G, core, u);
                core->at(u) = pre_core_value;
                auto u_vertex = G->get_vertex(u);
                for (const auto &[v, e]:*u_vertex->get_edge_map()) {
                    if (init_core_number <= core->at(v) && core->at(v) <= pre_core_value) {
                        mark_mutex->lock();
                        if (!EU->count(core->at(v))) {
                            EU->insert({core->at(v), make_shared<unordered_set<shared_ptr<abstract_edge>>>()});
                        }
                        EU->at(core->at(v))->insert(e);
                        e->swap(v, u);
                        marked_edge_set->insert(e);
                        mark->insert(v);
                        mark_mutex->unlock();
                    }
                }
            });
        }
        pool->barrier();


        /**
         * @brief parallel classify associated edges in EI
         */
        for (const auto &e:*EI) {
            if(marked_edge_set->count(e)){
                continue;
            }
            pool->submit_task([=] {
                auto x = e->get_source_vertex_id();
                auto y = e->get_destination_vertex_id();
                /**
                 * @remarks enable core->at(y)<=core->at(x)
                 * @details x is root vertex
                 */
                if (core->at(x) < core->at(y)) {
                    swap(x, y);
                }

                if (Vc->count(y) || core->at(y) < core->at(x)) {
                    mark_mutex->lock();
                    if (!mark->count(y)) {
                        if (!EU->count(core->at(y))) {
                            EU->insert({core->at(y), make_shared<unordered_set<shared_ptr<abstract_edge>>>()});
                        }
                        EU->at(core->at(y))->insert(e);
                        e->swap(y, x);
                        mark->insert(y);
                    }
                    mark_mutex->unlock();
                }

                if (core->at(x) == core->at(y)) {
                    mark_mutex->lock();
                    if (!mark->count(x) && !mark->count(y)) {
                        if (!EU->count(core->at(x))) {
                            EU->insert({core->at(x), make_shared<unordered_set<shared_ptr<abstract_edge>>>()});
                        }
                        EU->at(core->at(x))->insert(e);
                        mark->insert(x);
                        mark->insert(y);
                    }
                    mark_mutex->unlock();
                }
            });
        }
        pool->barrier();


        /**
          * @brief get G-EI
          * @remarks  next loop will operate on previous graph
        */
        G->remove_edge_collection(EI);

        return make_pair(EU, Vc);
    }

    /**
     * @details compute pre-core number of a given vertex
     * @remarks the function is similar to compute core number since they share similar idea
     * @param G the given graph
     * @param core a map of vertices and core numbers
     * @param v the given vertex
     * @return
     */
    uint32_t
    jes_core_maintenance::compute_pre_core_number(const shared_ptr<abstract_graph> &G,
                                                      const shared_ptr<unordered_map<uint32_t, uint32_t>> &core,
                                                      uint32_t v) {
        return compute_core_number(G,core,v);
    }

    /**
     * @details get maximal 3-hop independent set of a given set
     * @param G the given graph
     * @param superior_vertex_set  a set of superior vertices
     * @return
     */


    shared_ptr<unordered_set<uint32_t>> jes_core_maintenance::get_maximal_three_hop_independent_set(
            const shared_ptr<abstract_graph> &G,
            const shared_ptr<unordered_set<uint32_t>> &superior_vertex_set) {
        auto Vc = make_shared<unordered_set<uint32_t>>();
        /**
         * @brief when the superior vertex set is empty, directly return
         */
        if(superior_vertex_set->empty())
        {
            return Vc;
        }

        while(!superior_vertex_set->empty()){
            auto p = superior_vertex_set->begin();
            auto u = *p;

            superior_vertex_set->erase(u);
            Vc->insert(u);
            auto u_vertex = G->get_vertex(u);
            for (const auto&[v, e]:*u_vertex->get_edge_map()) {
                superior_vertex_set->erase(v);
                auto v_vertex = G->get_vertex(v);
                for (const auto&[w, e]:*v_vertex->get_edge_map()) {
                    superior_vertex_set->erase(w);
                }
            }
        }
//        while(p!=superior_vertex_set->end())
//        {
//            auto u = *p;
//            auto u_vertex = G->get_vertex(u);
//            bool flag = false;
//            for (const auto &v:*Vc) {
//                if (u_vertex->get_edge(v)) {
//                    flag = true;
//                    break;
//                }
//            }
//            if(!flag){
//                for(const auto&v:*Vc){
//                    auto v_vertex = G->get_vertex(v);
//                    for (auto &[w, e]:*v_vertex->get_edge_map()) {
//                        if (u_vertex->get_edge(w)) {
//                            flag = true;
//                            break;
//                        }
//                    }
//                    if(flag){
//                        break;
//                    }
//                }
//            }
//            if(!flag)
//            {
//                Vc->insert(u);
//            }
//            ++p;
//        }
        return Vc;
    }

    shared_ptr<unordered_set<uint32_t>>
    jes_core_maintenance::get_k_path_tree(const shared_ptr<abstract_graph> &G,
                                              const shared_ptr<unordered_map<uint32_t, uint32_t>> &core,
                                              uint32_t w,
                                              uint32_t k
    ) {
        /**
         * @brief initialise
         */
        queue<uint32_t> vertex_queue;
        unordered_set<uint32_t> visited_set;
        vertex_queue.push(w);
        visited_set.insert(w);
        auto result_set = make_shared<unordered_set<uint32_t>>();
        while (!vertex_queue.empty()) {
            auto u = vertex_queue.front();
            vertex_queue.pop();
            result_set->insert(u);
            auto u_vertex = G->get_vertex(u);
            if(!u_vertex)
            {
                continue;
            }
            for (const auto &[v,e]:*u_vertex->get_edge_map()) {
                if (core->at(v) == k && !visited_set.count(v)) {
                    vertex_queue.push(v);
                    visited_set.insert(v);
                }
            }
        }
        return result_set;
    }

    shared_ptr<unordered_set<uint32_t>>
    jes_core_maintenance::get_k_path_tree(const shared_ptr<abstract_graph> &G,
                                              const shared_ptr<unordered_map<uint32_t, uint32_t>> &core,
                                              const shared_ptr<unordered_set<uint32_t>>& vertex_set,
                                              uint32_t k
    ) {
        /**
         * @brief initialise
         */
        queue<uint32_t> vertex_queue;
        unordered_set<uint32_t> visited_set;
        for (const auto &w:*vertex_set) {
            vertex_queue.push(w);
            visited_set.insert(w);
        }
        auto result_set = make_shared<unordered_set<uint32_t>>();
        while (!vertex_queue.empty()) {
            auto u = vertex_queue.front();
            vertex_queue.pop();
            result_set->insert(u);
            auto u_vertex = G->get_vertex(u);
            if (!u_vertex) {
                continue;
            }
            for (const auto &[v, e]:*u_vertex->get_edge_map()) {
                if (core->at(v) == k && !visited_set.count(v)) {
                    vertex_queue.push(v);
                    visited_set.insert(v);
                }
            }
        }
        return result_set;
    }

    /**
     * @details get the superior degree of a given vertex, it is equal to the number of superior edges
     * @param G
     * @param core
     * @param u
     * @return superior degree
     */
    uint32_t jes_core_maintenance::get_superior_degree(const shared_ptr<abstract_graph> &G,
                                                           const shared_ptr<unordered_map<uint32_t, uint32_t>> &core,
                                                           uint32_t u) {
        uint32_t superior_degree = 0;
        auto u_vertex = G->get_vertex(u);
        if(!u_vertex)
        {
            return 0;
        }
        for (const auto &[v,e]:*u_vertex->get_edge_map()) {
            if (core->at(u) <= core->at(v)) {
                ++superior_degree;
            }
        }
        return superior_degree;
    }


    /**
     * @details get a set of vertices having at least two superior edges
     * @param E a set of inserted (removed) edges
     * @param core a map of vertices and core numbers
     * @return
     */
    shared_ptr<unordered_set<uint32_t>>
    jes_core_maintenance::get_superior_vertex_set(const shared_ptr<unordered_set<shared_ptr<abstract_edge>>> &E,
                                                      const shared_ptr<unordered_map<uint32_t, uint32_t>> &core) {
        auto superior_vertex_set = make_shared<unordered_set<uint32_t>>();
        auto result_map = make_shared<unordered_map<uint32_t,uint32_t>>();
        for (const auto &e:*E) {
            auto u = e->get_source_vertex_id();
            auto v = e->get_destination_vertex_id();
            if(!core->count(u))
            {
                core->insert({u,1});
            }
            if(!core->count(v))
            {
                core->insert({v,1});
            }
            if (core->at(u)<=core->at(v)) {
                if(!result_map->count(u)){
                    result_map->insert({u,0});
                }
                ++result_map->at(u);
            }
            if (core->at(v) <= core->at(u))
            {
                if(!result_map->count(v)){
                    result_map->insert({v,0});
                }
                ++result_map->at(v);
            }
        }
        for(const auto& [u,superior_edge_number]:*result_map)
        {
            if(superior_edge_number >= 2)
            {
                superior_vertex_set->insert(u);
            }
        }
        return superior_vertex_set;
    }

    /**
     * @details a generalized traversal method, it can deal with a set of deleted edges simultaneously
     * @param G the given graph
     * @param ED a set of deleted edges
     * @param core a map of vertices and their core numbers
     * @param thread_number the number of threads
     */
    void jes_core_maintenance:: joint_delete(const shared_ptr<abstract_graph> &G,
                                                const shared_ptr<unordered_set<shared_ptr<abstract_edge>>> &ED,
                                                const shared_ptr<unordered_map<uint32_t, uint32_t>> &core,
                                                uint32_t thread_number) {
        auto pool = make_shared<thread_pool>(thread_number);
        auto isolated_vertex_set = make_shared<unordered_set<uint32_t>>();
        while (!ED->empty()) {
            auto [EU, Vc] = compute_delete_edge_set(G, ED, core, pool);

            for (const auto &[k,EUk]:*EU) {
                for (const auto &e : *EU->at(k)) {
                    if(!ED->count(e)){
                        continue;
                    }
                    auto sub_isolated_vertex_set = G->remove_edge(e);
                    for (const auto &u:*sub_isolated_vertex_set) {
                        isolated_vertex_set->insert(u);
                    }
                    ED->erase(e);
                }
            }

            /**
             * @brief parallel part
             */
            auto V = make_shared<safe_unordered_map<uint32_t,shared_ptr<unordered_set<uint32_t>>>>();
            for (const auto &k_pair:*EU) {
                pool->submit_task([=]() {
                    auto k = k_pair.first;
                    auto Ek = k_pair.second;
                    auto Vk = k_joint_delete(G, Ek, core, k);
                    V->insert({k, Vk});
                });
            }
            pool->barrier();

            for(const auto&[k,Vk]:*V)
            {
                for (const auto &v:*Vk) {
                    --core->at(v);
                }
            }

            for (const auto &v:*Vc) {
                core->at(v) = compute_core_number(G, core, v);
            }
        }
        for(const auto&v:*isolated_vertex_set)
        {
            core->erase(v);
        }
    }


    /**
     * @details a generalized traversal maintenance method, it can deal with a set of inserted edges simultaneously
     * @param G the given graph
     * @param EI a set of inserted edges
     * @param core a map of vertices and core numbers
     * @param thread_number the number of threads
     */
    void jes_core_maintenance::joint_insert(const shared_ptr<abstract_graph> &G,
                                                const shared_ptr<unordered_set<shared_ptr<abstract_edge>>> &EI,
                                                const shared_ptr<unordered_map<uint32_t, uint32_t>> &core,
                                                uint32_t thread_number) {
        auto pool = make_shared<thread_pool>(thread_number);
        while (!EI->empty()) {
            auto [EU, Vc] = compute_insert_edge_set(G, EI, core, pool);
            for (const auto &[k,Ek]:*EU) {
                for (const auto &e : *Ek) {
                    G->insert_edge(e);
                    EI->erase(e);
                }
            }

            /**
              * @brief parallel part
              */
            auto V = make_shared<safe_unordered_map<uint32_t, shared_ptr<unordered_set<uint32_t>>>>();
            for (const auto &k_pair:*EU) {
                auto k = k_pair.first;
                auto Ek = k_pair.second;
                pool->submit_task([=] {
                    auto Vk = k_joint_insert(G, Ek, core, k);
                    V->insert({k, Vk});
                });
            }
            pool->barrier();

            for(const auto& k_pair:*V)
            {
                auto Vk = k_pair.second;
                for(const auto& v:*Vk)
                {
                    core->at(v) = core->at(v) + 1;
                }
            }

            for (const auto &v:*Vc) {
                auto core_number = compute_core_number(G, core, v);
                core->at(v) = core_number;
            }
        }
    }

    /**
     * @details find a influenced vertex set by using k-joint deleted edge set
     * @param G the given graph
     * @param Ek a k-joint set of deleted edges
     * @param core a map of vertices and their core numbers
     * @param k the given core number
     * @return a set of affected vertices
     */
    shared_ptr<unordered_set<uint32_t>>
    jes_core_maintenance::k_joint_delete(const shared_ptr<abstract_graph> &G,
                                             const shared_ptr<unordered_set<shared_ptr<abstract_edge>>> &Ek,
                                             const shared_ptr<unordered_map<uint32_t, uint32_t>> &core,
                                             uint32_t k) {
        auto cd = make_shared<unordered_map<uint32_t, long long>>();
        auto removed_edge_set = make_shared<unordered_set<shared_ptr<abstract_edge>>>();


        auto vertex_set = make_shared<unordered_set<uint32_t>>();
        for (const auto &e : *Ek) {
            auto u = e->get_source_vertex_id();
            auto v = e->get_destination_vertex_id();
            vertex_set->insert(u);
            if(core->at(u) == core->at(v)){
                vertex_set->insert(v);
            }
        }
        auto KPT = get_k_path_tree(G,core,vertex_set,k);
        for(const auto&w:*KPT)
        {
            cd->insert({w,get_superior_degree(G,core,w)});
        }


       auto evicted = make_shared<unordered_set<uint32_t>>();
        for (const auto &e : *Ek) {
            auto u = e->get_source_vertex_id();
            auto v = e->get_destination_vertex_id();
            unordered_set<uint32_t> C;
            if (!evicted->count(u) && cd->at(u) < k) {
                C.insert(u);
            }

            if (core->at(u) == core->at(v) && !evicted->count(v) && cd->at(v) < k) {
                C.insert(v);
            }
            for(const auto w:C) {
                stack<uint32_t> Stk;
                Stk.push(w);
                evicted->insert(w);
                while (!Stk.empty()) {
                    auto x = Stk.top();
                    Stk.pop();
                    auto x_vertex = G->get_vertex(x);
                    if (!x_vertex) {
                        continue;
                    }
                    for (const auto&[y, e]:*x_vertex->get_edge_map()) {
                        if (core->at(y) == k && cd->at(y) >= k) {
                            --cd->at(y);
                            if (cd->at(y) < k && !evicted->count(y)) {
                                Stk.push(y);
                                evicted->insert(y);
                            }
                        }
                    }
                }
            }
        }
        return evicted;
    }

    /**
     * @details find a influenced vertex set by using k-joint inserted edge set
     * @param G the given graph
     * @param Ek a k-joint set of inserted edges
     * @param core a map of vertices and their core numbers
     * @param k the given core number
     * @return a set of affected vertices
     */
    shared_ptr<unordered_set<uint32_t>>
    jes_core_maintenance::k_joint_insert(const shared_ptr<abstract_graph> &G,
                                             const shared_ptr<unordered_set<shared_ptr<abstract_edge>>> &Ek,
                                             const shared_ptr<unordered_map<uint32_t, uint32_t>> &core,
                                             uint32_t k) {
        /**
         * @brief initialize variables
         */
        stack<uint32_t> Stk;
        auto Vt = make_shared<unordered_set<uint32_t>>();
        auto evicted = make_shared<unordered_set<uint32_t>>();
        unordered_set<uint32_t> visited;
        /**
         * @remark the value in cd may be less than 0
         */
        auto cd = make_shared<unordered_map<uint32_t, long long>>();
//        for(const auto&[u,core_number]:*core){
//            if(core_number == k){
//                cd->insert({u,get_superior_degree(G, core, u)});
//            }
//        }

        auto vertex_set = make_shared<unordered_set<uint32_t>>();
        for(const auto&e:*Ek)
        {
            auto u = e->get_source_vertex_id();
            //auto v = e->get_destination_vertex_id();
            vertex_set->insert(u);
        }
        auto KPT = get_k_path_tree(G,core,vertex_set,k);
        for(const auto&w:*KPT){
            cd->insert({w,get_superior_degree(G,core,w)});
        }

        for (const auto &e : *Ek) {
            auto u = e->get_source_vertex_id();
            //auto v = e->get_destination_vertex_id();
            if(!visited.count(u))
            {
                if(cd->at(u) > k){
                    Stk.push(u);
                    visited.insert(u);
                }
            }
        }

        while (!Stk.empty()) {
            auto u = Stk.top();
            Stk.pop();
            if (cd->at(u) > k) {
                auto u_vertex = G->get_vertex(u);
                for (const auto &[v,e]:*u_vertex->get_edge_map()) {
                     if (!visited.count(v) && core->at(v) == k
                     //&& cd->at(v)>k
                     ) {
                         Stk.push(v);
                         visited.insert(v);
                     }
                }
            } else if (!evicted->count(u)) {
                stack<uint32_t> S2;
                S2.push(u);
                evicted->insert(u);
                while (!S2.empty()) {
                    auto x = S2.top();
                    S2.pop();
                    auto x_vertex = G->get_vertex(x);
                    for (const auto &[y,e]:*x_vertex->get_edge_map()) {
                        if (core->at(y)==k) {
                            --cd->at(y);
                            if (!evicted->count(y) && cd->at(y) <= k) {
                                S2.push(y);
                                evicted->insert(y);
                            }
                        }
                    }
                }
            }
        }

        for (const auto &v:visited) {
            if (!evicted->count(v)) {
                Vt->insert(v);
            }
        }
        return Vt;
    }
}