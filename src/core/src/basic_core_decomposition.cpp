
#include "core/basic_core_decomposition.h"

namespace scnu {
    /**
     * @details base core decomposition to get the core number of each vertex
     * @remarks the vertex id may not consecutive
     * @param G
     * @param vertex_core_map
     * @return
     */
    uint32_t basic_core_decomposition::decompose(const shared_ptr<abstract_graph> &G,
                                                 const shared_ptr<unordered_map<uint32_t, uint32_t>> &vertex_core_map) {
        auto vertex_set = make_shared<unordered_set<uint32_t>>();
        auto vertex_degree_map = make_shared<unordered_map<uint32_t,uint32_t>>();
        for(auto [v,v_vertex]:*G->get_vertex_map()){
            vertex_set->insert(v);
            vertex_degree_map->insert({v, v_vertex->get_degree()});
        }

        uint32_t k_max = 1;
        uint32_t  k = 1;
        auto evicted_set = make_shared<unordered_set<uint32_t>>();
        while(!vertex_set->empty()){
            k_max = k;

            for(const auto&u:*vertex_set){
                if(vertex_degree_map->at(u) <= k){
                    evicted_set->insert(u);
                }
            }

            while(!evicted_set->empty()){
                auto u = *evicted_set->begin();
                evicted_set->erase(u);
                vertex_set->erase(u);
                vertex_core_map->insert({u,k});

                for(const auto&[v,e]:*G->get_vertex(u)->get_edge_map()){
                    if(!vertex_set->count(v) || evicted_set->count(v)){
                        continue;
                    }

                    --vertex_degree_map->at(v);
                    if (vertex_degree_map->at(v) <= k) {
                        evicted_set->insert(v);
                    }
                }
            }
            ++k;
        }
        return k_max;
    }

    uint32_t basic_core_decomposition::decompose(const shared_ptr<abstract_graph> &G,
                                                 const shared_ptr<unordered_map<uint32_t, uint32_t>> &vertex_core_map,
                                                 uint32_t thread_number) {
        uint32_t k_max = 1;
        uint32_t  k = 1;

        auto vertex_set = make_shared<unordered_set<uint32_t>>();
        auto vertex_degree_map = make_shared<unordered_map<uint32_t, uint32_t>>();
        auto vertex_mutex_map = make_shared<unordered_map<uint32_t, shared_ptr<mutex>>>(G->get_vertex_number()+1);

        for(auto [v,v_vertex]:*G->get_vertex_map()){
            vertex_set->insert(v);
            vertex_degree_map->insert({v, v_vertex->get_degree()});
            vertex_mutex_map->insert({v, make_shared<mutex>()});
        }

        auto current_task_vector = make_shared<vector<shared_ptr<unordered_set<uint32_t>>>>(thread_number);
        for(uint32_t i = 0;i<thread_number;++i){
            current_task_vector->at(i) = make_shared<unordered_set<uint32_t>>();
        }

        auto pool = make_shared<thread_pool>(thread_number);
        auto evicted_set = make_shared<unordered_set<uint32_t>>();
        while(!vertex_set->empty()){
            k_max = k;

            for(const auto&u:*vertex_set){
                if(vertex_degree_map->at(u) <= k){
                    evicted_set->insert(u);
                }
            }

            while(!evicted_set->empty()){
                for(const auto&u:*evicted_set){
                    vertex_core_map->insert({u,k});

                    pool->submit_task([=]{
                        auto thread_id = pool->get_thread_id(std::this_thread::get_id());
                        for(const auto&[v,e]:*G->get_vertex(u)->get_edge_map()){
                            if(!vertex_set->count(v) || evicted_set->count(v)){
                                continue;
                            }

                            vertex_mutex_map->at(v)->lock();
                            --vertex_degree_map->at(v);
                            if(vertex_degree_map->at(v) <= k){
                                current_task_vector->at(thread_id)->insert(v);
                            }
                            vertex_mutex_map->at(v)->unlock();
                        }
                    });
                }
                pool->barrier();
                for(const auto u:*evicted_set){
                    vertex_set->erase(u);
                }
                evicted_set->clear();
                merger_set(current_task_vector, evicted_set, pool);
            }
            ++k;
        }
        return k_max;
    }

    /**
    * @details basic core decomposition to get core number and core degree of each vertex
    * @remarks the vertex id may not consecutive
    * @param G
    * @param vertex_core_map
    * @param CD
    * @return
    */
    uint32_t basic_core_decomposition::decompose(const shared_ptr<abstract_graph> &G,
                                                 const shared_ptr<unordered_map<uint32_t,uint32_t>>& vertex_core_map,
                                                 const shared_ptr<unordered_map<uint32_t,shared_ptr<unordered_map<uint32_t,uint32_t>>>>& CD) {
        auto vertex_set = make_shared<unordered_set<uint32_t>>();
        auto vertex_degree_map = make_shared<unordered_map<uint32_t,uint32_t>>();
        for(auto [v,v_vertex]:*G->get_vertex_map()){
            vertex_set->insert(v);
            vertex_degree_map->insert({v, v_vertex->get_degree()});
        }

        uint32_t k_max = 1;
        uint32_t  k = 1;
        auto evicted_set = make_shared<unordered_set<uint32_t>>();
        while(!vertex_set->empty()){
            k_max = k;

            CD->insert({k, make_shared<unordered_map<uint32_t,uint32_t>>()});
            for(const auto&u:*vertex_set){
                CD->at(k)->insert({u, vertex_degree_map->at(u)});
                if(vertex_degree_map->at(u) <= k){
                    evicted_set->insert(u);
                }
            }

            while(!evicted_set->empty()){
                auto u = *evicted_set->begin();
                evicted_set->erase(u);

                vertex_set->erase(u);

                vertex_core_map->insert({u,k});

                for(const auto&[v,e]:*G->get_vertex(u)->get_edge_map()){
                    if(!vertex_set->count(v) || evicted_set->count(v)){
                        continue;
                    }

                    --vertex_degree_map->at(v);
                    if (vertex_degree_map->at(v) <= k) {
                        evicted_set->insert(v);
                    }
                }
            }
            ++k;
        }
        return k_max;
    }

    /**
     * @details basic core decomposition to get a core order list
     * @param G
     * @param vertex_core_map
     * @param core_order
     * @return
     */
    uint32_t basic_core_decomposition::decompose(const shared_ptr<abstract_graph> &G,
                                                 const shared_ptr<unordered_map<uint32_t,uint32_t>> &vertex_core_map,
                                                 const shared_ptr<unordered_map<uint32_t,shared_ptr<extend_list<double,uint32_t>>>> &core_order) {
//        uint32_t k_max = 1;
//        uint32_t  k = 1;
//
//        auto vertex_degree_map = make_shared<unordered_map<uint32_t,uint32_t>>();
//        for(auto [v,v_vertex]:*G->get_vertex_map()){
//            vertex_degree_map->insert({v, v_vertex->get_neighbor_size()});
//        }
//
//        auto vertex_set = make_shared<unordered_set<uint32_t>>();
//        while(!vertex_degree_map->empty()){
//            k_max = k;
//            /**
//             * @brief prepare core order list
//             */
//            if (!core_order->count(k)) {
//                core_order->insert({k, make_shared<extend_list<double, uint32_t>>()});
//            }
//
//            for(const auto&[u,u_degree]:*vertex_degree_map){
//                if(u_degree <= k){
//                    vertex_set->insert(u);
//                }
//            }
//            while(!vertex_set->empty()){
//                auto u = *vertex_set->begin();
//                vertex_set->erase(u);
//
//                vertex_core_map->insert({u,k});
//                core_order->at(k)->push_back(u);
//
//                for(const auto&[v,e]:*G->get_vertex(u)->get_edge_map()){
//                    if(!vertex_degree_map->count(v) || vertex_degree_map->at(v) <= k){
//                        continue;
//                    }
//                    --vertex_degree_map->at(v);
//                    if(vertex_degree_map->at(v) == k){
//                        vertex_set->insert(v);
//                    }
//                }
//                vertex_degree_map->erase(u);
//            }
//            ++k;
//        }
//        return k_max;
        auto vertex_set = make_shared<unordered_set<uint32_t>>();
        auto vertex_degree_map = make_shared<unordered_map<uint32_t,uint32_t>>();
        for(auto [v,v_vertex]:*G->get_vertex_map()){
            vertex_set->insert(v);
            vertex_degree_map->insert({v, v_vertex->get_degree()});
        }

        uint32_t k_max = 1;
        uint32_t  k = 1;
        auto evicted_set = make_shared<unordered_set<uint32_t>>();
        while(!vertex_set->empty()){
            k_max = k;

            if (!core_order->count(k)) {
                core_order->insert({k, make_shared<extend_list<double, uint32_t>>()});
            }
            for(const auto&u:*vertex_set){
                if(vertex_degree_map->at(u) <= k){
                    evicted_set->insert(u);
                }
            }

            while(!evicted_set->empty()){
                auto u = *evicted_set->begin();
                evicted_set->erase(u);
                vertex_set->erase(u);
                vertex_core_map->insert({u,k});

                core_order->at(k)->push_back(u);

                for(const auto&[v,e]:*G->get_vertex(u)->get_edge_map()){
                    if(!vertex_set->count(v) || evicted_set->count(v)){
                        continue;
                    }

                    --vertex_degree_map->at(v);
                    if (vertex_degree_map->at(v) <= k) {
                        evicted_set->insert(v);
                    }
                }
            }
            ++k;
        }
        return k_max;
    }


    uint32_t basic_core_decomposition::decompose(const shared_ptr<abstract_graph> &G,
                                                 const shared_ptr<unordered_map<uint32_t,uint32_t>> &vertex_core_map,
                                                 const shared_ptr<unordered_map<uint32_t,shared_ptr<extend_list<double,uint32_t>>>> &core_order,
                                                 uint32_t thread_number) {
        uint32_t k_max = 1;
        uint32_t  k = 1;

        auto vertex_set = make_shared<unordered_set<uint32_t>>();
        auto vertex_degree_map = make_shared<unordered_map<uint32_t, uint32_t>>();
        auto vertex_mutex_map = make_shared<unordered_map<uint32_t, shared_ptr<mutex>>>();

        for(auto [v,v_vertex]:*G->get_vertex_map()){
            vertex_set->insert(v);
            vertex_degree_map->insert({v, v_vertex->get_degree()});
            vertex_mutex_map->insert({v, make_shared<mutex>()});
        }

        auto current_task_vector = make_shared<vector<shared_ptr<unordered_set<uint32_t>>>>(thread_number);
        for(uint32_t i = 0;i<thread_number;++i){
            current_task_vector->at(i) = make_shared<unordered_set<uint32_t>>();
        }

        auto pool = make_shared<thread_pool>(thread_number);
        auto evicted_set = make_shared<unordered_set<uint32_t>>();
        while(!vertex_set->empty()){
            k_max = k;

            if (!core_order->count(k)) {
                core_order->insert({k, make_shared<extend_list<double, uint32_t>>()});
            }
            for(const auto&u:*vertex_set){
                if(vertex_degree_map->at(u) <= k){
                    evicted_set->insert(u);
                }
            }

            while(!evicted_set->empty()){
                for(const auto&u:*evicted_set){
                    vertex_core_map->insert({u,k});
                    core_order->at(k)->push_back(u);

                    pool->submit_task([=]{
                        auto thread_id = pool->get_thread_id(std::this_thread::get_id());
                        for(const auto&[v,e]:*G->get_vertex(u)->get_edge_map()){
                            if(!vertex_set->count(v) || evicted_set->count(v)){
                                continue;
                            }

                            vertex_mutex_map->at(v)->lock();
                            --vertex_degree_map->at(v);
                            if(vertex_degree_map->at(v) <= k){
                                current_task_vector->at(thread_id)->insert(v);
                            }
                            vertex_mutex_map->at(v)->unlock();
                        }
                    });
                }
                pool->barrier();
                for(const auto u:*evicted_set){
                    vertex_set->erase(u);
                }
                evicted_set->clear();
                merger_set(current_task_vector, evicted_set, pool);
            }
            ++k;
        }
        return k_max;
    }

    /**
     * @details produce a double list and an order
     * @param G
     * @param vertex_core_map
     * @param k_order
     * @param node_map
     * @param tree
     * @return
     */
    uint32_t basic_core_decomposition::decompose(const shared_ptr<abstract_graph> &G,
                                                           const shared_ptr<unordered_map<uint32_t,uint32_t>> &vertex_core_map,
                                                           const shared_ptr<unordered_map<uint32_t,shared_ptr<double_list<uint32_t>>>>& k_order,
                                                           const shared_ptr<unordered_map<uint32_t, shared_ptr<double_node<uint32_t>>>> &node_map,
                                                           const shared_ptr<unordered_map<uint32_t,shared_ptr<extend_list<double,uint32_t>>>> &tree) {
        auto vertex_set = make_shared<unordered_set<uint32_t>>();
        auto vertex_degree_map = make_shared<unordered_map<uint32_t, uint32_t>>();
        for(auto [v,v_vertex]:*G->get_vertex_map()){
            vertex_set->insert(v);
            vertex_degree_map->insert({v, v_vertex->get_degree()});
        }

        uint32_t k_max = 1;
        uint32_t  k = 1;
        auto evicted_set = make_shared<unordered_set<uint32_t>>();
        while(!vertex_set->empty()){
            k_max = k;

            if (!tree->count(k)) {
                k_order->insert({k, make_shared<double_list<uint32_t>>()});
                tree->insert({k, make_shared<extend_list<double, uint32_t>>()});
            }
            for(const auto&u:*vertex_set){
                if(vertex_degree_map->at(u) <= k){
                    evicted_set->insert(u);
                }
            }

            while(!evicted_set->empty()){
                auto u = *evicted_set->begin();
                evicted_set->erase(u);

                vertex_set->erase(u);

                vertex_core_map->insert({u,k});

                auto u_node = make_shared<double_node<uint32_t>>(u);
                node_map->insert({u,u_node});
                k_order->at(k)->right_insert(u_node);

                tree->at(k)->push_back(u);

                for(const auto&[v,e]:*G->get_vertex(u)->get_edge_map()){
                    if(!vertex_set->count(v) || evicted_set->count(v)){
                        continue;
                    }

                    --vertex_degree_map->at(v);
                    if (vertex_degree_map->at(v) <= k) {
                        evicted_set->insert(v);
                    }
                }
            }
            ++k;
        }
        return k_max;
    }

    uint32_t basic_core_decomposition::decompose(const shared_ptr<abstract_graph> &G,
                                                 const shared_ptr<unordered_map<uint32_t,uint32_t>> &vertex_core_map,
                                                 const shared_ptr<unordered_map<uint32_t,shared_ptr<double_list<uint32_t>>>>& k_order,
                                                 const shared_ptr<unordered_map<uint32_t, shared_ptr<double_node<uint32_t>>>> &node_map,
                                                 const shared_ptr<unordered_map<uint32_t,shared_ptr<extend_list<double,uint32_t>>>> &tree,
                                                 uint32_t thread_number){
        uint32_t k_max = 1;
        uint32_t  k = 1;

        auto vertex_set = make_shared<unordered_set<uint32_t>>();
        auto vertex_degree_map = make_shared<unordered_map<uint32_t, uint32_t>>();
        auto vertex_mutex_map = make_shared<unordered_map<uint32_t, shared_ptr<mutex>>>();

        for(auto [v,v_vertex]:*G->get_vertex_map()){
            vertex_set->insert(v);
            vertex_degree_map->insert({v, v_vertex->get_degree()});
            vertex_mutex_map->insert({v, make_shared<mutex>()});
        }

        auto current_task_vector = make_shared<vector<shared_ptr<unordered_set<uint32_t>>>>(thread_number);
        for(uint32_t i = 0;i<thread_number;++i){
            current_task_vector->at(i) = make_shared<unordered_set<uint32_t>>();
        }

        auto pool = make_shared<thread_pool>(thread_number);
        auto evicted_set = make_shared<unordered_set<uint32_t>>();
        while(!vertex_set->empty()){
            k_max = k;

            if (!tree->count(k)) {
                k_order->insert({k, make_shared<double_list<uint32_t>>()});
                tree->insert({k, make_shared<extend_list<double, uint32_t>>()});
            }
            for(const auto&u:*vertex_set){
                if(vertex_degree_map->at(u) <= k){
                    evicted_set->insert(u);
                }
            }

            while(!evicted_set->empty()){
                for(const auto&u:*evicted_set){
                    vertex_core_map->insert({u,k});

                    auto u_node = make_shared<double_node<uint32_t>>(u);
                    node_map->insert({u,u_node});
                    k_order->at(k)->right_insert(u_node);

                    tree->at(k)->push_back(u);

                    pool->submit_task([=]{
                        auto thread_id = pool->get_thread_id(std::this_thread::get_id());
                        for(const auto&[v,e]:*G->get_vertex(u)->get_edge_map()){
                            if(!vertex_set->count(v) || evicted_set->count(v)){
                                continue;
                            }

                            vertex_mutex_map->at(v)->lock();
                            --vertex_degree_map->at(v);
                            if(vertex_degree_map->at(v) <= k){
                                current_task_vector->at(thread_id)->insert(v);
                            }
                            vertex_mutex_map->at(v)->unlock();
                        }
                    });
                }
                pool->barrier();
                for(const auto u:*evicted_set){
                    vertex_set->erase(u);
                }
                evicted_set->clear();
                merger_set(current_task_vector, evicted_set, pool);
            }
            ++k;
        }
        return k_max;
    }

    /**
     * @details base core decomposition to get the core number of each vertex
     * @param G
     * @param vertex_core_map
     * @return
     */
    uint32_t basic_core_decomposition::decompose(const shared_ptr<abstract_graph> &G,
                                                 const shared_ptr<unordered_map<uint32_t,uint32_t>> &vertex_core_map,
                                                 const shared_ptr<unordered_map<uint32_t,shared_ptr<double_list<uint32_t>>>>& k_order,
                                                 const shared_ptr<unordered_map<uint32_t, shared_ptr<double_node<uint32_t>>>> &node_map,
                                                 const shared_ptr<gadget::Treap> &tree,
                                                 const shared_ptr<vector<long long>>& root) {
        auto vertex_set = make_shared<unordered_set<uint32_t>>();
        auto vertex_degree_map = make_shared<unordered_map<uint32_t, uint32_t>>();
        for(auto [v,v_vertex]:*G->get_vertex_map()){
            vertex_set->insert(v);
            vertex_degree_map->insert({v, v_vertex->get_degree()});
        }

        uint32_t k_max = 1;
        uint32_t  k = 1;
        auto evicted_set = make_shared<unordered_set<uint32_t>>();
        while(!vertex_set->empty()){
            k_max = k;

            if (!k_order->count(k)) {
                k_order->insert({k,make_shared<double_list<uint32_t>>()});
            }
            for(const auto&u:*vertex_set){
                if(vertex_degree_map->at(u) <= k){
                    evicted_set->insert(u);
                }
            }

            while(!evicted_set->empty()){
                auto u = *evicted_set->begin();
                evicted_set->erase(u);
                vertex_set->erase(u);

                vertex_core_map->insert({u,k});
                auto u_node = make_shared<double_node<uint32_t>>(u);
                k_order->at(k)->right_insert(u_node);
                tree->Insert(u, false,root->at(k));
                node_map->insert({u,u_node});

                for(const auto&[v,e]:*G->get_vertex(u)->get_edge_map()){
                    if(!vertex_set->count(v) || evicted_set->count(v)){
                        continue;
                    }

                    --vertex_degree_map->at(v);
                    if (vertex_degree_map->at(v) <= k) {
                        evicted_set->insert(v);
                    }
                }
            }
            ++k;
        }
        return k_max;
    }


    uint32_t basic_core_decomposition::decompose(const shared_ptr<abstract_graph> &G,
                                                 const shared_ptr<unordered_map<uint32_t,uint32_t>> &vertex_core_map,
                                                 const shared_ptr<unordered_map<uint32_t,shared_ptr<double_list<uint32_t>>>>& k_order,
                                                 const shared_ptr<unordered_map<uint32_t, shared_ptr<double_node<uint32_t>>>> &node_map,
                                                 const shared_ptr<gadget::Treap> &tree,
                                                 const shared_ptr<vector<long long>>& root,
                                                 uint32_t thread_number){
        uint32_t k_max = 1;
        uint32_t  k = 1;

        auto vertex_set = make_shared<unordered_set<uint32_t>>();
        auto vertex_degree_map = make_shared<unordered_map<uint32_t, uint32_t>>();
        auto vertex_mutex_map = make_shared<unordered_map<uint32_t, shared_ptr<mutex>>>();

        for(auto [v,v_vertex]:*G->get_vertex_map()){
            vertex_set->insert(v);
            vertex_degree_map->insert({v, v_vertex->get_degree()});
            vertex_mutex_map->insert({v, make_shared<mutex>()});
        }

        auto current_task_vector = make_shared<vector<shared_ptr<unordered_set<uint32_t>>>>(thread_number);
        for(uint32_t i = 0;i<thread_number;++i){
            current_task_vector->at(i) = make_shared<unordered_set<uint32_t>>();
        }

        auto pool = make_shared<thread_pool>(thread_number);
        auto evicted_set = make_shared<unordered_set<uint32_t>>();
        while(!vertex_set->empty()){
            k_max = k;

            if (!k_order->count(k)) {
                k_order->insert({k,make_shared<double_list<uint32_t>>()});
            }
            for(const auto&u:*vertex_set){
                if(vertex_degree_map->at(u) <= k){
                    evicted_set->insert(u);
                }
            }

            while(!evicted_set->empty()){
                for(const auto&u:*evicted_set){
                    vertex_core_map->insert({u,k});

                    auto u_node = make_shared<double_node<uint32_t>>(u);
                    k_order->at(k)->right_insert(u_node);
                    tree->Insert(u, false,root->at(k));
                    node_map->insert({u,u_node});

                    pool->submit_task([=]{
                        auto thread_id = pool->get_thread_id(std::this_thread::get_id());
                        for(const auto&[v,e]:*G->get_vertex(u)->get_edge_map()){
                            if(!vertex_set->count(v) || evicted_set->count(v)){
                                continue;
                            }

                            vertex_mutex_map->at(v)->lock();
                            --vertex_degree_map->at(v);
                            if(vertex_degree_map->at(v) <= k){
                                current_task_vector->at(thread_id)->insert(v);
                            }
                            vertex_mutex_map->at(v)->unlock();
                        }
                    });
                }
                pool->barrier();
                for(const auto u:*evicted_set){
                    vertex_set->erase(u);
                }
                evicted_set->clear();
                merger_set(current_task_vector, evicted_set, pool);
            }
            ++k;
        }
        return k_max;
    }

    /**
     * @details basic core decomposition to get core degree of each vertex and a core order list
     * @param G
     * @param vertex_core_map
     * @param CD
     * @param core_order
     * @return
     */
    uint32_t basic_core_decomposition::decompose(const shared_ptr<abstract_graph> &G,
                                                 const shared_ptr<unordered_map<uint32_t, uint32_t>> &vertex_core_map,
                                                 const shared_ptr<unordered_map<uint32_t,shared_ptr<extend_list<double,uint32_t>>>> &core_order,
                                                 const shared_ptr<unordered_map<uint32_t, shared_ptr<unordered_map<uint32_t, uint32_t>>>> &CD) {
        auto vertex_set = make_shared<unordered_set<uint32_t>>();
        auto vertex_degree_map = make_shared<unordered_map<uint32_t, uint32_t>>();
        for(auto [v,v_vertex]:*G->get_vertex_map()){
            vertex_set->insert(v);
            vertex_degree_map->insert({v, v_vertex->get_degree()});
        }

        uint32_t k_max = 1;
        uint32_t  k = 1;
        auto evicted_set = make_shared<unordered_set<uint32_t>>();
        while(!vertex_set->empty()){
            k_max = k;

            if (!CD->count(k)) {
                CD->insert({k, make_shared<unordered_map<uint32_t,uint32_t>>()});
                core_order->insert({k, make_shared<extend_list<double, uint32_t>>()});
            }
            for(const auto&u:*vertex_set){
                CD->at(k)->insert({u, vertex_degree_map->at(u)});
                if(vertex_degree_map->at(u) <= k){
                    evicted_set->insert(u);
                }
            }

            while(!evicted_set->empty()){
                auto u = *evicted_set->begin();
                evicted_set->erase(u);

                vertex_set->erase(u);

                vertex_core_map->insert({u,k});
                core_order->at(k)->push_back(u);

                for(const auto&[v,e]:*G->get_vertex(u)->get_edge_map()){
                    if(!vertex_set->count(v) || evicted_set->count(v)){
                        continue;
                    }

                    --vertex_degree_map->at(v);
                    if (vertex_degree_map->at(v) <= k) {
                        evicted_set->insert(v);
                    }
                }
            }
            ++k;
        }
        return k_max;
    }

    uint32_t basic_core_decomposition::decompose(const shared_ptr<abstract_graph> &G,
                                                 const shared_ptr<unordered_map<uint32_t, uint32_t>> &vertex_core_map,
                                                 const shared_ptr<unordered_map<uint32_t,shared_ptr<extend_list<double,uint32_t>>>> &core_order,
                                                 const shared_ptr<unordered_map<uint32_t, shared_ptr<unordered_map<uint32_t, uint32_t>>>> &CD,
                                                 uint32_t thread_number){
        uint32_t k_max = 1;
        uint32_t  k = 1;

        auto vertex_set = make_shared<unordered_set<uint32_t>>();
        auto vertex_degree_map = make_shared<unordered_map<uint32_t, uint32_t>>(G->get_vertex_number() + 1);
        auto vertex_mutex_map = make_shared<unordered_map<uint32_t, shared_ptr<mutex>>>(G->get_vertex_number()+1);

        for(auto [v,v_vertex]:*G->get_vertex_map()){
            vertex_set->insert(v);
            vertex_degree_map->insert({v, v_vertex->get_degree()});
            vertex_mutex_map->insert({v,  make_shared<mutex>()});
        }

        auto current_task_vector = make_shared<vector<shared_ptr<unordered_set<uint32_t>>>>(thread_number);
        for(uint32_t i = 0;i<thread_number;++i){
            current_task_vector->at(i) = make_shared<unordered_set<uint32_t>>();
        }

        auto pool = make_shared<thread_pool>(thread_number);
        auto evicted_set = make_shared<unordered_set<uint32_t>>();
        while(!vertex_set->empty()){
            k_max = k;

            if (!CD->count(k)) {
                CD->insert({k, make_shared<unordered_map<uint32_t,uint32_t>>()});
                core_order->insert({k, make_shared<extend_list<double, uint32_t>>()});
            }
            for(const auto&u:*vertex_set){
                if(vertex_degree_map->at(u) <= k){
                    evicted_set->insert(u);
                }
            }

            while(!evicted_set->empty()){
                for(const auto&u:*evicted_set){
                    vertex_core_map->insert({u,k});
                    core_order->at(k)->push_back(u);

                    pool->submit_task([=]{
                        auto thread_id = pool->get_thread_id(std::this_thread::get_id());
                        for(const auto&[v,e]:*G->get_vertex(u)->get_edge_map()){
                            if(!vertex_set->count(v) || evicted_set->count(v)){
                                continue;
                            }

                            vertex_mutex_map->at(v)->lock();
                            --vertex_degree_map->at(v);
                            if(vertex_degree_map->at(v) <= k){
                                current_task_vector->at(thread_id)->insert(v);
                            }
                            vertex_mutex_map->at(v)->unlock();
                        }
                    });
                }
                pool->barrier();
                for(const auto u:*evicted_set){
                    vertex_set->erase(u);
                }
                evicted_set->clear();
                merger_set(current_task_vector, evicted_set, pool);
            }
            ++k;
        }
        return k_max;

    }

    uint32_t basic_core_decomposition::decompose(const shared_ptr<abstract_graph> &G,
                                                 const shared_ptr<unordered_map<uint32_t, uint32_t>> &vertex_core_map,
                                                 const shared_ptr<unordered_map<uint32_t,shared_ptr<extend_list<double,uint32_t>>>> &core_order,
                                                 const shared_ptr<unordered_map<uint32_t,shared_ptr<extend_node<double,uint32_t>>>>& node_map,
                                                 const shared_ptr<unordered_map<uint32_t, shared_ptr<unordered_map<uint32_t, uint32_t>>>> &CD) {
        auto vertex_set = make_shared<unordered_set<uint32_t>>();
        auto vertex_degree_map = make_shared<unordered_map<uint32_t, uint32_t>>();
        for(auto [v,v_vertex]:*G->get_vertex_map()){
            vertex_set->insert(v);
            vertex_degree_map->insert({v, v_vertex->get_degree()});
        }

        uint32_t k_max = 1;
        uint32_t  k = 1;
        auto evicted_set = make_shared<unordered_set<uint32_t>>();
        while(!vertex_set->empty()){
            k_max = k;

            if (!CD->count(k)) {
                CD->insert({k, make_shared<unordered_map<uint32_t, uint32_t>>()});
                core_order->insert({k, make_shared<extend_list<double, uint32_t>>()});
            }
            for(const auto&u:*vertex_set){
                if(vertex_degree_map->at(u) <= k){
                    evicted_set->insert(u);
                }
            }

            while(!evicted_set->empty()){
                auto u = *evicted_set->begin();
                evicted_set->erase(u);
                vertex_set->erase(u);
                vertex_core_map->insert({u,k});

                auto u_node = make_shared<extend_node<double, uint32_t>>(0, u);
                node_map->insert({u, u_node});
                core_order->at(k)->right_insert(u_node);

                for(const auto&[v,e]:*G->get_vertex(u)->get_edge_map()){
                    if(!vertex_set->count(v) || evicted_set->count(v)){
                        continue;
                    }

                    --vertex_degree_map->at(v);
                    if (vertex_degree_map->at(v) <= k) {
                        evicted_set->insert(v);
                    }
                }
            }
            ++k;
        }
        return k_max;
    }

    uint32_t basic_core_decomposition::decompose(const shared_ptr<abstract_graph> &G,
                                                 const shared_ptr<unordered_map<uint32_t, uint32_t>> &vertex_core_map,
                                                 const shared_ptr<unordered_map<uint32_t,shared_ptr<extend_list<double,uint32_t>>>> &core_order,
                                                 const shared_ptr<unordered_map<uint32_t,shared_ptr<extend_node<double,uint32_t>>>>& node_map,
                                                 const shared_ptr<unordered_map<uint32_t, shared_ptr<unordered_map<uint32_t, uint32_t>>>> &CD,
                                                 uint32_t thread_number){
        uint32_t k_max = 1;
        uint32_t  k = 1;

        auto vertex_set = make_shared<unordered_set<uint32_t>>();
        auto vertex_degree_map = make_shared<unordered_map<uint32_t, uint32_t>>();
        auto vertex_mutex_map = make_shared<unordered_map<uint32_t, shared_ptr<mutex>>>();

        for(auto [v,v_vertex]:*G->get_vertex_map()){
            vertex_set->insert(v);
            vertex_degree_map->insert({v, v_vertex->get_degree()});
            vertex_mutex_map->insert({v, make_shared<mutex>()});
        }

        auto current_task_vector = make_shared<vector<shared_ptr<unordered_set<uint32_t>>>>(thread_number);
        for(uint32_t i = 0;i<thread_number;++i){
            current_task_vector->at(i) = make_shared<unordered_set<uint32_t>>();
        }

        auto pool = make_shared<thread_pool>(thread_number);
        auto evicted_set = make_shared<unordered_set<uint32_t>>();
        while(!vertex_set->empty()){
            k_max = k;

            if (!CD->count(k)) {
                CD->insert({k, make_shared<unordered_map<uint32_t, uint32_t>>()});
                core_order->insert({k, make_shared<extend_list<double, uint32_t>>()});
            }
            for(const auto&u:*vertex_set){
                if(vertex_degree_map->at(u) <= k){
                    evicted_set->insert(u);
                }
            }

            while(!evicted_set->empty()){
                for(const auto&u:*evicted_set){
                    vertex_core_map->insert({u,k});

                    auto u_node = make_shared<extend_node<double, uint32_t>>(0, u);
                    node_map->insert({u, u_node});
                    core_order->at(k)->right_insert(u_node);

                    pool->submit_task([=]{
                        auto thread_id = pool->get_thread_id(std::this_thread::get_id());
                        for(const auto&[v,e]:*G->get_vertex(u)->get_edge_map()){
                            if(!vertex_set->count(v) || evicted_set->count(v)){
                                continue;
                            }

                            vertex_mutex_map->at(v)->lock();
                            --vertex_degree_map->at(v);
                            if(vertex_degree_map->at(v) <= k){
                                current_task_vector->at(thread_id)->insert(v);
                            }
                            vertex_mutex_map->at(v)->unlock();
                        }
                    });
                }
                pool->barrier();
                for(const auto u:*evicted_set){
                    vertex_set->erase(u);
                }
                evicted_set->clear();
                merger_set(current_task_vector, evicted_set, pool);
            }
            ++k;
        }
        return k_max;

    }
    /**
     * @brief give the core degree of a given vertex with k
     * @param G
     * @param vertex_core_map
     * @param v
     * @param k
     * @return
     */
    uint32_t basic_core_decomposition::get_core_degree(const shared_ptr<abstract_graph> &G,
                                                                 const shared_ptr<unordered_map<uint32_t, uint32_t>> &vertex_core_map,
                                                                 uint32_t v,
                                                                 uint32_t k) {
        uint32_t core_degree = 0;
        auto v_Vertex = G->get_vertex(v);
        for (const auto &[w,e]:*v_Vertex->get_edge_map()) {
            if (vertex_core_map->at(w) >= k) {
                ++core_degree;
            }
        }
        return core_degree;
    }


    void basic_core_decomposition::merger_set(const shared_ptr<vector<shared_ptr<unordered_set<uint32_t>>>> &input_set_vector,
                                            const shared_ptr<unordered_set<uint32_t>> &output_set,
                                            const shared_ptr<thread_pool> &pool){
        for(auto i = 1; i < input_set_vector->size();i *= 2){
            for(uint32_t j = 0; j + i < input_set_vector->size();j += i * 2){
                pool->submit_task([=]{
                    auto set1 = input_set_vector->at(j);
                    auto set2 = input_set_vector->at(j + i);
                    set1->merge(*set2);
                    set2->clear();
                });
            }
            pool->barrier();
        }
        swap(*output_set,*input_set_vector->at(0));
    }
}



