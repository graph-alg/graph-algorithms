#include "wing/basic_wing_decomposition.h"

namespace scnu
{
    void basic_wing_decomposition::init(const shared_ptr<abstract_bipartite_graph> &G,
                                        const shared_ptr<unordered_map<shared_ptr<abstract_bipartite_edge>, shared_ptr<mutex>>> &edge_mutex_map,
                                        const shared_ptr<unordered_map<shared_ptr<abstract_bipartite_edge>, uint32_t>> &edge_rank_map,
                                        const shared_ptr<unordered_map<shared_ptr<abstract_bipartite_edge>, uint32_t>> &edge_support_map,
                                        const shared_ptr<unordered_map<shared_ptr<abstract_bipartite_edge>, uint32_t>> &edge_wing_map,
                                        const shared_ptr<thread_pool> &pool)
    {
        auto edge_set = G->get_edge_set();
        pool->submit_task([=]{
            for(const auto &e:*edge_set){
                edge_mutex_map->insert({e, make_shared<mutex>()});
            }
        });

        pool->submit_task([=]{
            for(const auto &e:*edge_set){
                edge_rank_map->insert({e, UINT32_MAX});
            }
        });

        pool->submit_task([=]{
            for(const auto &e:*edge_set){
                edge_support_map->insert({e, 0});
            }
        });

        pool->submit_task([=]{
            for(const auto &e:*edge_set){
                edge_wing_map->insert({e, 0});
            }
        });

        pool->barrier();
    }

    void basic_wing_decomposition::init(const shared_ptr<abstract_bipartite_graph> &G,
                                        const shared_ptr<unordered_map<shared_ptr<abstract_bipartite_edge>, shared_ptr<mutex>>> &edge_mutex_map,
                                        const shared_ptr<unordered_map<shared_ptr<abstract_bipartite_edge>, uint32_t>> &edge_rank_map,
                                        const shared_ptr<unordered_map<shared_ptr<abstract_bipartite_edge>, uint32_t>> &edge_support_map,
                                        const shared_ptr<unordered_map<shared_ptr<abstract_bipartite_edge>, uint32_t>> &edge_wing_map,
                                        const shared_ptr<unordered_map<shared_ptr<abstract_bipartite_edge>, uint32_t>> &WS,
                                        const shared_ptr<thread_pool> &pool){
        auto edge_set = G->get_edge_set();
        pool->submit_task([=]{
            for(const auto &e:*edge_set){
                edge_mutex_map->insert({e, make_shared<mutex>()});
            }
        });

        pool->submit_task([=]{
            for(const auto &e:*edge_set){
                edge_rank_map->insert({e, UINT32_MAX});
            }
        });

        pool->submit_task([=]{
            for(const auto &e:*edge_set){
                edge_support_map->insert({e, 0});
            }
        });

        pool->submit_task([=]{
            for(const auto &e:*edge_set){
                edge_wing_map->insert({e, 0});
            }
        });

        pool->submit_task([=]{
            for(const auto &e:*edge_set){
                WS->insert({e, 0});
            }
        });

        pool->barrier();
    }

    void basic_wing_decomposition::init(const shared_ptr<abstract_bipartite_graph> &G,
                                        const shared_ptr<unordered_map<shared_ptr<abstract_bipartite_edge>, shared_ptr<mutex>>> &edge_mutex_map,
                                        const shared_ptr<unordered_map<shared_ptr<abstract_bipartite_edge>, uint32_t>> &edge_rank_map,
                                        const shared_ptr<unordered_map<shared_ptr<abstract_bipartite_edge>, uint32_t>> &edge_support_map,
                                        const shared_ptr<unordered_map<shared_ptr<abstract_bipartite_edge>, uint32_t>> &edge_wing_map,
                                        const shared_ptr<unordered_map<shared_ptr<abstract_bipartite_edge>, uint32_t>> &WS,
                                        const shared_ptr<unordered_map<shared_ptr<abstract_bipartite_edge>, uint32_t>> &rem,
                                        const shared_ptr<thread_pool> &pool){
        auto edge_set = G->get_edge_set();
        pool->submit_task([=]{
            for(const auto &e:*edge_set){
                edge_mutex_map->insert({e, make_shared<mutex>()});
            }
        });

        pool->submit_task([=]{
            for(const auto &e:*edge_set){
                edge_rank_map->insert({e, UINT32_MAX});
            }
        });

        pool->submit_task([=]{
            for(const auto &e:*edge_set){
                edge_support_map->insert({e, 0});
            }
        });

        pool->submit_task([=]{
            for(const auto &e:*edge_set){
                edge_wing_map->insert({e, 0});
            }
        });

        pool->submit_task([=]{
            for(const auto &e:*edge_set){
                WS->insert({e, 0});
            }
        });

        pool->submit_task([=]{
            for(const auto &e:*edge_set){
                rem->insert({e, 0});
            }
        });

        pool->barrier();
    }

    /**
   * @details decompose a given G
   * @param G
   * @param edge_wing_map
   * @return
   */
    uint32_t basic_wing_decomposition::decompose(const shared_ptr<abstract_bipartite_graph> &G,
                                                 const shared_ptr<unordered_map<shared_ptr<abstract_bipartite_edge>,shared_ptr<mutex>>>& edge_mutex_map,
                                                 const shared_ptr<unordered_map<shared_ptr<abstract_bipartite_edge>, uint32_t>> &edge_rank_map,
                                                 const shared_ptr<unordered_map<shared_ptr<abstract_bipartite_edge>, uint32_t>> &edge_support_map,
                                                 const shared_ptr<unordered_map<shared_ptr<abstract_bipartite_edge>, uint32_t>> &edge_wing_map,
                                                 const shared_ptr<thread_pool>& pool) {

        auto edge_set = G->get_edge_set();
        auto vertex_priority_map = make_shared<unordered_map<uint32_t, uint32_t>>();
        vertex_priority_computation(G,vertex_priority_map);
        edge_support_computation(G, edge_set, edge_mutex_map,edge_rank_map, edge_support_map, vertex_priority_map,pool);
        
        uint32_t k = 0;
        uint32_t k_max = 0;

        auto rank_id = make_shared<uint32_t>(0);
        
        auto thread_number = pool->get_thread_number();
        auto global_mutex = make_shared<mutex>();

        while (!edge_set->empty()) {
            auto current_edge_set = make_shared<unordered_set<shared_ptr<abstract_bipartite_edge>>>();
            {
                auto location_vector = pool->split_task(edge_set);
                for(uint32_t i = 0; i < thread_number; ++i){
                    pool->submit_task([=]{
                        auto sub_current_task_set = make_shared<unordered_set<shared_ptr<abstract_bipartite_edge>>>();
                        
                        auto &sub_begin = *location_vector->at(i);
                        auto &sub_end = *location_vector->at(i + 1);
                        
                        for(auto iter = sub_begin; iter!=sub_end;++iter){
                            auto &e = *iter;

                            if (edge_support_map->at(e) <= k) {
                                sub_current_task_set->insert(e);
                            }
                        }
                        
                        global_mutex->lock();
                        current_edge_set->merge(*sub_current_task_set);
                        global_mutex->unlock();
                    });
                }
                 pool->barrier();
            }
            
            if(!current_edge_set->empty())
            {
                k_max = k;
            }

            while (!current_edge_set->empty()) {
                set_edge_rank(current_edge_set, edge_rank_map, rank_id);
                auto next_edge_set = make_shared<unordered_set<shared_ptr<abstract_bipartite_edge>>>();
                auto location_vector = pool->split_task(current_edge_set);
                for(uint32_t i = 0; i < thread_number; ++i){
                    pool->submit_task([=]{
                        auto sub_next_edge_set = make_shared<unordered_set<shared_ptr<abstract_bipartite_edge>>>();
                        
                        auto &sub_begin = *location_vector->at(i);
                        auto &sub_end = *location_vector->at(i + 1);
                        
                        for (auto iter = sub_begin; iter != sub_end; ++iter) {
                            auto &e1 = *iter;

                            auto l1 = e1->get_left_vertex_id();
                            auto r1 = e1->get_right_vertex_id();

                            auto l1_vertex = G->get_left_vertex(l1);
                            auto r1_vertex = G->get_right_vertex(r1);

                            for (const auto&[r2, e2]: *l1_vertex->get_edge_map()) {
                                if (r2 == r1 || edge_rank_map->at(e2) < edge_rank_map->at(e1)) {
                                    continue;
                                }

                                for (const auto&[l2, e3]: *r1_vertex->get_edge_map()) {
                                    if (l2 == l1 || edge_rank_map->at(e3) < edge_rank_map->at(e1)) {
                                        continue;
                                    }


                                    auto e4 = G->get_edge(l2, r2);
                                    if (!e4 || edge_rank_map->at(e4) < edge_rank_map->at(e1)) {
                                        continue;
                                    }

                                    if(edge_support_map->at(e2)>k){
                                        edge_mutex_map->at(e2)->lock();
                                        edge_support_map->at(e2)--;
                                        edge_mutex_map->at(e2)->unlock();

                                        if (edge_support_map->at(e2) <= k) {
                                            sub_next_edge_set->insert(e2);
                                        }
                                    }

                                    if(edge_support_map->at(e3)>k){
                                        edge_mutex_map->at(e3)->lock();
                                        edge_support_map->at(e3)--;
                                        edge_mutex_map->at(e3)->unlock();

                                        if (edge_support_map->at(e3) <= k) {
                                            sub_next_edge_set->insert(e3);
                                        }
                                    }


                                    if(edge_support_map->at(e4)>k){
                                        edge_mutex_map->at(e4)->lock();
                                        edge_support_map->at(e4)--;
                                        edge_mutex_map->at(e4)->unlock();

                                        if (edge_support_map->at(e4) <= k) {
                                            sub_next_edge_set->insert(e4);
                                        }
                                    }
                                }
                            }
                        }
                        
                        global_mutex->lock();
                        next_edge_set->merge(*sub_next_edge_set);
                        global_mutex->unlock();
                    });
                }
                for(const auto&e:*current_edge_set){
                    edge_set->erase(e);
                    edge_wing_map->at(e) = k;
                }
                pool->barrier();
                swap(*current_edge_set, *next_edge_set);
            }
            ++k;
        }
        return k_max;
    }

    uint32_t basic_wing_decomposition::decompose(const shared_ptr<abstract_bipartite_graph> &G,
                                                 const shared_ptr<unordered_map<shared_ptr<abstract_bipartite_edge>,shared_ptr<mutex>>>& edge_mutex_map,
                                                 const shared_ptr<unordered_map<shared_ptr<abstract_bipartite_edge>, uint32_t>> &edge_rank_map,
                                                 const shared_ptr<unordered_map<shared_ptr<abstract_bipartite_edge>, uint32_t>> &edge_support_map,
                                                 const shared_ptr<unordered_map<shared_ptr<abstract_bipartite_edge>, uint32_t>> &edge_wing_map,
                                                 const shared_ptr<unordered_map<shared_ptr<abstract_bipartite_edge>, uint32_t>> &edge_wing_support_map,
                                                 const shared_ptr<thread_pool>& pool) {
        auto thread_number = pool->get_thread_number();
        auto global_mutex = make_shared<mutex>();

        auto edge_set = G->get_edge_set();
        auto vertex_priority_map = make_shared<unordered_map<uint32_t, uint32_t>>();
        vertex_priority_computation(G, vertex_priority_map);
        edge_support_computation(G, edge_set, edge_mutex_map,edge_rank_map,  edge_support_map, vertex_priority_map,pool);

        uint32_t k = 0;
        uint32_t k_max = 0;

        auto rank_id = make_shared<uint32_t>(0);

        while (!edge_set->empty()) {
            auto current_edge_set = make_shared<unordered_set<shared_ptr<abstract_bipartite_edge>>>();
            {
                auto location_vector = pool->split_task(edge_set);
                for(uint32_t i = 0; i < thread_number; ++i){
                    pool->submit_task([=]{
                        auto sub_current_task_set = make_shared<unordered_set<shared_ptr<abstract_bipartite_edge>>>();

                        auto &sub_begin = *location_vector->at(i);
                        auto &sub_end = *location_vector->at(i + 1);

                        for(auto iter = sub_begin; iter!=sub_end;++iter){
                            auto &e = *iter;

                            edge_wing_support_map->at(e) = edge_support_map->at(e);

                            if (edge_support_map->at(e) <= k) {
                                sub_current_task_set->insert(e);
                            }
                        }

                        global_mutex->lock();
                        current_edge_set->merge(*sub_current_task_set);
                        global_mutex->unlock();
                    });
                }
                pool->barrier();
            }

            if(!current_edge_set->empty())
            {
                k_max = k;
            }

            while (!current_edge_set->empty()) {
               set_edge_rank(current_edge_set, edge_rank_map, rank_id);
                auto next_edge_set = make_shared<unordered_set<shared_ptr<abstract_bipartite_edge>>>();
                auto location_vector = pool->split_task(current_edge_set);
                for (uint32_t i = 0; i < thread_number; ++i) {
                    pool->submit_task([=] {
                        auto sub_next_edge_set = make_shared<unordered_set<shared_ptr<abstract_bipartite_edge>>>();

                        auto &sub_begin = *location_vector->at(i);
                        auto &sub_end = *location_vector->at(i + 1);

                        for (auto iter = sub_begin; iter != sub_end; ++iter) {
                            auto &e1 = *iter;

                            auto l1 = e1->get_left_vertex_id();
                            auto r1 = e1->get_right_vertex_id();

                            auto l1_vertex = G->get_left_vertex(l1);
                            auto r1_vertex = G->get_right_vertex(r1);

                            for (const auto&[r2, e2]: *l1_vertex->get_edge_map()) {
                                if (r2 == r1 || edge_rank_map->at(e2) < edge_rank_map->at(e1)) {
                                    continue;
                                }

                                for (const auto&[l2, e3]: *r1_vertex->get_edge_map()) {
                                    if (l2 == l1 || edge_rank_map->at(e3) < edge_rank_map->at(e1)) {
                                        continue;
                                    }


                                    auto e4 = G->get_edge(l2, r2);
                                    if (!e4 || edge_rank_map->at(e4) < edge_rank_map->at(e1)) {
                                        continue;
                                    }

                                    if(edge_support_map->at(e2)>k){
                                        edge_mutex_map->at(e2)->lock();
                                        edge_support_map->at(e2)--;
                                        edge_mutex_map->at(e2)->unlock();

                                        if (edge_support_map->at(e2) <= k) {
                                            sub_next_edge_set->insert(e2);
                                        }
                                    }

                                    if(edge_support_map->at(e3)>k){
                                        edge_mutex_map->at(e3)->lock();
                                        edge_support_map->at(e3)--;
                                        edge_mutex_map->at(e3)->unlock();

                                        if (edge_support_map->at(e3) <= k) {
                                            sub_next_edge_set->insert(e3);
                                        }
                                    }


                                    if(edge_support_map->at(e4)>k){
                                        edge_mutex_map->at(e4)->lock();
                                        edge_support_map->at(e4)--;
                                        edge_mutex_map->at(e4)->unlock();

                                        if (edge_support_map->at(e4) <= k) {
                                            sub_next_edge_set->insert(e4);
                                        }
                                    }
                                }
                            }
                        }
                        global_mutex->lock();
                        next_edge_set->merge(*sub_next_edge_set);
                        global_mutex->unlock();
                    });
                }
                for (const auto &e: *current_edge_set) {
                    edge_set->erase(e);
                    edge_wing_map->at(e) = k;
                }
                pool->barrier();
                for (const auto &e: *current_edge_set) {
                    edge_rank_map->at(e) = 0;
                }
                swap(*current_edge_set, *next_edge_set);
            }
            ++k;
        }
        return k_max;
    }


    uint32_t basic_wing_decomposition::decompose(const shared_ptr<abstract_bipartite_graph> &G,
                                                 const shared_ptr<unordered_map<shared_ptr<abstract_bipartite_edge>,shared_ptr<mutex>>>& edge_mutex_map,
                                                 const shared_ptr<unordered_map<shared_ptr<abstract_bipartite_edge>, uint32_t>> &edge_rank_map,
                                                 const shared_ptr<unordered_map<shared_ptr<abstract_bipartite_edge>, uint32_t>> &edge_support_map,
                                                 const shared_ptr<unordered_map<shared_ptr<abstract_bipartite_edge>, uint32_t>> &edge_wing_map,
                                                 const shared_ptr<unordered_map<shared_ptr<abstract_bipartite_edge>, uint32_t>> &WS,
                                                 const shared_ptr<unordered_map<uint32_t, shared_ptr<extend_list<double, shared_ptr<abstract_bipartite_edge>>>>> &wing_order_map,
                                                 const shared_ptr<unordered_map<shared_ptr<abstract_bipartite_edge>, uint32_t>> &rem,
                                                 const shared_ptr<thread_pool>& pool) {
        auto thread_number = pool->get_thread_number();
        auto global_mutex = make_shared<mutex>();

        auto edge_set = G->get_edge_set();
        auto vertex_priority_map = make_shared<unordered_map<uint32_t, uint32_t>>();
        vertex_priority_computation(G, vertex_priority_map);
        edge_support_computation(G, edge_set, edge_mutex_map,edge_rank_map,  edge_support_map, vertex_priority_map,pool);

        uint32_t k = 0;
        uint32_t k_max = 0;

        auto rank_id = make_shared<uint32_t>(0);

        while (!edge_set->empty()) {
            auto current_edge_set = make_shared<unordered_set<shared_ptr<abstract_bipartite_edge>>>();
            {
                auto location_vector = pool->split_task(edge_set);
                for(uint32_t i = 0; i < thread_number; ++i){
                    pool->submit_task([=]{
                        auto sub_current_task_set = make_shared<unordered_set<shared_ptr<abstract_bipartite_edge>>>();

                        auto &sub_begin = *location_vector->at(i);
                        auto &sub_end = *location_vector->at(i + 1);

                        for(auto iter = sub_begin; iter!=sub_end;++iter){
                            auto &e = *iter;

                            WS->at(e) = edge_support_map->at(e);

                            if (edge_support_map->at(e) <= k) {
                                sub_current_task_set->insert(e);
                            }
                        }

                        global_mutex->lock();
                        current_edge_set->merge(*sub_current_task_set);
                        global_mutex->unlock();
                    });
                }
                pool->barrier();
            }

            if(!current_edge_set->empty())
            {
                k_max = k;
            }

            if(!wing_order_map->count(k)){
                wing_order_map->insert({k, make_shared<extend_list<double, shared_ptr<abstract_bipartite_edge>>>()});
            }

            while (!current_edge_set->empty()) {
                set_edge_rank(current_edge_set, edge_rank_map, rank_id);
                auto next_edge_set = make_shared<unordered_set<shared_ptr<abstract_bipartite_edge>>>();
                auto location_vector = pool->split_task(current_edge_set);
                for (uint32_t i = 0; i < thread_number; ++i) {
                    pool->submit_task([=] {
                        auto sub_next_edge_set = make_shared<unordered_set<shared_ptr<abstract_bipartite_edge>>>();

                        auto &sub_begin = *location_vector->at(i);
                        auto &sub_end = *location_vector->at(i + 1);

                        for (auto iter = sub_begin; iter != sub_end; ++iter) {
                            auto &e1 = *iter;

                            auto l1 = e1->get_left_vertex_id();
                            auto r1 = e1->get_right_vertex_id();

                            auto l1_vertex = G->get_left_vertex(l1);
                            auto r1_vertex = G->get_right_vertex(r1);

                            for (const auto &[r2, e2]: *l1_vertex->get_edge_map()) {
                                if (r2 == r1 || edge_rank_map->at(e2) < edge_rank_map->at(e1)) {
                                    continue;
                                }

                                for (const auto &[l2, e3]: *r1_vertex->get_edge_map()) {
                                    if (l2 == l1 || edge_rank_map->at(e3) < edge_rank_map->at(e1)) {
                                        continue;
                                    }


                                    auto e4 = G->get_edge(l2, r2);
                                    if (!e4 || edge_rank_map->at(e4) < edge_rank_map->at(e1)) {
                                        continue;
                                    }

                                    ++rem->at(e1);

                                    if(edge_support_map->at(e2)>k){
                                        edge_mutex_map->at(e2)->lock();
                                        edge_support_map->at(e2)--;
                                        edge_mutex_map->at(e2)->unlock();

                                        if (edge_support_map->at(e2) <= k) {
                                            sub_next_edge_set->insert(e2);
                                        }
                                    }

                                    if(edge_support_map->at(e3)>k){
                                        edge_mutex_map->at(e3)->lock();
                                        edge_support_map->at(e3)--;
                                        edge_mutex_map->at(e3)->unlock();

                                        if (edge_support_map->at(e3) <= k) {
                                            sub_next_edge_set->insert(e3);
                                        }
                                    }


                                    if(edge_support_map->at(e4)>k){
                                        edge_mutex_map->at(e4)->lock();
                                        edge_support_map->at(e4)--;
                                        edge_mutex_map->at(e4)->unlock();

                                        if (edge_support_map->at(e4) <= k) {
                                            sub_next_edge_set->insert(e4);
                                        }
                                    }
                                }
                            }
                        }

                        global_mutex->lock();
                        next_edge_set->merge(*sub_next_edge_set);
                        global_mutex->unlock();
                    });
                }
                pool->barrier();
                for (const auto &e: *current_edge_set) {
                    edge_set->erase(e);
                    edge_wing_map->at(e) = k;
                    wing_order_map->at(k)->push_back(e);
                }
                swap(*current_edge_set, *next_edge_set);
            }
            ++k;
        }
        return k_max;
    }

    void basic_wing_decomposition::edge_support_computation(const shared_ptr<abstract_bipartite_graph> &G,
                                                            const shared_ptr<unordered_set<shared_ptr<abstract_bipartite_edge>>> &edge_set,
                                                            const shared_ptr<unordered_map<shared_ptr<abstract_bipartite_edge>, uint32_t>> &edge_support_map,
                                                            const shared_ptr<unordered_map<uint32_t, uint32_t>> &vertex_priority_map) {
        /**
         * @brief compute wedges from left vertices
         */
        for(const auto&[l1,l1_vertex]:*G->get_left_vertex_map()){
            auto wedge_count_map = make_shared<unordered_map<uint32_t, uint32_t>>();

            for (const auto&[r1, e1]:*l1_vertex->get_edge_map()) {
                if (vertex_priority_map->at(r1) < vertex_priority_map->at(l1)) {
                    auto r1_vertex = G->get_right_vertex(r1);
                    for (const auto&[l2, e2]:*r1_vertex->get_edge_map()) {
                        if (vertex_priority_map->at(l2) < vertex_priority_map->at(l1)) {
                            if (!wedge_count_map->count(l2)) {
                                wedge_count_map->insert({l2, 0});
                            }
                            ++wedge_count_map->at(l2);
                        }
                    }
                }
            }
            for (const auto&[r1, e1]:*l1_vertex->get_edge_map()) {
                if (vertex_priority_map->at(r1) < vertex_priority_map->at(l1)) {
                    auto r1_vertex = G->get_right_vertex(r1);
                    for (const auto&[l2, e2]:*r1_vertex->get_edge_map()) {
                        if (vertex_priority_map->at(l2) < vertex_priority_map->at(l1)) {
                            if (wedge_count_map->at(l2) > 1) {
                                auto delta = wedge_count_map->at(l2) - 1;

                                edge_support_map->at(e1) += delta;
                                edge_support_map->at(e2) += delta;
                            }
                        }
                    }
                }
            }
        }
        /**
         * @brief compute wedges from right vertices
         */
        for(const auto&[r1, r1_vertex]:*G->get_right_vertex_map()){
            auto wedge_count_map = make_shared<unordered_map<uint32_t, uint32_t>>();

            for (const auto &[l1, e1]: *r1_vertex->get_edge_map()) {
                if (vertex_priority_map->at(l1) < vertex_priority_map->at(r1)) {
                    auto l1_vertex = G->get_left_vertex(l1);
                    for (const auto &[r2, e2]: *l1_vertex->get_edge_map()) {
                        if (vertex_priority_map->at(r2) < vertex_priority_map->at(r1)) {
                            if (!wedge_count_map->count(r2)) {
                                wedge_count_map->insert({r2, 0});
                            }
                            ++wedge_count_map->at(r2);
                        }
                    }
                }
            }
            for (const auto &[l1, e1]: *r1_vertex->get_edge_map()) {
                if (vertex_priority_map->at(l1) < vertex_priority_map->at(r1)) {
                    auto l1_vertex = G->get_left_vertex(l1);
                    for (const auto &[r2, e2]: *l1_vertex->get_edge_map()) {
                        if (vertex_priority_map->at(r2) < vertex_priority_map->at(r1)) {
                            if (wedge_count_map->at(r2) > 1) {
                                auto delta = wedge_count_map->at(r2) - 1;

                                edge_support_map->at(e1) += delta;
                                edge_support_map->at(e2) += delta;
                            }
                        }
                    }
                }
            }
        }
    }

    /**
    * @details compute the support (the number of butterflies) of each edge in the given graph
    * @param G
    * @param vertex_priority_map
    * @param thread_count
    * @return
    */
    void basic_wing_decomposition::edge_support_computation(const shared_ptr<abstract_bipartite_graph> &G,
                                                            const shared_ptr<unordered_set<shared_ptr<abstract_bipartite_edge>>> &edge_set,
                                                            const shared_ptr<unordered_map<shared_ptr<abstract_bipartite_edge>, shared_ptr<mutex>>> &edge_mutex_map,
                                                            const shared_ptr<unordered_map<shared_ptr<abstract_bipartite_edge>, uint32_t>> &edge_rank_map,
                                                            const shared_ptr<unordered_map<shared_ptr<abstract_bipartite_edge>, uint32_t>> &edge_support_map,
                                                            const shared_ptr<unordered_map<uint32_t, uint32_t>> &vertex_priority_map,
                                                            const shared_ptr<thread_pool> &pool) {

        auto thread_number = pool->get_thread_number();
        {
            /**
             * @brief compute wedges from left vertices
             */
            auto left_vertex_map = G->get_left_vertex_map();
            auto location_vector = pool->split_task(left_vertex_map);
            for(uint32_t i = 0; i < thread_number; ++i){
                pool->submit_task([=]{
                    auto &sub_begin = *location_vector->at(i);
                    auto &sub_end = *location_vector->at(i + 1);

                    for(auto iter = sub_begin; iter!=sub_end; ++iter){
                        auto& [l1,l1_vertex] = *iter;
                        auto wedge_count_map = make_shared<unordered_map<uint32_t,uint32_t>>();

                        for(const auto&[r1,e1]:*l1_vertex->get_edge_map()){
                            if(vertex_priority_map->at(r1) < vertex_priority_map->at(l1)){
                                auto r1_vertex = G->get_right_vertex(r1);
                                for(const auto&[l2,e2]:*r1_vertex->get_edge_map()){
                                    if(vertex_priority_map->at(l2) < vertex_priority_map->at(l1)){
                                        if(!wedge_count_map->count(l2)){
                                            wedge_count_map->insert({l2, 0});
                                        }
                                        ++wedge_count_map->at(l2);
                                    }
                                }
                            }
                        }
                        for(const auto&[r1,e1]:*l1_vertex->get_edge_map()){
                            if(vertex_priority_map->at(r1) < vertex_priority_map->at(l1)) {
                                auto r1_vertex = G->get_right_vertex(r1);
                                for (const auto&[l2, e2]:*r1_vertex->get_edge_map()) {
                                    if(vertex_priority_map->at(l2) < vertex_priority_map->at(l1)){
                                        if(wedge_count_map->at(l2) > 1){
                                            auto delta = wedge_count_map->at(l2) - 1;

                                            edge_mutex_map->at(e1)->lock();
                                            edge_support_map->at(e1) += delta;
                                            edge_mutex_map->at(e1)->unlock();

                                            edge_mutex_map->at(e2)->lock();
                                            edge_support_map->at(e2) += delta;
                                            edge_mutex_map->at(e2)->unlock();
                                        }
                                    }
                                }
                            }
                        }
                    }
                });
            }
        }


        /**
         * @brief compute wedges from right vertices
         */
        {
            auto right_vertex_map = G->get_right_vertex_map();
            auto location_vector = pool->split_task(right_vertex_map);
            for(uint32_t i = 0; i < thread_number; ++i) {
                pool->submit_task([=] {
                    auto &sub_begin = *location_vector->at(i);
                    auto &sub_end = *location_vector->at(i + 1);

                    for (auto iter = sub_begin; iter != sub_end; ++iter) {
                        auto& [r1, r1_vertex] = *iter;
                        auto wedge_count_map = make_shared<unordered_map<uint32_t,uint32_t>>();

                        for(const auto&[l1, e1]:*r1_vertex->get_edge_map()){
                            if(vertex_priority_map->at(l1) < vertex_priority_map->at(r1)){
                                auto l1_vertex = G->get_left_vertex(l1);
                                for(const auto&[r2,e2]:*l1_vertex->get_edge_map()){
                                    if(vertex_priority_map->at(r2) < vertex_priority_map->at(r1)){
                                        if(!wedge_count_map->count(r2)){
                                            wedge_count_map->insert({r2, 0});
                                        }
                                        ++wedge_count_map->at(r2);
                                    }
                                }
                            }
                        }

                        for(const auto&[l1,e1]:*r1_vertex->get_edge_map()){
                            if(vertex_priority_map->at(l1) < vertex_priority_map->at(r1)) {
                                auto l1_vertex = G->get_left_vertex(l1);
                                for (const auto&[r2, e2]:*l1_vertex->get_edge_map()) {
                                    if(vertex_priority_map->at(r2) < vertex_priority_map->at(r1)){
                                        if(wedge_count_map->at(r2) > 1){
                                            auto delta = wedge_count_map->at(r2)-1;

                                            edge_mutex_map->at(e1)->lock();
                                            edge_support_map->at(e1) += delta;
                                            edge_mutex_map->at(e1)->unlock();

                                            edge_mutex_map->at(e2)->lock();
                                            edge_support_map->at(e2) += delta;
                                            edge_mutex_map->at(e2)->unlock();
                                        }
                                    }
                                }
                            }
                        }
                    }
                });
            }
        }
        pool->barrier();
    }

    void basic_wing_decomposition::vertex_priority_computation(const shared_ptr<abstract_bipartite_graph> &G,
                                                            const shared_ptr<unordered_map<uint32_t, uint32_t>>& vertex_priority_map) {
        auto pair_vector = make_shared<vector<pair<uint32_t,uint32_t>>>();
        for(const auto&[l,l_vertex]:*G->get_left_vertex_map()){
            pair_vector->emplace_back(make_pair(l,l_vertex->get_degree()));
        }
        for(const auto&[r,r_vertex]:*G->get_right_vertex_map()){
            pair_vector->emplace_back(make_pair(r,r_vertex->get_degree()));
        }

        auto pair_compare=[](const pair<uint32_t,uint32_t>&p1,const pair<uint32_t,uint32_t>&p2){
            if (p1.second > p2.second) {
                return false;
            } else if (p1.second == p2.second && p1.first > p2.first) {
                return false;
            }
            return true;
        };
        sort(pair_vector->begin(), pair_vector->end(), pair_compare);
        for (uint32_t i = 0; i < pair_vector->size(); ++i) {
            vertex_priority_map->insert({pair_vector->at(i).first, i + 1});
        }
    }

    void basic_wing_decomposition::edge_support_computation(const shared_ptr<abstract_bipartite_graph> &G,
                                                            const shared_ptr<unordered_set<shared_ptr<abstract_bipartite_edge>>> &edge_set,
                                                            const shared_ptr<unordered_map<shared_ptr<abstract_bipartite_edge>, uint32_t>> &edge_support_map) {
        auto visited_edge_set = make_shared<unordered_set<shared_ptr<abstract_bipartite_edge>>>();
        for (const auto &e1:*edge_set) {
            auto l1 = e1->get_left_vertex_id();
            auto r1 = e1->get_right_vertex_id();
            visited_edge_set->insert(e1);

            for (const auto&[r2, e2]:*G->get_left_vertex(l1)->get_edge_map()) {
                if (r2 == r1 || visited_edge_set->count(e2)) {
                    continue;
                }
                for (const auto&[l2, e3]:*G->get_right_vertex(r1)->get_edge_map()) {
                    if (l2 == l1 || visited_edge_set->count(e3)) {
                        continue;
                    }
                    auto e4 = G->get_edge(l2, r2);
                    if (!e4 || visited_edge_set->count(e4)) {
                        continue;
                    }
                    ++edge_support_map->at(e1);
                    ++edge_support_map->at(e2);
                    ++edge_support_map->at(e3);
                    ++edge_support_map->at(e4);
                }
            }
        }
    }

    void basic_wing_decomposition::set_edge_rank(
            const shared_ptr<unordered_set<shared_ptr<scnu::abstract_bipartite_edge>>> &edge_set,
            const shared_ptr<unordered_map<shared_ptr<scnu::abstract_bipartite_edge>, uint32_t>> &edge_rank_map,
            const shared_ptr<uint32_t> &rank_id) {

        for(const auto &e:*edge_set){
            *rank_id = *rank_id + 1;
            edge_rank_map->at(e) =*rank_id;
        }
    }
}