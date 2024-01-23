
#include "truss/basic_truss_decomposition.h"

namespace scnu{

    void basic_truss_decomposition::init(const shared_ptr<scnu::abstract_graph> &G,
                                         const shared_ptr<unordered_map<shared_ptr<abstract_edge>, uint32_t>>& edge_support_map,
                                         const shared_ptr<unordered_map<shared_ptr<scnu::abstract_edge>, uint32_t>> &edge_truss_map) {
        auto edge_set = G->get_edge_set();
        for(const auto &e:*edge_set){
            edge_support_map->insert({e, 0});
            edge_truss_map->insert({e, 0});
        }
    }

    void basic_truss_decomposition::init(const shared_ptr<scnu::abstract_graph> &G,
                                         const shared_ptr<unordered_map<shared_ptr<abstract_edge>, uint32_t>>& edge_support_map,
                                         const shared_ptr<unordered_map<shared_ptr<scnu::abstract_edge>, uint32_t>> &edge_truss_map,
                                         const shared_ptr<unordered_map<shared_ptr<abstract_edge>, uint32_t>>& edge_truss_support_map) {
        auto edge_set = G->get_edge_set();
        for(const auto &e:*edge_set){
            edge_support_map->insert({e, 0});
            edge_truss_map->insert({e, 0});
            edge_truss_support_map->insert({e, 0});
        }
    }

    void basic_truss_decomposition::init(const shared_ptr<scnu::abstract_graph> &G,
                                         const shared_ptr<unordered_map<shared_ptr<abstract_edge>, shared_ptr<mutex>>>& edge_mutex_map,
                                         const shared_ptr<unordered_map<shared_ptr<abstract_edge>, uint32_t>>& edge_rank_map,
                                         const shared_ptr<unordered_map<shared_ptr<abstract_edge>, uint32_t>>& edge_support_map,
                                         const shared_ptr<unordered_map<shared_ptr<abstract_edge>, uint32_t>> &edge_truss_map,
                                         const shared_ptr<thread_pool>& pool) {
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
            for(const auto &e:*edge_set) {
                edge_truss_map->insert({e, 0});
            }
        });

        pool->barrier();
    }

    void basic_truss_decomposition::init(const shared_ptr<scnu::abstract_graph> &G,
                                         const shared_ptr<unordered_map<shared_ptr<abstract_edge>, shared_ptr<mutex>>>& edge_mutex_map,
                                         const shared_ptr<unordered_map<shared_ptr<abstract_edge>, uint32_t>>& edge_rank_map,
                                         const shared_ptr<unordered_map<shared_ptr<abstract_edge>, uint32_t>>& edge_support_map,
                                         const shared_ptr<unordered_map<shared_ptr<scnu::abstract_edge>, uint32_t>> &edge_truss_map,
                                         const shared_ptr<unordered_map<shared_ptr<abstract_edge>, uint32_t>>& edge_truss_support_map,
                                         const shared_ptr<unordered_map<shared_ptr<abstract_edge>, uint32_t>>& rem,
                                         const shared_ptr<thread_pool>& pool) {
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
            for(const auto &e:*edge_set) {
                edge_truss_map->insert({e, 0});
            }
        });

        pool->submit_task([=]{
            for(const auto &e:*edge_set){
                edge_truss_support_map->insert({e, 0});
            }
        });

        pool->submit_task([=]{
            for(const auto &e:*edge_set){
                rem->insert({e, 0});
            }
        });

        pool->barrier();
    }

    uint32_t basic_truss_decomposition::decompose(const shared_ptr<abstract_graph>& G,
                                                  const shared_ptr<unordered_map<shared_ptr<abstract_edge>, uint32_t>>& edge_support_map,
                                                  const shared_ptr<unordered_map<shared_ptr<scnu::abstract_edge>, uint32_t>> &edge_truss_map) {
        auto edge_set = G->get_edge_set();
        edge_support_computation(G,edge_set,edge_support_map);

        uint32_t k = 2;
        uint32_t k_max = 2;
        while (!edge_set->empty()) {
            k_max = k;
            auto current_edge_set = make_shared<unordered_set<shared_ptr<abstract_edge>>>();
            for (const auto& e:*edge_set) {
                if(edge_support_map->at(e) <= k- 2){
                    current_edge_set->insert(e);
                }
            }
            while (!current_edge_set->empty()) {
                auto next_edge_set = make_shared<unordered_set<shared_ptr<abstract_edge>>>();
                for(const auto &e1:*current_edge_set){

                    edge_set->erase(e1);
                    edge_truss_map->at(e1) = k;

                    auto u = e1->get_source_vertex_id();
                    auto v = e1->get_destination_vertex_id();
                    if(G->get_vertex(u)->get_degree() > G->get_vertex(v)->get_degree()){
                        swap(u,v);
                    }
                    for(const auto&[w,e2]:* G->get_vertex(u)->get_edge_map()){
                        if(v == w ||!edge_set->count(e2)){
                            continue;
                        }
                        auto e3 = G->get_edge(v, w);
                        if(!e3 || !edge_set->count(e3)){
                            continue;
                        }

                        if (edge_support_map->at(e2) > k - 2) {
                            --edge_support_map->at(e2);
                            if (edge_support_map->at(e2) <= k - 2) {
                                next_edge_set->insert(e2);
                            }
                        }

                        if (edge_support_map->at(e3) > k - 2) {
                            --edge_support_map->at(e3);
                            if (edge_support_map->at(e3) <= k - 2) {
                                next_edge_set->insert(e3);
                            }
                        }
                    }
                }
                swap(*current_edge_set, *next_edge_set);
            }
            ++k;
        }
        return k_max;
    }

    uint32_t basic_truss_decomposition::decompose(const shared_ptr<abstract_graph>& G,
                                                  const shared_ptr<unordered_map<shared_ptr<abstract_edge>, uint32_t>>& edge_truss_map,
                                                  const shared_ptr<unordered_map<shared_ptr<abstract_edge>, uint32_t>>& edge_truss_support_map,
                                                  const shared_ptr<unordered_map<uint32_t, shared_ptr<extend_list<int, shared_ptr<abstract_edge>>>>>& truss_order_map) {
        auto edge_set = G->get_edge_set();
        auto edge_support_map = make_shared<unordered_map<shared_ptr<abstract_edge>, uint32_t>>();
        edge_support_computation(G,edge_set,edge_support_map);

        uint32_t k = 2;
        uint32_t k_max = 2;
        while (!edge_set->empty()) {
            k_max = k;
            if(!truss_order_map->count(k)){
                truss_order_map->insert({k, make_shared<extend_list<int, shared_ptr<abstract_edge>>>()});
            }
            auto evicted_set = make_shared<unordered_set<shared_ptr<abstract_edge>>>();
            for (const auto& e:*edge_set) {
                edge_truss_support_map->insert({e,edge_support_map->at(e)});
                if(edge_support_map->at(e) <= k- 2){
                    evicted_set->insert(e);
                }
            }
            while (!evicted_set->empty()) {
                auto e = *evicted_set->begin();
                evicted_set->erase(e);
                edge_set->erase(e);

                edge_truss_map->insert({e, k});

                auto e_node  = make_shared<extend_node<int, shared_ptr<abstract_edge>>>(0, e);
                truss_order_map->at(k)->right_insert(e_node);

                auto u = e->get_source_vertex_id();
                auto v = e->get_destination_vertex_id();
                if(G->get_vertex(u)->get_degree() > G->get_vertex(v)->get_degree()){
                    swap(u,v);
                }
                for(const auto&[w,e1]:* G->get_vertex(u)->get_edge_map()){
                    if(v == w ||!edge_set->count(e1)){
                        continue;
                    }
                    auto e2 = G->get_edge(v,w);
                    if(!e2 || !edge_set->count(e2)){
                        continue;
                    }

                    if (edge_support_map->at(e1) > k - 2) {
                        --edge_support_map->at(e1);
                        if (edge_support_map->at(e1) <= k - 2) {
                            evicted_set->insert(e1);
                        }
                    }

                    if (edge_support_map->at(e2) > k - 2) {
                        --edge_support_map->at(e2);
                        if (edge_support_map->at(e2) <= k - 2) {
                            evicted_set->insert(e2);
                        }
                    }
                }
            }
            ++k;
        }
        return k_max;
    }

    uint32_t basic_truss_decomposition::decompose(const shared_ptr<abstract_graph>& G,
                                                  const shared_ptr<unordered_map<shared_ptr<abstract_edge>, shared_ptr<mutex>>>& edge_mutex_map,
                                                  const shared_ptr<unordered_map<shared_ptr<abstract_edge>, uint32_t>>& edge_rank_map,
                                                  const shared_ptr<unordered_map<shared_ptr<abstract_edge>, uint32_t>>& edge_support_map,
                                                  const shared_ptr<unordered_map<shared_ptr<abstract_edge>, uint32_t>>& edge_truss_map,
                                                  const shared_ptr<thread_pool>& pool) {
        auto thread_number = pool->get_thread_number();
        auto global_mutex = make_shared<mutex>();

        auto edge_set = G->get_edge_set();
        edge_support_computation(G, edge_set, edge_mutex_map, edge_rank_map, edge_support_map,pool);

        uint32_t k_max = 2;
        auto base_rank = make_shared<uint32_t>(0);
        for (uint32_t k = 2;!edge_set->empty(); ++k) {
            k_max = k;
            auto current_edge_set = make_shared<unordered_set<shared_ptr<abstract_edge>>>();
            {
                auto location_vector = pool->split_task(edge_set);
                for(uint32_t i = 0;i < thread_number;++i){
                    pool->submit_task([=]{
                        auto sub_current_edge_set = make_shared<unordered_set<shared_ptr<abstract_edge>>>();

                        auto &sub_begin = *location_vector->at(i);
                        auto &sub_end = *location_vector->at(i + 1);

                        for(auto iter = sub_begin; iter!=sub_end; ++iter){
                            auto &e = *iter;
                            if (edge_support_map->at(e) <= k - 2) {
                                sub_current_edge_set->insert(e);
                            }
                        }

                        global_mutex->lock();
                        current_edge_set->merge(*sub_current_edge_set);
                        global_mutex->unlock();
                    });
                }
                pool->barrier();
            }
            while (!current_edge_set->empty()) {
                assign_rank(current_edge_set, edge_rank_map,base_rank);
                auto next_edge_set = make_shared<unordered_set<shared_ptr<abstract_edge>>>();
                auto location_vector = pool->split_task(current_edge_set);
                for(uint32_t i = 0; i < thread_number; ++i){
                    pool->submit_task([=]{
                        auto sub_next_edge_set = make_shared<unordered_set<shared_ptr<abstract_edge>>>();

                        auto &sub_begin = *location_vector->at(i);
                        auto &sub_end = *location_vector->at(i + 1);

                        for(auto iter = sub_begin; iter!=sub_end; ++iter) {
                            auto &e1 = *iter;

                            auto u = e1->get_source_vertex_id();
                            auto v = e1->get_destination_vertex_id();
                            if (G->get_vertex(u)->get_degree() > G->get_vertex(v)->get_degree()) {
                                swap(u, v);
                            }
                            for (const auto &[w, e2]: *G->get_vertex(u)->get_edge_map()) {
                                if (v == w || edge_rank_map->at(e2) < edge_rank_map->at(e1)) {
                                    continue;
                                }
                                auto e3 = G->get_edge(v, w);
                                if (!e3 || edge_rank_map->at(e3) < edge_rank_map->at(e1)) {
                                    continue;
                                }

                                if (edge_support_map->at(e2) > k - 2) {
                                    edge_mutex_map->at(e2)->lock();
                                    --edge_support_map->at(e2);
                                    edge_mutex_map->at(e2)->unlock();

                                    if (edge_support_map->at(e2) <= k - 2) {
                                        sub_next_edge_set->insert(e2);
                                    }
                                }

                                if (edge_support_map->at(e3) > k - 2) {
                                    edge_mutex_map->at(e3)->lock();
                                    --edge_support_map->at(e3);
                                    edge_mutex_map->at(e3)->unlock();

                                    if (edge_support_map->at(e3) <= k - 2) {
                                        sub_next_edge_set->insert(e3);
                                    }
                                }
                            }
                        }

                        global_mutex->lock();
                        next_edge_set->merge(*sub_next_edge_set);
                        global_mutex->unlock();
                    });
                }
                for(const auto &e:*current_edge_set) {
                    edge_truss_map->at(e) = k;
                    edge_set->erase(e);
                }
                pool->barrier();
                swap(*current_edge_set, *next_edge_set);
            }
        }
        return k_max;
    }

    uint32_t basic_truss_decomposition::decompose(const shared_ptr<abstract_graph>& G,
                                                  const shared_ptr<unordered_map<shared_ptr<abstract_edge>, shared_ptr<mutex>>>& edge_mutex_map,
                                                  const shared_ptr<unordered_map<shared_ptr<abstract_edge>, uint32_t>>& edge_rank_map,
                                                  const shared_ptr<unordered_map<shared_ptr<abstract_edge>, uint32_t>>& edge_support_map,
                                                  const shared_ptr<unordered_map<shared_ptr<abstract_edge>, uint32_t>>& edge_truss_map,
                                                  const shared_ptr<unordered_map<shared_ptr<abstract_edge>, uint32_t>>& edge_truss_support_map,
                                                  const shared_ptr<thread_pool>& pool) {
        auto thread_number = pool->get_thread_number();
        auto global_mutex = make_shared<mutex>();

        auto edge_set = G->get_edge_set();
        edge_support_computation(G, edge_set, edge_mutex_map, edge_rank_map, edge_support_map,pool);

        uint32_t k_max = 2;
        auto base_rank = make_shared<uint32_t>(0);
        for (uint32_t k = 2;!edge_set->empty(); ++k) {
            k_max = k;
            auto current_edge_set = make_shared<unordered_set<shared_ptr<abstract_edge>>>();
            {
                auto location_vector = pool->split_task(edge_set);
                for(uint32_t i = 0;i < thread_number;++i){
                    pool->submit_task([=]{
                        auto sub_current_edge_set = make_shared<unordered_set<shared_ptr<abstract_edge>>>();

                       auto &sub_begin = *location_vector->at(i);
                       auto &sub_end = *location_vector->at(i + 1);

                       for(auto iter = sub_begin; iter!=sub_end; ++iter){
                           auto &e = *iter;
                           if (edge_support_map->at(e) <= k - 2) {
                               sub_current_edge_set->insert(e);
                           }
                       }

                       global_mutex->lock();
                       current_edge_set->merge(*sub_current_edge_set);
                       global_mutex->unlock();
                    });
                }
                pool->barrier();
            }

            while (!current_edge_set->empty()) {
                assign_rank(current_edge_set, edge_rank_map, base_rank);
                auto next_edge_set = make_shared<unordered_set<shared_ptr<abstract_edge>>>();
                auto location_vector = pool->split_task(current_edge_set);
                for(uint32_t i = 0; i < thread_number; ++i){
                    pool->submit_task([=]{
                        auto sub_next_edge_set = make_shared<unordered_set<shared_ptr<abstract_edge>>>();

                        auto &sub_begin = *location_vector->at(i);
                        auto &sub_end = *location_vector->at(i + 1);

                        for(auto iter = sub_begin; iter!=sub_end; ++iter) {
                            auto &e1 = *iter;

                            auto u = e1->get_source_vertex_id();
                            auto v = e1->get_destination_vertex_id();
                            if (G->get_vertex(u)->get_degree() > G->get_vertex(v)->get_degree()) {
                                swap(u, v);
                            }
                            for (const auto &[w, e2]: *G->get_vertex(u)->get_edge_map()) {
                                if (v == w || edge_rank_map->at(e2) < edge_rank_map->at(e1)) {
                                    continue;
                                }
                                auto e3 = G->get_edge(v, w);
                                if (!e3 || edge_rank_map->at(e3) < edge_rank_map->at(e1)) {
                                    continue;
                                }

                                if (edge_support_map->at(e2) > k - 2) {
                                    edge_mutex_map->at(e2)->lock();
                                    --edge_support_map->at(e2);
                                    edge_mutex_map->at(e2)->unlock();

                                    if (edge_support_map->at(e2) <= k - 2 && !current_edge_set->count(e2)) {
                                        sub_next_edge_set->insert(e2);
                                    }
                                }

                                if (edge_support_map->at(e3) > k - 2) {
                                    edge_mutex_map->at(e3)->lock();
                                    --edge_support_map->at(e3);
                                    edge_mutex_map->at(e3)->unlock();

                                    if (edge_support_map->at(e3) <= k - 2 && !current_edge_set->count(e3)) {
                                        sub_next_edge_set->insert(e3);
                                    }
                                }
                            }
                        }

                        global_mutex->lock();
                        next_edge_set->merge(*sub_next_edge_set);
                        global_mutex->unlock();
                    });
                }
                for(const auto &e:*current_edge_set){
                    edge_truss_map->at(e) = k;
                    edge_set->erase(e);
                }
                pool->barrier();
                swap(*current_edge_set, *next_edge_set);
            }
        }
        return k_max;
    }

    uint32_t basic_truss_decomposition::decompose(const shared_ptr<abstract_graph>& G,
                                                  const shared_ptr<unordered_map<shared_ptr<abstract_edge>, shared_ptr<mutex>>>& edge_mutex_map,
                                                  const shared_ptr<unordered_map<shared_ptr<abstract_edge>, uint32_t>>& edge_rank_map,
                                                  const shared_ptr<unordered_map<shared_ptr<abstract_edge>, uint32_t>>& edge_support_map,
                                                  const shared_ptr<unordered_map<shared_ptr<abstract_edge>, uint32_t>>& edge_truss_map,
                                                  const shared_ptr<unordered_map<shared_ptr<abstract_edge>, uint32_t>>& edge_truss_support_map,
                                                  const shared_ptr<unordered_map<uint32_t, shared_ptr<extend_list<int, shared_ptr<abstract_edge>>>>>& truss_order_map,
                                                  const shared_ptr<unordered_map<shared_ptr<abstract_edge>, uint32_t>>& rem) {
        auto edge_set = G->get_edge_set();
        edge_support_computation(G, edge_set, edge_support_map);

        uint32_t k_max = 2;
        auto base_rank = make_shared<uint32_t>(0);
        for (uint32_t k = 2;!edge_set->empty(); ++k) {
            k_max = k;
            if(!truss_order_map->count(k)){
                truss_order_map->insert({k, make_shared<extend_list<int, shared_ptr<abstract_edge>>>()});
            }
            auto current_edge_set = make_shared<unordered_set<shared_ptr<abstract_edge>>>();
            for (const auto &e: *edge_set) {
                edge_truss_support_map->at(e) = edge_support_map->at(e);
                if (edge_support_map->at(e) <= k - 2) {
                    current_edge_set->insert(e);
                }
            }

            while (!current_edge_set->empty()) {
                assign_rank(current_edge_set, edge_rank_map, base_rank);
                auto next_edge_set = make_shared<unordered_set<shared_ptr<abstract_edge>>>();

                for (const auto &e1: *current_edge_set) {

                    auto u = e1->get_source_vertex_id();
                    auto v = e1->get_destination_vertex_id();
                    if (G->get_vertex(u)->get_degree() > G->get_vertex(v)->get_degree()) {
                        swap(u, v);
                    }
                    for (const auto &[w, e2]: *G->get_vertex(u)->get_edge_map()) {
                        if (v == w || edge_rank_map->at(e2) < edge_rank_map->at(e1)) {
                            continue;
                        }
                        auto e3 = G->get_edge(v, w);
                        if (!e3 || edge_rank_map->at(e3) < edge_rank_map->at(e1)) {
                            continue;
                        }
                        /**
                         * @brief compute rem
                         */
                        ++rem->at(e1);

                        if (edge_support_map->at(e2) > k - 2) {
                            edge_mutex_map->at(e2)->lock();
                            --edge_support_map->at(e2);
                            edge_mutex_map->at(e2)->unlock();

                            if (edge_support_map->at(e2) <= k - 2 && !current_edge_set->count(e2)) {
                                next_edge_set->insert(e2);
                            }
                        }

                        if (edge_support_map->at(e3) > k - 2) {
                            edge_mutex_map->at(e3)->lock();
                            --edge_support_map->at(e3);
                            edge_mutex_map->at(e3)->unlock();

                            if (edge_support_map->at(e3) <= k - 2 && !current_edge_set->count(e3)) {
                                next_edge_set->insert(e3);
                            }
                        }
                    }
                }
                for(const auto &e:*current_edge_set){
                    edge_truss_map->at(e) = k;
                    edge_set->erase(e);

                    truss_order_map->at(k)->push_back(e);
                }
                swap(*current_edge_set, *next_edge_set);
            }
        }

        return k_max;
    }


    uint32_t basic_truss_decomposition::decompose(const shared_ptr<abstract_graph>& G,
                                                  const shared_ptr<unordered_map<shared_ptr<abstract_edge>, shared_ptr<mutex>>>& edge_mutex_map,
                                                  const shared_ptr<unordered_map<shared_ptr<abstract_edge>, uint32_t>>& edge_rank_map,
                                                  const shared_ptr<unordered_map<shared_ptr<abstract_edge>, uint32_t>>& edge_support_map,
                                                  const shared_ptr<unordered_map<shared_ptr<abstract_edge>, uint32_t>>& edge_truss_map,
                                                  const shared_ptr<unordered_map<shared_ptr<abstract_edge>, uint32_t>>& edge_truss_support_map,
                                                  const shared_ptr<unordered_map<uint32_t, shared_ptr<extend_list<int, shared_ptr<abstract_edge>>>>>& truss_order_map,
                                                  const shared_ptr<unordered_map<shared_ptr<abstract_edge>, uint32_t>>& rem,
                                                  const shared_ptr<thread_pool>& pool) {
        auto thread_number = pool->get_thread_number();
        auto global_mutex = make_shared<mutex>();

        auto edge_set = G->get_edge_set();
        edge_support_computation(G, edge_set, edge_mutex_map, edge_rank_map, edge_support_map,pool);

        uint32_t k_max = 2;
        auto base_rank = make_shared<uint32_t>(0);
        for (uint32_t k = 2;!edge_set->empty(); ++k) {
            k_max = k;
            if(!truss_order_map->count(k)){
                truss_order_map->insert({k, make_shared<extend_list<int, shared_ptr<abstract_edge>>>()});
            }
            auto current_edge_set = make_shared<unordered_set<shared_ptr<abstract_edge>>>();
            {
                auto location_vector = pool->split_task(edge_set);
                for(uint32_t i = 0;i < thread_number;++i){
                    pool->submit_task([=] {
                        auto sub_current_edge_set = make_shared<unordered_set<shared_ptr<abstract_edge>>>();

                        auto &sub_begin = *location_vector->at(i);
                        auto &sub_end = *location_vector->at(i + 1);

                        for (auto iter = sub_begin; iter != sub_end; ++iter) {
                            auto &e = *iter;
                            edge_truss_support_map->at(e) = edge_support_map->at(e);
                            if (edge_support_map->at(e) <= k - 2) {
                                sub_current_edge_set->insert(e);
                            }
                        }

                        global_mutex->lock();
                        current_edge_set->merge(*sub_current_edge_set);
                        global_mutex->unlock();
                    });
                }
                pool->barrier();
            }

            while (!current_edge_set->empty()) {
                assign_rank(current_edge_set, edge_rank_map, base_rank);
                auto next_edge_set = make_shared<unordered_set<shared_ptr<abstract_edge>>>();
                {
                    auto location_vector = pool->split_task(current_edge_set);
                    for(uint32_t i = 0; i < thread_number; ++i){
                        pool->submit_task([=] {
                            auto sub_next_edge_set = make_shared<unordered_set<shared_ptr<abstract_edge>>>();

                            auto &sub_begin = *location_vector->at(i);
                            auto &sub_end = *location_vector->at(i + 1);

                            for (auto iter = sub_begin; iter != sub_end; ++iter) {
                                auto &e1 = *iter;

                                auto u = e1->get_source_vertex_id();
                                auto v = e1->get_destination_vertex_id();
                                if (G->get_vertex(u)->get_degree() > G->get_vertex(v)->get_degree()) {
                                    swap(u, v);
                                }
                                for (const auto &[w, e2]: *G->get_vertex(u)->get_edge_map()) {
                                    if (v == w || edge_rank_map->at(e2) < edge_rank_map->at(e1)) {
                                        continue;
                                    }
                                    auto e3 = G->get_edge(v, w);
                                    if (!e3 || edge_rank_map->at(e3) < edge_rank_map->at(e1)) {
                                        continue;
                                    }
                                    /**
                                     * @brief compute rem
                                     */
                                    ++rem->at(e1);

                                    if (edge_support_map->at(e2) > k - 2) {
                                        edge_mutex_map->at(e2)->lock();
                                        --edge_support_map->at(e2);
                                        edge_mutex_map->at(e2)->unlock();

                                        if (edge_support_map->at(e2) <= k - 2) {
                                            sub_next_edge_set->insert(e2);
                                        }
                                    }

                                    if (edge_support_map->at(e3) > k - 2) {
                                        edge_mutex_map->at(e3)->lock();
                                        --edge_support_map->at(e3);
                                        edge_mutex_map->at(e3)->unlock();

                                        if (edge_support_map->at(e3) <= k - 2) {
                                            sub_next_edge_set->insert(e3);
                                        }
                                    }
                                }
                            }

                            global_mutex->lock();
                            next_edge_set->merge(*sub_next_edge_set);
                            global_mutex->unlock();
                        });
                    }
                }
                pool->barrier();
                for(const auto &e:*current_edge_set){
                    edge_truss_map->at(e) = k;
                    edge_set->erase(e);

                    truss_order_map->at(k)->push_back(e);
                }
                swap(*current_edge_set, *next_edge_set);
            }
        }

        return k_max;
    }

    void basic_truss_decomposition::assign_rank(const shared_ptr<unordered_set<shared_ptr<abstract_edge>>>& edge_set,
                                                const shared_ptr<unordered_map<shared_ptr<abstract_edge>, uint32_t>> &edge_rank_map,
                                                const shared_ptr<uint32_t>& base_rank) {
        for (const auto &e: *edge_set) {
            *base_rank = *base_rank + 1;
            edge_rank_map->at(e) = *base_rank;
        }
    }

    void basic_truss_decomposition::edge_support_computation(const shared_ptr<abstract_graph> &G,
                                                             const shared_ptr<unordered_set<shared_ptr<abstract_edge>>> &edge_set,
                                                             const shared_ptr<unordered_map<shared_ptr<abstract_edge>, uint32_t>> &edge_support_map) {
        auto visited_set = make_shared<unordered_set<shared_ptr<abstract_edge>>>();
        for(const auto &e1:*edge_set){
            visited_set->insert(e1);

            auto u = e1->get_source_vertex_id();
            auto v = e1->get_destination_vertex_id();

            if(G->get_vertex(u)->get_degree() > G->get_vertex(v)->get_degree()){
                swap(u,v);
            }

            for(const auto &[w,e2]:*G->get_vertex(u)->get_edge_map()){
                if(w == v || visited_set->count(e2)){
                    continue;
                }

                auto e3 = G->get_edge(v,w);
                if(!e3 || visited_set->count(e3))
                {
                    continue;
                }
                ++edge_support_map->at(e1);
                ++edge_support_map->at(e2);
                ++edge_support_map->at(e3);
            }
        }
    }

    void basic_truss_decomposition::edge_support_computation(const shared_ptr<abstract_graph> &G,
                                                             const shared_ptr<unordered_set<shared_ptr<abstract_edge>>> &edge_set,
                                                             const shared_ptr<unordered_map<shared_ptr<abstract_edge>, shared_ptr<mutex>>> &edge_mutex_map,
                                                             const shared_ptr<unordered_map<shared_ptr<abstract_edge>, uint32_t>> &edge_rank_map,
                                                             const shared_ptr<unordered_map<shared_ptr<abstract_edge>, uint32_t>> &edge_support_map,
                                                             const shared_ptr<thread_pool> &pool) {
        uint32_t thread_number = pool->get_thread_number();
        auto base_rank = make_shared<uint32_t>(0);
        assign_rank(edge_set, edge_rank_map, base_rank);

        auto location_vector = pool->split_task(edge_set);

        for(uint32_t i = 0;i < thread_number;++i){
            //pool->submit_task([=]{
               auto &sub_begin = *location_vector->at(i);
               auto &sub_end = *location_vector->at(i + 1);

                for (auto iter = sub_begin; iter != sub_end; ++iter) {
                    auto &e1 = *iter;

                    auto u = e1->get_source_vertex_id();
                    auto v = e1->get_destination_vertex_id();

                    if (G->get_vertex(u)->get_degree() > G->get_vertex(v)->get_degree()) {
                        swap(u, v);
                    }

                    for (const auto &[w, e2]: *G->get_vertex(u)->get_edge_map()) {
                        if (w == v || edge_rank_map->at(e2) < edge_rank_map->at(e1)) {
                            continue;
                        }

                        auto e3 = G->get_edge(v, w);
                        if (!e3 || edge_rank_map->at(e3) < edge_rank_map->at(e1)) {
                            continue;
                        }

                        edge_mutex_map->at(e1)->lock();
                        ++edge_support_map->at(e1);
                        edge_mutex_map->at(e1)->unlock();

                        edge_mutex_map->at(e2)->lock();
                        ++edge_support_map->at(e2);
                        edge_mutex_map->at(e2)->unlock();

                        edge_mutex_map->at(e3)->lock();
                        ++edge_support_map->at(e3);
                        edge_mutex_map->at(e3)->unlock();
                    }
                }
            //});
        }
        pool->barrier();

        for(const auto&e:*edge_set){
            edge_rank_map->at(e) = UINT32_MAX;
        }
    }
}