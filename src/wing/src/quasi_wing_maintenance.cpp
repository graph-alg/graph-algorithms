
#include "wing/quasi_wing_maintenance.h"

namespace scnu {

    void quasi_wing_maintenance::candidate_graph_finding(const shared_ptr<abstract_bipartite_graph> &G,
                                                         const shared_ptr<unordered_map<shared_ptr<abstract_bipartite_edge>,shared_ptr<mutex>>> &edge_mutex_map,
                                                         const shared_ptr<unordered_map<shared_ptr<abstract_bipartite_edge>,uint32_t>>& edge_rank_map,
                                                         const shared_ptr<unordered_map<shared_ptr<abstract_bipartite_edge>, uint32_t>> &edge_wing_map,
                                                         const shared_ptr<unordered_set<shared_ptr<abstract_bipartite_edge>>> &edge_set,
                                                         const shared_ptr<unordered_set<shared_ptr<abstract_bipartite_edge>>> &candidate_edge_set,
                                                         const shared_ptr<unordered_map<shared_ptr<abstract_bipartite_edge>, uint32_t>> &candidate_edge_support_map,
                                                         const shared_ptr<unordered_map<shared_ptr<abstract_bipartite_edge>, uint32_t>> &WS,
                                                         uint32_t k,
                                                         const shared_ptr<thread_pool>& pool) {
        auto clear_set = make_shared<unordered_set<shared_ptr<abstract_bipartite_edge>>>();
        auto thread_number = pool->get_thread_number();
        auto global_mutex = make_shared<mutex>();
        /**
         * @brief find the k-insert graph
         */
        auto current_edge_set = make_shared<unordered_set<shared_ptr<abstract_bipartite_edge>>>();
        {
            auto removed_edge_set = make_shared<unordered_set<shared_ptr<abstract_bipartite_edge>>>();
            auto location_vector = pool->split_task(edge_set);
            for(uint32_t i = 0; i < thread_number; ++i){
                pool->submit_task([=]{
                    auto sub_current_edge_set = make_shared<unordered_set<shared_ptr<abstract_bipartite_edge>>>();
                    auto sub_removed_edge_set = make_shared<unordered_set<shared_ptr<abstract_bipartite_edge>>>();
                    auto sub_clear_set = make_shared<unordered_set<shared_ptr<abstract_bipartite_edge>>>();
                    for(auto iter = *location_vector->at(i); iter != *location_vector->at(i + 1); ++iter){
                        auto &e = *iter;

                        if(edge_wing_map->at(e) < k-1){
                            sub_removed_edge_set->insert(e);
                        }
                        else
                        {
                            if(edge_wing_map->at(e) == k-1){
                                sub_clear_set->insert(e);
                                WS->at(e) = edge_support_computation(G,e,edge_wing_map,k-1);
                                if(WS->at(e) >= k){
                                    sub_current_edge_set->insert(e);
                                }
                            }
                        }
                    }

                    global_mutex->lock();
                    current_edge_set->merge(*sub_current_edge_set);
                    removed_edge_set->merge(*sub_removed_edge_set);
                    clear_set->merge(*sub_clear_set);
                    global_mutex->unlock();
                });
            }
            pool->barrier();
            for(const auto &e:*removed_edge_set){
                edge_set->erase(e);
            }
        }

        /**
         * @brief find a set of candidate edges
         */
        auto rank_id = make_shared<uint32_t>(0);

        auto rectangle_edge_set = make_shared<unordered_set<shared_ptr<abstract_bipartite_edge>>>();
        auto evicted_edge_set =  make_shared<unordered_set<shared_ptr<abstract_bipartite_edge>>>();
        auto invalid_edge_set = make_shared<unordered_set<shared_ptr<abstract_bipartite_edge>>>();
        auto visited_set = make_shared<unordered_set<shared_ptr<abstract_bipartite_edge>>>();

        while (!current_edge_set->empty()) {
            auto next_edge_set = make_shared<unordered_set<shared_ptr<abstract_bipartite_edge>>>();
            auto next_edge_map = make_shared<unordered_map<shared_ptr<abstract_bipartite_edge>, shared_ptr<unordered_set<shared_ptr<abstract_bipartite_edge>>>>>();
            set_edge_rank(current_edge_set, edge_rank_map, rank_id);
            auto location_vector = pool->split_task(current_edge_set);
            for(uint32_t i = 0; i < thread_number; ++i){
                pool->submit_task([=]{
                    auto sub_next_edge_set = make_shared<unordered_set<shared_ptr<abstract_bipartite_edge>>>();
                    auto sub_next_edge_map = make_shared<unordered_map<shared_ptr<abstract_bipartite_edge>, shared_ptr<unordered_set<shared_ptr<abstract_bipartite_edge>>>>>();
                    auto sub_rectangle_set = make_shared<unordered_set<shared_ptr<abstract_bipartite_edge>>>();

                    for(auto iter = *location_vector->at(i); iter != *location_vector->at(i + 1); ++iter){
                        auto &e1 = *iter;
                        sub_rectangle_set->insert(e1);

                        auto l1 = e1->get_left_vertex_id();
                        auto r1 = e1->get_right_vertex_id();


                        auto e1_set = make_shared<unordered_set<shared_ptr<abstract_bipartite_edge>>>();

                        for (const auto&[r2, e2]:*G->get_left_vertex(l1)->get_edge_map()) {
                            if (r2 == r1 || evicted_edge_set->count(e2)
                                || edge_rank_map->at(e2) < edge_rank_map->at(e1)
                                || edge_wing_map->at(e2) < k - 1 ) {
                                continue;
                            }

                            if(edge_wing_map->at(e2) == k - 1){
                                edge_mutex_map->at(e2)->lock();
                                if(WS->at(e2) == 0){
                                    WS->at(e2) = edge_support_computation(G, e2, edge_wing_map, k - 1);
                                }
                                edge_mutex_map->at(e2)->unlock();
                                if(WS->at(e2) < k){
                                    continue;
                                }
                            }

                            for (const auto&[l2, e3]:*G->get_right_vertex(r1)->get_edge_map()) {
                                if (l2 == l1 || evicted_edge_set->count(e3)
                                    || edge_rank_map->at(e3) < edge_rank_map->at(e1)
                                    || edge_wing_map->at(e3) < k - 1 ) {
                                    continue;
                                }

                                if(edge_wing_map->at(e3) == k - 1){
                                    edge_mutex_map->at(e3)->lock();
                                    if(WS->at(e3) == 0){
                                        WS->at(e3) = edge_support_computation(G, e3, edge_wing_map, k - 1);
                                    }
                                    edge_mutex_map->at(e3)->unlock();
                                    if(WS->at(e3) < k){
                                        continue;
                                    }
                                }

                                auto e4 = G->get_edge(l2, r2);
                                if (!e4  || evicted_edge_set->count(e4)
                                    || edge_rank_map->at(e4) < edge_rank_map->at(e1)
                                    || edge_wing_map->at(e4) < k - 1 ) {
                                    continue;
                                }

                                if(edge_wing_map->at(e4) == k - 1){
                                    edge_mutex_map->at(e4)->lock();
                                    if(WS->at(e4) == 0){
                                        WS->at(e4) = edge_support_computation(G, e4, edge_wing_map, k - 1);
                                    }
                                    edge_mutex_map->at(e4)->unlock();
                                    if(WS->at(e4) < k){
                                        continue;
                                    }
                                }


                                edge_mutex_map->at(e1)->lock();
                                ++candidate_edge_support_map->at(e1);
                                edge_mutex_map->at(e1)->unlock();

                                if(current_edge_set->count(e2)){
                                    edge_mutex_map->at(e2)->lock();
                                    ++candidate_edge_support_map->at(e2);
                                    edge_mutex_map->at(e2)->unlock();
                                }

                                if(current_edge_set->count(e3)){
                                    edge_mutex_map->at(e3)->lock();
                                    ++candidate_edge_support_map->at(e3);
                                    edge_mutex_map->at(e3)->unlock();
                                }

                                if(current_edge_set->count(e4)){
                                    edge_mutex_map->at(e4)->lock();
                                    ++candidate_edge_support_map->at(e4);
                                    edge_mutex_map->at(e4)->unlock();
                                }

                                if(edge_wing_map->at(e2) == k-1){
                                    sub_rectangle_set->insert(e2);

                                    if(!visited_set->count(e2) && !current_edge_set->count(e2)){
                                        e1_set->insert(e2);
                                    }
                                }

                                if(edge_wing_map->at(e3) == k-1){
                                    sub_rectangle_set->insert(e3);

                                    if(!visited_set->count(e3) && !current_edge_set->count(e3)){
                                        e1_set->insert(e3);
                                    }
                                }

                                if(edge_wing_map->at(e4) == k-1){
                                    sub_rectangle_set->insert(e4);

                                    if(!visited_set->count(e4) && !current_edge_set->count(e4)){
                                        e1_set->insert(e4);
                                    }
                                }
                            }
                        }
                        sub_next_edge_map->insert({e1, e1_set});
                    }

                    global_mutex->lock();
                    next_edge_map->merge(*sub_next_edge_map);
                    rectangle_edge_set->merge(*sub_rectangle_set);
                    global_mutex->unlock();
                });
            }
            pool->barrier();
            for(const auto &e1:*current_edge_set){
                edge_rank_map->at(e1) = UINT32_MAX;
                candidate_edge_set->insert(e1);
                if (candidate_edge_support_map->at(e1) < k) {
                    next_edge_map->erase(e1);
                    invalid_edge_set->insert(e1);
                }
            }

            remove_unsatisfied_edges(G, edge_mutex_map, edge_rank_map, edge_wing_map, invalid_edge_set, candidate_edge_set, candidate_edge_support_map, rectangle_edge_set, next_edge_map, evicted_edge_set, k, pool);
            for(const auto &[e, e_set]:*next_edge_map){
                next_edge_set->merge(*e_set);
            }

            visited_set->merge(*current_edge_set);
            swap(*current_edge_set, *next_edge_set);
        }


        for(const auto&e:*clear_set){
            WS->at(e) = 0;
        }
    }


    void quasi_wing_maintenance::candidate_graph_finding2(const shared_ptr<abstract_bipartite_graph> &G,
                                                          const shared_ptr<unordered_map<shared_ptr<abstract_bipartite_edge>,shared_ptr<mutex>>> &edge_mutex_map,
                                                          const shared_ptr<unordered_map<shared_ptr<abstract_bipartite_edge>,uint32_t>>& edge_rank_map,
                                                          const shared_ptr<unordered_map<shared_ptr<abstract_bipartite_edge>, uint32_t>> &edge_wing_map,
                                                          const shared_ptr<unordered_set<shared_ptr<abstract_bipartite_edge>>> &edge_set,
                                                          const shared_ptr<unordered_set<shared_ptr<abstract_bipartite_edge>>> &candidate_edge_set,
                                                          const shared_ptr<unordered_map<shared_ptr<abstract_bipartite_edge>, uint32_t>> &candidate_edge_support_map,
                                                          const shared_ptr<unordered_map<shared_ptr<abstract_bipartite_edge>, uint32_t>> &WS,
                                                          uint32_t k,
                                                          const shared_ptr<thread_pool>& pool) {
        auto thread_number = pool->get_thread_number();
        auto global_mutex = make_shared<mutex>();
        /**
         * @brief find the k-insert graph
         */
        auto current_edge_set = make_shared<unordered_set<shared_ptr<abstract_bipartite_edge>>>();
        {
            auto removed_edge_set = make_shared<unordered_set<shared_ptr<abstract_bipartite_edge>>>();
            auto location_vector = pool->split_task(edge_set);
            for(uint32_t i = 0; i < thread_number; ++i){
                pool->submit_task([=]{
                    auto sub_current_edge_set = make_shared<unordered_set<shared_ptr<abstract_bipartite_edge>>>();
                    auto sub_removed_edge_set = make_shared<unordered_set<shared_ptr<abstract_bipartite_edge>>>();

                    auto &sub_begin = *location_vector->at(i);
                    auto &sub_end = *location_vector->at(i + 1);

                    for(auto iter = sub_begin; iter!=sub_end; ++iter){
                        auto &e = *iter;

                        if(edge_wing_map->at(e) < k-1){
                            sub_removed_edge_set->insert(e);
                        }
                        else
                        {
                            if(edge_wing_map->at(e) == k-1 && WS->at(e) >= k){
                                sub_current_edge_set->insert(e);
                            }
                        }
                    }

                    global_mutex->lock();
                    current_edge_set->merge(*sub_current_edge_set);
                    removed_edge_set->merge(*sub_removed_edge_set);
                    global_mutex->unlock();
                });
            }
            pool->barrier();
            for(const auto &e:*removed_edge_set){
                edge_set->erase(e);
            }
        }




        /**
         * @brief find a set of candidate edges
         */
        auto rank_id = make_shared<uint32_t>(0);

        auto rectangle_edge_set = make_shared<unordered_set<shared_ptr<abstract_bipartite_edge>>>();
        auto evicted_edge_set =  make_shared<unordered_set<shared_ptr<abstract_bipartite_edge>>>();
        auto invalid_edge_set = make_shared<unordered_set<shared_ptr<abstract_bipartite_edge>>>();
        auto visited_set = make_shared<unordered_set<shared_ptr<abstract_bipartite_edge>>>();
        while (!current_edge_set->empty()) {
            auto next_edge_set = make_shared<unordered_set<shared_ptr<abstract_bipartite_edge>>>();
            auto next_edge_map = make_shared<unordered_map<shared_ptr<abstract_bipartite_edge>, shared_ptr<unordered_set<shared_ptr<abstract_bipartite_edge>>>>>();
            set_edge_rank(current_edge_set, edge_rank_map, rank_id);
            auto location_vector = pool->split_task(current_edge_set);
            for(uint32_t i = 0; i < thread_number; ++i){
                pool->submit_task([=]{
                    auto sub_next_edge_set = make_shared<unordered_set<shared_ptr<abstract_bipartite_edge>>>();
                    auto sub_next_edge_map = make_shared<unordered_map<shared_ptr<abstract_bipartite_edge>, shared_ptr<unordered_set<shared_ptr<abstract_bipartite_edge>>>>>();
                    auto sub_rectangle_set = make_shared<unordered_set<shared_ptr<abstract_bipartite_edge>>>();

                    for(auto iter = *location_vector->at(i); iter != *location_vector->at(i + 1); ++iter){
                        auto &e1 = *iter;
                        sub_rectangle_set->insert(e1);

                        auto l1 = e1->get_left_vertex_id();
                        auto r1 = e1->get_right_vertex_id();


                        auto e1_set = make_shared<unordered_set<shared_ptr<abstract_bipartite_edge>>>();

                        for (const auto&[r2, e2]:*G->get_left_vertex(l1)->get_edge_map()) {
                            if (r2 == r1 || evicted_edge_set->count(e2)
                                || edge_rank_map->at(e2) < edge_rank_map->at(e1)
                                || edge_wing_map->at(e2) < k - 1 || (edge_wing_map->at(e2) == k -1 && WS->at(e2) < k)) {
                                continue;
                            }

                            for (const auto&[l2, e3]:*G->get_right_vertex(r1)->get_edge_map()) {
                                if (l2 == l1 || evicted_edge_set->count(e3)
                                    || edge_rank_map->at(e3) < edge_rank_map->at(e1)
                                    || edge_wing_map->at(e3) < k - 1 || (edge_wing_map->at(e3) == k -1 && WS->at(e3) < k)) {
                                    continue;
                                }

                                auto e4 = G->get_edge(l2, r2);
                                if (!e4  || evicted_edge_set->count(e4)
                                    || edge_rank_map->at(e4) < edge_rank_map->at(e1)
                                    || edge_wing_map->at(e4) < k - 1 || (edge_wing_map->at(e4) == k -1 && WS->at(e4) < k)) {
                                    continue;
                                }


                                edge_mutex_map->at(e1)->lock();
                                ++candidate_edge_support_map->at(e1);
                                edge_mutex_map->at(e1)->unlock();

                                if(current_edge_set->count(e2)){
                                    edge_mutex_map->at(e2)->lock();
                                    ++candidate_edge_support_map->at(e2);
                                    edge_mutex_map->at(e2)->unlock();
                                }

                                if(current_edge_set->count(e3)){
                                    edge_mutex_map->at(e3)->lock();
                                    ++candidate_edge_support_map->at(e3);
                                    edge_mutex_map->at(e3)->unlock();
                                }

                                if(current_edge_set->count(e4)){
                                    edge_mutex_map->at(e4)->lock();
                                    ++candidate_edge_support_map->at(e4);
                                    edge_mutex_map->at(e4)->unlock();
                                }

                                if(edge_wing_map->at(e2) == k-1){
                                    sub_rectangle_set->insert(e2);

                                    if(!visited_set->count(e2) && !current_edge_set->count(e2)){
                                        e1_set->insert(e2);
                                    }
                                }

                                if(edge_wing_map->at(e3) == k-1){
                                    sub_rectangle_set->insert(e3);

                                    if(!visited_set->count(e3) && !current_edge_set->count(e3)){
                                        e1_set->insert(e3);
                                    }
                                }

                                if(edge_wing_map->at(e4) == k-1){
                                    sub_rectangle_set->insert(e4);

                                    if(!visited_set->count(e4) && !current_edge_set->count(e4)){
                                        e1_set->insert(e4);
                                    }
                                }
                            }
                        }
                        sub_next_edge_map->insert({e1, e1_set});
                    }

                    global_mutex->lock();
                    next_edge_map->merge(*sub_next_edge_map);
                    rectangle_edge_set->merge(*sub_rectangle_set);
                    global_mutex->unlock();
                });
            }
            pool->barrier();
            for(const auto &e1:*current_edge_set){
                edge_rank_map->at(e1) = UINT32_MAX;
                candidate_edge_set->insert(e1);
                if (candidate_edge_support_map->at(e1) < k) {
                    next_edge_map->erase(e1);
                    invalid_edge_set->insert(e1);
                }
            }

            remove_unsatisfied_edges(G, edge_mutex_map, edge_rank_map, edge_wing_map, invalid_edge_set, candidate_edge_set, candidate_edge_support_map, rectangle_edge_set, next_edge_map, evicted_edge_set, k, pool);
            for(const auto &[e, e_set]:*next_edge_map){
                next_edge_set->merge(*e_set);
            }

            visited_set->merge(*current_edge_set);
            swap(*current_edge_set, *next_edge_set);
        }
    }

    void quasi_wing_maintenance::candidate_graph_finding3(const shared_ptr<abstract_bipartite_graph> &G,
                                                          const shared_ptr<BE_index>& bloom_index,
                                                          const shared_ptr<unordered_map<shared_ptr<abstract_bipartite_edge>,shared_ptr<mutex>>> &edge_mutex_map,
                                                          const shared_ptr<unordered_map<shared_ptr<abstract_bipartite_edge>,uint32_t>>& edge_rank_map,
                                                          const shared_ptr<unordered_map<shared_ptr<abstract_bipartite_edge>, uint32_t>> &edge_wing_map,
                                                          const shared_ptr<unordered_set<shared_ptr<abstract_bipartite_edge>>> &edge_set,
                                                          const shared_ptr<unordered_set<shared_ptr<abstract_bipartite_edge>>> &candidate_edge_set,
                                                          const shared_ptr<unordered_map<shared_ptr<abstract_bipartite_edge>, uint32_t>> &candidate_edge_support_map,
                                                          const shared_ptr<unordered_map<shared_ptr<abstract_bipartite_edge>, uint32_t>> &WS,
                                                          uint32_t k,
                                                          const shared_ptr<thread_pool>& pool) {
        auto thread_number = pool->get_thread_number();
        auto global_mutex = make_shared<mutex>();
        /**
         * @brief find the k-insert graph
         */
        auto current_edge_set = make_shared<unordered_set<shared_ptr<abstract_bipartite_edge>>>();
        {
            auto removed_edge_set = make_shared<unordered_set<shared_ptr<abstract_bipartite_edge>>>();
            auto location_vector = pool->split_task(edge_set);
            for(uint32_t i = 0; i < thread_number; ++i){
                pool->submit_task([=]{
                    auto sub_current_edge_set = make_shared<unordered_set<shared_ptr<abstract_bipartite_edge>>>();
                    auto sub_removed_edge_set = make_shared<unordered_set<shared_ptr<abstract_bipartite_edge>>>();

                    for(auto iter = *location_vector->at(i); iter!= *location_vector->at(i + 1); ++iter){
                        auto &e = *iter;

                        if(edge_wing_map->at(e) < k-1){
                            sub_removed_edge_set->insert(e);
                        }
                        else
                        {
                            if(edge_wing_map->at(e) == k-1 && WS->at(e) >= k){
                                sub_current_edge_set->insert(e);
                            }
                        }
                    }

                    global_mutex->lock();
                    current_edge_set->merge(*sub_current_edge_set);
                    removed_edge_set->merge(*sub_removed_edge_set);
                    global_mutex->unlock();
                });
            }
            pool->barrier();
            for(const auto &e:*removed_edge_set){
                edge_set->erase(e);
            }
        }

        /**
         * @brief find a set of candidate edges
         */
        auto rank_id = make_shared<uint32_t>(0);

        auto rectangle_edge_set = make_shared<unordered_set<shared_ptr<abstract_bipartite_edge>>>();
        auto evicted_edge_set =  make_shared<unordered_set<shared_ptr<abstract_bipartite_edge>>>();
        auto invalid_edge_set = make_shared<unordered_set<shared_ptr<abstract_bipartite_edge>>>();
        auto visited_set = make_shared<unordered_set<shared_ptr<abstract_bipartite_edge>>>();
        while (!current_edge_set->empty()) {
            auto next_edge_set = make_shared<unordered_set<shared_ptr<abstract_bipartite_edge>>>();
            auto next_edge_map = make_shared<unordered_map<shared_ptr<abstract_bipartite_edge>, shared_ptr<unordered_set<shared_ptr<abstract_bipartite_edge>>>>>();
            set_edge_rank(current_edge_set, edge_rank_map, rank_id);
            auto location_vector = pool->split_task(current_edge_set);
            for(uint32_t i = 0; i < thread_number; ++i){
                pool->submit_task([=]{
                    auto sub_next_edge_set = make_shared<unordered_set<shared_ptr<abstract_bipartite_edge>>>();
                    auto sub_next_edge_map = make_shared<unordered_map<shared_ptr<abstract_bipartite_edge>, shared_ptr<unordered_set<shared_ptr<abstract_bipartite_edge>>>>>();
                    auto sub_rectangle_set = make_shared<unordered_set<shared_ptr<abstract_bipartite_edge>>>();
                    for(auto iter = *location_vector->at(i); iter != *location_vector->at(i + 1); ++iter){
                        auto &e1 = *iter;
                        sub_rectangle_set->insert(e1);

                        auto l1 = e1->get_left_vertex_id();
                        auto r1 = e1->get_right_vertex_id();
                        auto e1_set = make_shared<unordered_set<shared_ptr<abstract_bipartite_edge>>>();
                        auto bloom_set = bloom_index->get_bloom_set(e1);

                        if(bloom_set){
                            for (const auto &B: *bloom_set) {
                                if (B->get_butterfly_count() > 0) {
                                    auto e2 = B->get_twin(e1);

                                    if ( edge_wing_map->at(e2) < k - 1 || evicted_edge_set->count(e2)
                                         || edge_rank_map->at(e2) < edge_rank_map->at(e1)
                                         || (edge_wing_map->at(e2) == k - 1 && WS->at(e2) < k) ) {
                                        continue;
                                    }

                                    auto sub_visited_set = make_shared<unordered_set<shared_ptr<abstract_bipartite_edge>>>();
                                    for (const auto &[e3, e4]: *B->get_edge_map()) {
                                        if (e3 == e1 || e3 == e2||
                                            evicted_edge_set->count(e3) || evicted_edge_set->count(e4) || sub_visited_set->count(e3)) {
                                            continue;
                                        }

                                        sub_visited_set->insert(e4);

//                                        if ((e3->get_left_vertex_id() == e4->get_left_vertex_id() &&
//                                             e3->get_right_vertex_id() < e4->get_right_vertex_id())
//                                            || (e3->get_right_vertex_id() == e4->get_right_vertex_id() &&
//                                                e3->get_left_vertex_id() < e4->get_left_vertex_id())) {
//                                            continue;
//                                        }

                                        if ( edge_wing_map->at(e3) < k - 1
                                             || edge_rank_map->at(e3) < edge_rank_map->at(e1)
                                             || (edge_wing_map->at(e3) == k - 1 && WS->at(e3) < k)) {
                                            continue;
                                        }

                                        if ( edge_wing_map->at(e4) < k - 1
                                             || edge_rank_map->at(e4) < edge_rank_map->at(e1)
                                             || (edge_wing_map->at(e4) == k - 1 && WS->at(e4) < k)) {
                                            continue;
                                        }

//                                        if(e3->get_left_vertex_id()!=l1 && e3->get_right_vertex_id()!=r1){
//                                            if(sub_visited_set->count(e3)){
//                                                continue;
//                                            }
//                                            sub_visited_set->insert(e3);
//                                        }
//                                        else{
//                                            if(sub_visited_set->count(e4)){
//                                                continue;
//                                            }
//                                            sub_visited_set->insert(e4);
//                                        }


                                        edge_mutex_map->at(e1)->lock();
                                        ++candidate_edge_support_map->at(e1);
                                        edge_mutex_map->at(e1)->unlock();

                                        if(current_edge_set->count(e2)){
                                            edge_mutex_map->at(e2)->lock();
                                            ++candidate_edge_support_map->at(e2);
                                            edge_mutex_map->at(e2)->unlock();
                                        }

                                        if(current_edge_set->count(e3)){
                                            edge_mutex_map->at(e3)->lock();
                                            ++candidate_edge_support_map->at(e3);
                                            edge_mutex_map->at(e3)->unlock();
                                        }

                                        if(current_edge_set->count(e4)){
                                            edge_mutex_map->at(e4)->lock();
                                            ++candidate_edge_support_map->at(e4);
                                            edge_mutex_map->at(e4)->unlock();
                                        }

                                        if(edge_wing_map->at(e2) == k-1){
                                            sub_rectangle_set->insert(e2);

                                            if(!visited_set->count(e2) && !current_edge_set->count(e2)){
                                                e1_set->insert(e2);
                                            }
                                        }

                                        if(edge_wing_map->at(e3) == k-1){
                                            sub_rectangle_set->insert(e3);

                                            if(!visited_set->count(e3) && !current_edge_set->count(e3)){
                                                e1_set->insert(e3);
                                            }
                                        }

                                        if(edge_wing_map->at(e4) == k-1){
                                            sub_rectangle_set->insert(e4);

                                            if(!visited_set->count(e4) && !current_edge_set->count(e4)){
                                                e1_set->insert(e4);
                                            }
                                        }

                                    }
                                }
                            }
                        }
                        sub_next_edge_map->insert({e1, e1_set});
                    }

                    global_mutex->lock();
                    next_edge_map->merge(*sub_next_edge_map);
                    rectangle_edge_set->merge(*sub_rectangle_set);
                    global_mutex->unlock();
                });
            }
            pool->barrier();
            for(const auto &e1:*current_edge_set){
                edge_rank_map->at(e1) = UINT32_MAX;
                candidate_edge_set->insert(e1);
                if (candidate_edge_support_map->at(e1) < k) {
                    next_edge_map->erase(e1);
                    invalid_edge_set->insert(e1);
                }
            }
            remove_unsatisfied_edges(G, bloom_index, edge_mutex_map, edge_rank_map, edge_wing_map, invalid_edge_set, candidate_edge_set, candidate_edge_support_map, rectangle_edge_set, next_edge_map, evicted_edge_set, k, pool);
            for(const auto &[e, e_set]:*next_edge_map){
                next_edge_set->merge(*e_set);
            }

            visited_set->merge(*current_edge_set);
            swap(*current_edge_set, *next_edge_set);
        }

    }

    void quasi_wing_maintenance::candidate_graph_finding4(const shared_ptr<abstract_bipartite_graph> &G,
                                                          const shared_ptr<unordered_map<shared_ptr<abstract_bipartite_edge>,shared_ptr<mutex>>> &edge_mutex_map,
                                                          const shared_ptr<unordered_map<shared_ptr<abstract_bipartite_edge>,uint32_t>>& edge_rank_map,
                                                          const shared_ptr<unordered_map<shared_ptr<abstract_bipartite_edge>, uint32_t>> &edge_wing_map,
                                                          const shared_ptr<unordered_set<shared_ptr<abstract_bipartite_edge>>> &edge_set,
                                                          const shared_ptr<unordered_set<shared_ptr<abstract_bipartite_edge>>> &candidate_edge_set,
                                                          const shared_ptr<unordered_map<shared_ptr<abstract_bipartite_edge>, uint32_t>> &candidate_edge_support_map,
                                                          const shared_ptr<unordered_map<shared_ptr<abstract_bipartite_edge>, uint32_t>> &WS,
                                                          const shared_ptr<unordered_map<uint32_t,shared_ptr<map<uint32_t,shared_ptr<unordered_set<uint32_t>>>>>> &WL,
                                                          uint32_t k,
                                                          const shared_ptr<thread_pool>& pool) {
        auto thread_number = pool->get_thread_number();
        auto global_mutex = make_shared<mutex>();
        /**
         * @brief find the k-insert graph
         */
        auto current_edge_set = make_shared<unordered_set<shared_ptr<abstract_bipartite_edge>>>();
        {
            auto removed_edge_set = make_shared<unordered_set<shared_ptr<abstract_bipartite_edge>>>();
            auto location_vector = pool->split_task(edge_set);
            for(uint32_t i = 0; i < thread_number; ++i){
                pool->submit_task([=]{
                    auto sub_current_edge_set = make_shared<unordered_set<shared_ptr<abstract_bipartite_edge>>>();
                    auto sub_removed_edge_set = make_shared<unordered_set<shared_ptr<abstract_bipartite_edge>>>();

                    for(auto iter = *location_vector->at(i); iter!=*location_vector->at(i + 1); ++iter){
                        auto &e = *iter;

                        if(edge_wing_map->at(e) < k-1){
                            sub_removed_edge_set->insert(e);
                        }
                        else
                        {
                            if(edge_wing_map->at(e) == k-1 && WS->at(e) >= k){
                                sub_current_edge_set->insert(e);
                            }
                        }
                    }

                    global_mutex->lock();
                    current_edge_set->merge(*sub_current_edge_set);
                    removed_edge_set->merge(*sub_removed_edge_set);
                    global_mutex->unlock();
                });
            }
            pool->barrier();
            for(const auto &e:*removed_edge_set){
                edge_set->erase(e);
            }
        }




        /**
         * @brief find a set of candidate edges
         */
        auto rank_id = make_shared<uint32_t>(0);

        auto rectangle_edge_set = make_shared<unordered_set<shared_ptr<abstract_bipartite_edge>>>();
        auto evicted_edge_set =  make_shared<unordered_set<shared_ptr<abstract_bipartite_edge>>>();
        auto invalid_edge_set = make_shared<unordered_set<shared_ptr<abstract_bipartite_edge>>>();
        auto visited_set = make_shared<unordered_set<shared_ptr<abstract_bipartite_edge>>>();
        while (!current_edge_set->empty()) {
            auto next_edge_set = make_shared<unordered_set<shared_ptr<abstract_bipartite_edge>>>();
            auto next_edge_map = make_shared<unordered_map<shared_ptr<abstract_bipartite_edge>, shared_ptr<unordered_set<shared_ptr<abstract_bipartite_edge>>>>>();
            set_edge_rank(current_edge_set, edge_rank_map, rank_id);
            auto location_vector = pool->split_task(current_edge_set);
            for(uint32_t i = 0; i < thread_number; ++i){
                pool->submit_task([=]{
                    auto sub_next_edge_set = make_shared<unordered_set<shared_ptr<abstract_bipartite_edge>>>();
                    auto sub_next_edge_map = make_shared<unordered_map<shared_ptr<abstract_bipartite_edge>, shared_ptr<unordered_set<shared_ptr<abstract_bipartite_edge>>>>>();
                    auto sub_rectangle_set = make_shared<unordered_set<shared_ptr<abstract_bipartite_edge>>>();

                    for(auto iter = *location_vector->at(i); iter != *location_vector->at(i + 1); ++iter){
                        auto &e1 = *iter;
                        sub_rectangle_set->insert(e1);

                        auto l1 = e1->get_left_vertex_id();
                        auto r1 = e1->get_right_vertex_id();


                        auto e1_set = make_shared<unordered_set<shared_ptr<abstract_bipartite_edge>>>();
                        for (const auto &r2: *WL->at(l1)->at(k - 1)) {
                            auto e2 = G->get_edge(l1, r2);
                            if (r2 == r1 || evicted_edge_set->count(e2)
                                || edge_rank_map->at(e2) < edge_rank_map->at(e1)
                                || (WS->at(e2) < k)) {
                                continue;
                            }
                            for (const auto &l2: *WL->at(r1)->at(k - 1)) {
                                auto e3 = G->get_edge(l2, r1);
                                if (l2 == l1 || evicted_edge_set->count(e3)
                                    || edge_rank_map->at(e3) < edge_rank_map->at(e1)
                                    || (WS->at(e3) < k)) {
                                    continue;
                                }
                                auto e4 = G->get_edge(l2, r2);
                                if (!e4 || evicted_edge_set->count(e4)
                                    || edge_rank_map->at(e4) < edge_rank_map->at(e1)
                                    || edge_wing_map->at(e4) < k - 1 ||
                                    (edge_wing_map->at(e4) == k - 1 && WS->at(e4) < k)) {
                                    continue;
                                }
                                edge_mutex_map->at(e1)->lock();
                                ++candidate_edge_support_map->at(e1);
                                edge_mutex_map->at(e1)->unlock();

                                if (current_edge_set->count(e2)) {
                                    edge_mutex_map->at(e2)->lock();
                                    ++candidate_edge_support_map->at(e2);
                                    edge_mutex_map->at(e2)->unlock();
                                }

                                if (current_edge_set->count(e3)) {
                                    edge_mutex_map->at(e3)->lock();
                                    ++candidate_edge_support_map->at(e3);
                                    edge_mutex_map->at(e3)->unlock();
                                }

                                if (current_edge_set->count(e4)) {
                                    edge_mutex_map->at(e4)->lock();
                                    ++candidate_edge_support_map->at(e4);
                                    edge_mutex_map->at(e4)->unlock();
                                }

                                //if (edge_wing_map->at(e2) == k - 1) {
                                sub_rectangle_set->insert(e2);

                                if (!visited_set->count(e2) && !current_edge_set->count(e2)) {
                                    e1_set->insert(e2);
                                }
                                //}

                                //if (edge_wing_map->at(e3) == k - 1) {
                                sub_rectangle_set->insert(e3);

                                if (!visited_set->count(e3) && !current_edge_set->count(e3)) {
                                    e1_set->insert(e3);
                                }
                                //}

                                if (edge_wing_map->at(e4) == k - 1) {
                                    sub_rectangle_set->insert(e4);

                                    if (!visited_set->count(e4) && !current_edge_set->count(e4)) {
                                        e1_set->insert(e4);
                                    }
                                }
                            }
                            for (auto iter3 = WL->at(r1)->lower_bound(k); iter3 != WL->at(r1)->end(); ++iter3) {
                                for (const auto &l2: *iter3->second) {
                                    auto e3 = G->get_edge(l2, r1);
//                                        if (r2 == r1 || evicted_edge_set->count(e2)
//                                            || edge_rank_map->at(e2) < edge_rank_map->at(e1)
//                                            || (WS->at(e2) < k)) {
//                                            continue;
//                                        }
                                    auto e4 = G->get_edge(l2, r2);
                                    if (!e4 || evicted_edge_set->count(e4)
                                        || edge_rank_map->at(e4) < edge_rank_map->at(e1)
                                        || edge_wing_map->at(e4) < k - 1 ||
                                        (edge_wing_map->at(e4) == k - 1 && WS->at(e4) < k)) {
                                        continue;
                                    }
                                    edge_mutex_map->at(e1)->lock();
                                    ++candidate_edge_support_map->at(e1);
                                    edge_mutex_map->at(e1)->unlock();

                                    if (current_edge_set->count(e2)) {
                                        edge_mutex_map->at(e2)->lock();
                                        ++candidate_edge_support_map->at(e2);
                                        edge_mutex_map->at(e2)->unlock();
                                    }

//                                        if (current_edge_set->count(e3)) {
//                                            edge_mutex_map->at(e3)->lock();
//                                            ++candidate_edge_support_map->at(e3);
//                                            edge_mutex_map->at(e3)->unlock();
//                                        }

                                    if (current_edge_set->count(e4)) {
                                        edge_mutex_map->at(e4)->lock();
                                        ++candidate_edge_support_map->at(e4);
                                        edge_mutex_map->at(e4)->unlock();
                                    }

                                    //if (edge_wing_map->at(e2) == k - 1) {
                                    sub_rectangle_set->insert(e2);

                                    if (!visited_set->count(e2) && !current_edge_set->count(e2)) {
                                        e1_set->insert(e2);
                                    }
                                    //}

//                                        if (edge_wing_map->at(e3) == k - 1) {
//                                            sub_rectangle_set->insert(e3);
//
//                                            if (!visited_set->count(e3) && !current_edge_set->count(e3)) {
//                                                e1_set->insert(e3);
//                                            }
//                                        }

                                    if (edge_wing_map->at(e4) == k - 1) {
                                        sub_rectangle_set->insert(e4);

                                        if (!visited_set->count(e4) && !current_edge_set->count(e4)) {
                                            e1_set->insert(e4);
                                        }
                                    }
                                }
                            }
                        }

                        for (auto iter2 = WL->at(l1)->lower_bound(k); iter2 != WL->at(l1)->end(); ++iter2) {
                            for (const auto &r2: *iter2->second) {
                                auto e2 = G->get_edge(l1, r2);
//                                if (r2 == r1 || evicted_edge_set->count(e2)
//                                    || edge_rank_map->at(e2) < edge_rank_map->at(e1)
//                                    || (WS->at(e2) < k)) {
//                                    continue;
//                                }
                                for (const auto &l2: *WL->at(r1)->at(k - 1)) {
                                    auto e3 = G->get_edge(l2, r1);
                                    if (l2 == l1 || evicted_edge_set->count(e3)
                                        || edge_rank_map->at(e3) < edge_rank_map->at(e1) || (WS->at(e3) < k)) {
                                        continue;
                                    }
                                    auto e4 = G->get_edge(l2, r2);
                                    if (!e4 || evicted_edge_set->count(e4)
                                        || edge_rank_map->at(e4) < edge_rank_map->at(e1)
                                        || edge_wing_map->at(e4) < k - 1 ||
                                        (edge_wing_map->at(e4) == k - 1 && WS->at(e4) < k)) {
                                        continue;
                                    }
                                    edge_mutex_map->at(e1)->lock();
                                    ++candidate_edge_support_map->at(e1);
                                    edge_mutex_map->at(e1)->unlock();

//                                        if (current_edge_set->count(e2)) {
//                                            edge_mutex_map->at(e2)->lock();
//                                            ++candidate_edge_support_map->at(e2);
//                                            edge_mutex_map->at(e2)->unlock();
//                                        }

                                    if (current_edge_set->count(e3)) {
                                        edge_mutex_map->at(e3)->lock();
                                        ++candidate_edge_support_map->at(e3);
                                        edge_mutex_map->at(e3)->unlock();
                                    }

                                    if (current_edge_set->count(e4)) {
                                        edge_mutex_map->at(e4)->lock();
                                        ++candidate_edge_support_map->at(e4);
                                        edge_mutex_map->at(e4)->unlock();
                                    }

//                                        if (edge_wing_map->at(e2) == k - 1) {
//                                            sub_rectangle_set->insert(e2);
//
//                                            if (!visited_set->count(e2) && !current_edge_set->count(e2)) {
//                                                e1_set->insert(e2);
//                                            }
//                                        }

                                    //if (edge_wing_map->at(e3) == k - 1) {
                                    sub_rectangle_set->insert(e3);

                                    if (!visited_set->count(e3) && !current_edge_set->count(e3)) {
                                        e1_set->insert(e3);
                                    }
                                    //}

                                    if (edge_wing_map->at(e4) == k - 1) {
                                        sub_rectangle_set->insert(e4);

                                        if (!visited_set->count(e4) && !current_edge_set->count(e4)) {
                                            e1_set->insert(e4);
                                        }
                                    }
                                }
                                for(auto iter3 = WL->at(r1)->lower_bound(k);iter3!=WL->at(r1)->end();++iter3){
                                    for (const auto &l2: *iter3->second) {
                                        auto e3 = G->get_edge(l2, r1);
//                                        if (l2 == l1 || evicted_edge_set->count(e3)
//                                            || edge_rank_map->at(e3) < edge_rank_map->at(e1) || (WS->at(e3) < k)) {
//                                            continue;
//                                        }
                                        auto e4 = G->get_edge(l2, r2);
                                        if (!e4 || evicted_edge_set->count(e4)
                                            || edge_rank_map->at(e4) < edge_rank_map->at(e1)
                                            || edge_wing_map->at(e4) < k - 1 ||
                                            (edge_wing_map->at(e4) == k - 1 && WS->at(e4) < k)) {
                                            continue;
                                        }
                                        edge_mutex_map->at(e1)->lock();
                                        ++candidate_edge_support_map->at(e1);
                                        edge_mutex_map->at(e1)->unlock();

//                                        if (current_edge_set->count(e2)) {
//                                            edge_mutex_map->at(e2)->lock();
//                                            ++candidate_edge_support_map->at(e2);
//                                            edge_mutex_map->at(e2)->unlock();
//                                        }
//
//                                        if (current_edge_set->count(e3)) {
//                                            edge_mutex_map->at(e3)->lock();
//                                            ++candidate_edge_support_map->at(e3);
//                                            edge_mutex_map->at(e3)->unlock();
//                                        }

                                        if (current_edge_set->count(e4)) {
                                            edge_mutex_map->at(e4)->lock();
                                            ++candidate_edge_support_map->at(e4);
                                            edge_mutex_map->at(e4)->unlock();
                                        }

//                                        if (edge_wing_map->at(e2) == k - 1) {
//                                            sub_rectangle_set->insert(e2);
//
//                                            if (!visited_set->count(e2) && !current_edge_set->count(e2)) {
//                                                e1_set->insert(e2);
//                                            }
//                                        }
//
//                                        if (edge_wing_map->at(e3) == k - 1) {
//                                            sub_rectangle_set->insert(e3);
//
//                                            if (!visited_set->count(e3) && !current_edge_set->count(e3)) {
//                                                e1_set->insert(e3);
//                                            }
//                                        }

                                        if (edge_wing_map->at(e4) == k - 1) {
                                            sub_rectangle_set->insert(e4);

                                            if (!visited_set->count(e4) && !current_edge_set->count(e4)) {
                                                e1_set->insert(e4);
                                            }
                                        }
                                    }
                                }
                            }
                        }
                        sub_next_edge_map->insert({e1, e1_set});
                    }

                    global_mutex->lock();
                    next_edge_map->merge(*sub_next_edge_map);
                    rectangle_edge_set->merge(*sub_rectangle_set);
                    global_mutex->unlock();
                });
            }
            pool->barrier();
            for(const auto &e1:*current_edge_set){
                edge_rank_map->at(e1) = UINT32_MAX;
                candidate_edge_set->insert(e1);
                if (candidate_edge_support_map->at(e1) < k) {
                    next_edge_map->erase(e1);
                    invalid_edge_set->insert(e1);
                }
            }

            remove_unsatisfied_edges(G, edge_mutex_map, edge_rank_map, edge_wing_map, invalid_edge_set, candidate_edge_set, candidate_edge_support_map, rectangle_edge_set, next_edge_map, evicted_edge_set, WL, k, pool);
            for(const auto &[e, e_set]:*next_edge_map){
                next_edge_set->merge(*e_set);
            }

            visited_set->merge(*current_edge_set);
            swap(*current_edge_set, *next_edge_set);
        }
    }


    uint32_t quasi_wing_maintenance::edge_support_computation(const shared_ptr<abstract_bipartite_graph> &G,
                                                              const shared_ptr<abstract_bipartite_edge> &e,
                                                              const shared_ptr<unordered_map<shared_ptr<abstract_bipartite_edge>, uint32_t>> &edge_wing_map,
                                                              uint32_t k) {

        const auto &e1 = e;
        uint32_t e1_support = 0;

        auto l1 = e1->get_left_vertex_id();
        auto r1 = e1->get_right_vertex_id();

        for (const auto&[r2, e2]:*G->get_left_vertex(l1)->get_edge_map()) {
            if (r2 == r1 || edge_wing_map->at(e2) < k) {
                continue;
            }
            for (const auto&[l2, e3]:*G->get_right_vertex(r1)->get_edge_map()) {
                if (l2 == l1 || edge_wing_map->at(e3) < k) {
                    continue;
                }
                auto e4 = G->get_edge(l2, r2);
                if (!e4 || edge_wing_map->at(e4) < k) {
                    continue;
                }

                ++e1_support;
            }
        }
        return e1_support;
    }

    /**
     * @details compute the support of edges in edge_set
     * @param G
     * @param edge_mutex_map
     * @param edge_rank_map
     * @param edge_set
     * @param WS
     * @param pool
     * @return
     */
//    void quasi_wing_maintenance::edge_support_computation(const shared_ptr<abstract_bipartite_graph> &G,
//                                                              const shared_ptr<unordered_map<shared_ptr<abstract_bipartite_edge>,shared_ptr<mutex>>> &edge_mutex_map,
//                                                              const shared_ptr<unordered_map<shared_ptr<abstract_bipartite_edge>,uint32_t>> &edge_rank_map,
//                                                              const shared_ptr<unordered_set<shared_ptr<abstract_bipartite_edge>>> &edge_set,
//                                                              const shared_ptr<unordered_map<shared_ptr<abstract_bipartite_edge>,uint32_t>> & WS,
//                                                              const shared_ptr<thread_pool>& pool) {
//        uint32_t rank_id = 0;
//        for(const auto&e1:*edge_set){
//            edge_rank_map->at(e1) = ++rank_id;
//            pool->submit_task([=]{
//                auto l1 = e1->get_left_vertex_id();
//                auto r1 = e1->get_right_vertex_id();
//
//                for (const auto&[r2, e2]:*G->get_left_vertex(l1)->get_edge_map()) {
//                    if (r2 == r1) {
//                        continue;
//                    }
//                    for (const auto&[l2, e3]:*G->get_right_vertex(r1)->get_edge_map()) {
//                        if (l2 == l1) {
//                            continue;
//                        }
//                        auto e4 = G->get_edge(l2, r2);
//                        if (!e4) {
//                            continue;
//                        }
//
//                        edge_mutex_map->at(e1)->lock();
//                        ++WS->at(e1);
//                        edge_mutex_map->at(e1)->unlock();
//
//                        edge_mutex_map->at(e2)->lock();
//                        ++WS->at(e2);
//                        edge_mutex_map->at(e2)->unlock();
//
//                        edge_mutex_map->at(e3)->lock();
//                        ++WS->at(e3);
//                        edge_mutex_map->at(e3)->unlock();
//
//                        edge_mutex_map->at(e4)->lock();
//                        ++WS->at(e4);
//                        edge_mutex_map->at(e4)->unlock();
//                    }
//                }
//            });
//        }
//        pool->barrier();
//    }


    /**
     * @details compute the support of a given edge in k-wing
     * @param G
     * @param e
     * @param edge_wing_map
     * @param WS
     * @param k
     * @return
     */
    uint32_t quasi_wing_maintenance::edge_support_computation(const shared_ptr<abstract_bipartite_graph> &G,
                                                              const shared_ptr<abstract_bipartite_edge> &e,
                                                              const shared_ptr<unordered_map<shared_ptr<abstract_bipartite_edge>, uint32_t>> &edge_wing_map,
                                                              const shared_ptr<unordered_map<shared_ptr<abstract_bipartite_edge>, uint32_t>> &WS,
                                                              uint32_t k) {
        if(WS->at(e) == 0){
            const auto & e1 = e;
            uint32_t e1_support = 0;

            auto l1 = e1->get_left_vertex_id();
            auto r1 = e1->get_right_vertex_id();

            for (const auto&[r2, e2]:*G->get_left_vertex(l1)->get_edge_map()) {
                if (r2 == r1 || edge_wing_map->at(e2) < k) {
                    continue;
                }
                for (const auto&[l2, e3]:*G->get_right_vertex(r1)->get_edge_map()) {
                    if (l2 == l1 || edge_wing_map->at(e3) < k) {
                        continue;
                    }
                    auto e4 = G->get_edge(l2, r2);
                    if (!e4 || edge_wing_map->at(e4) < k) {
                        continue;
                    }

                    ++e1_support;
                }
            }
            WS->at(e1) = e1_support;
        }
        return WS->at(e);
    }

    void quasi_wing_maintenance::edge_support_computation(const shared_ptr<abstract_bipartite_graph> &G,
                                                          const shared_ptr<unordered_set<shared_ptr<abstract_bipartite_edge>>> &edge_set,
                                                          const shared_ptr<unordered_map<shared_ptr<abstract_bipartite_edge>, shared_ptr<mutex>>> &edge_mutex_map,
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

                    for(auto iter = sub_begin; iter !=sub_end; ++iter){
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
                        for(const auto&[r1,e1_edge]:*l1_vertex->get_edge_map()){
                            if(vertex_priority_map->at(r1) < vertex_priority_map->at(l1)) {
                                auto r1_vertex = G->get_right_vertex(r1);
                                for (const auto&[l2, e2]:*r1_vertex->get_edge_map()) {
                                    if(vertex_priority_map->at(l2) < vertex_priority_map->at(l1)){
                                        if(wedge_count_map->at(l2) > 1){
                                            auto delta = wedge_count_map->at(l2) - 1;

                                            edge_mutex_map->at(e1_edge)->lock();
                                            edge_support_map->at(e1_edge) += delta;
                                            edge_mutex_map->at(e1_edge)->unlock();

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

        {
            /**
             * @brief compute wedges from right vertices
             */
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

    void quasi_wing_maintenance::init(const shared_ptr<unordered_map<shared_ptr<abstract_bipartite_edge>, uint32_t>> &edge_wing_map,
                                      const shared_ptr<unordered_map<shared_ptr<abstract_bipartite_edge>, shared_ptr<mutex>>> &edge_mutex_map,
                                      const shared_ptr<unordered_map<shared_ptr<abstract_bipartite_edge>, uint32_t>> &edge_rank_map,
                                      const shared_ptr<unordered_map<shared_ptr<abstract_bipartite_edge>, uint32_t>> &edge_support_map,
                                      const shared_ptr<thread_pool> &pool){
        pool->submit_task([=]{
            for(const auto &[e, wing_number]:*edge_wing_map){
                edge_mutex_map->insert({e , make_shared<mutex>()});
            }
        });

        pool->submit_task([=]{
            for(const auto &[e, wing_number]:*edge_wing_map){
                edge_rank_map->insert({e, UINT32_MAX});
            }
        });

        pool->submit_task([=]{
            for(const auto &[e, wing_number]:*edge_wing_map){
                edge_support_map->insert({e, 0});
            }
        });
        pool->barrier();
    }

    void quasi_wing_maintenance::init(const shared_ptr<unordered_map<shared_ptr<abstract_bipartite_edge>, uint32_t>> &edge_wing_map,
                                      const shared_ptr<unordered_map<shared_ptr<abstract_bipartite_edge>, shared_ptr<mutex>>> &edge_mutex_map,
                                      const shared_ptr<unordered_map<shared_ptr<abstract_bipartite_edge>, uint32_t>> &edge_rank_map,
                                      const shared_ptr<unordered_map<shared_ptr<abstract_bipartite_edge>, uint32_t>> &edge_support_map,
                                      const shared_ptr<unordered_map<shared_ptr<abstract_bipartite_edge>, uint32_t>> &WS,
                                      const shared_ptr<thread_pool> &pool){
        pool->submit_task([=]{
            for(const auto &[e, wing_number]:*edge_wing_map){
                edge_mutex_map->insert({e , make_shared<mutex>()});
            }
        });

        pool->submit_task([=]{
            for(const auto &[e, wing_number]:*edge_wing_map){
                edge_rank_map->insert({e, UINT32_MAX});
            }
        });

        pool->submit_task([=]{
            for(const auto &[e, wing_number]:*edge_wing_map){
                edge_support_map->insert({e, 0});
            }
        });

        pool->submit_task([=]{
            for(const auto &[e, wing_number]:*edge_wing_map){
                WS->insert({e, 0});
            }
        });

        pool->barrier();
    }



    /**
     * @details a basic insert method
     * @param G
     * @param edge_set
     * @param edge_wing_map
     * @param previous_k_max
     * @param thread_number
     */
    void quasi_wing_maintenance::insert(const shared_ptr<abstract_bipartite_graph> &G,
                                        const shared_ptr<unordered_map<shared_ptr<abstract_bipartite_edge>,shared_ptr<mutex>>>& edge_mutex_map,
                                        const shared_ptr<unordered_set<shared_ptr<abstract_bipartite_edge>>> &edge_set,
                                        const shared_ptr<unordered_map<shared_ptr<abstract_bipartite_edge>, uint32_t>>& edge_rank_map,
                                        const shared_ptr<unordered_map<shared_ptr<abstract_bipartite_edge>, uint32_t>>& candidate_edge_support_map,
                                        const shared_ptr<unordered_map<shared_ptr<abstract_bipartite_edge>, uint32_t>> &edge_wing_map,
                                        const shared_ptr<unordered_map<shared_ptr<abstract_bipartite_edge>, uint32_t>>& WS,
                                        const shared_ptr<uint32_t>& previous_k_max,
                                        const shared_ptr<thread_pool>& pool) {
        auto thread_number = pool->get_thread_number();
        auto global_mutex = make_shared<mutex>();

        G->insert_edge_collection(edge_set);
        for(const auto& e:*edge_set){
            edge_rank_map->insert({e,UINT32_MAX});
            edge_mutex_map->insert({e, make_shared<mutex>()});
            candidate_edge_support_map->insert({e,0});
            edge_wing_map->insert({e,0});
            WS->insert({e,0});
        }
        auto affected_edge_set = make_shared<unordered_set<shared_ptr<abstract_bipartite_edge>>>();
        {
            auto rank_id = make_shared<uint32_t>(0);
            set_edge_rank(edge_set, edge_rank_map, rank_id);
            auto location_vector = pool->split_task(edge_set);
            for(uint32_t i = 0; i < thread_number; ++i){
                pool->submit_task([=]{
                    auto sub_affected_set = make_shared<unordered_set<shared_ptr<abstract_bipartite_edge>>>();

                    for(auto iter = *location_vector->at(i); iter != *location_vector->at(i + 1); ++iter){
                        const auto &e1 = *iter;
                        auto l1 = e1->get_left_vertex_id();
                        auto r1 = e1->get_right_vertex_id();

                        for(const auto &[r2,e2]:*G->get_left_vertex(l1)->get_edge_map()){
                            if(r2 == r1 || edge_rank_map->at(e2) < edge_rank_map->at(e1)){
                                continue;
                            }
                            for(const auto&[l2,e3]:*G->get_right_vertex(r1)->get_edge_map()){
                                if(l2 == l1 || edge_rank_map->at(e3) < edge_rank_map->at(e1)){
                                    continue;
                                }

                                auto e4 = G->get_edge(l2,r2);
                                if(!e4 || edge_rank_map->at(e4) < edge_rank_map->at(e1)){
                                    continue;
                                }


                                sub_affected_set->insert(e1);
                                sub_affected_set->insert(e2);
                                sub_affected_set->insert(e3);
                                sub_affected_set->insert(e4);
                            }
                        }
                    }

                    global_mutex->lock();
                    affected_edge_set->merge(*sub_affected_set);
                    global_mutex->unlock();
                });
            }
            pool->barrier();
            for(const auto&e:*affected_edge_set){
                edge_rank_map->at(e) = UINT32_MAX;
                if (edge_wing_map->at(e) == 0) {
                    edge_wing_map->at(e) = 1;
                }
            }
        }


        /**
         * @brief update the remainder k-wings
         */
        uint32_t k = 2;
        while(!affected_edge_set->empty()){
            if (k <= *previous_k_max) {
                auto candidate_edge_set = make_shared<unordered_set<shared_ptr<abstract_bipartite_edge>>>();
                candidate_graph_finding(G, edge_mutex_map, edge_rank_map, edge_wing_map, affected_edge_set, candidate_edge_set, candidate_edge_support_map, WS, k, pool);

                for(const auto &e:*candidate_edge_set){
                    affected_edge_set->insert(e);
                    edge_wing_map->at(e) = k;
                    candidate_edge_support_map->at(e) = 0;
                }
                /**
                  * @brief continue the loop
                  */
                ++k;
            } else {
                auto candidate_edge_set = make_shared<unordered_set<shared_ptr<abstract_bipartite_edge>>>();
                candidate_graph_finding(G, edge_mutex_map, edge_rank_map, edge_wing_map, affected_edge_set, candidate_edge_set, candidate_edge_support_map, WS, k, pool);
                *previous_k_max = partial_wing_decomposition(G, edge_mutex_map, edge_rank_map, edge_wing_map, candidate_edge_set, candidate_edge_support_map, k, pool);
                break;
            }
        }
    }

    /**
    * @details a basic insert method
    * @param G
    * @param edge_set
    * @param edge_wing_map
    * @param previous_k_max
    * @param thread_number
    */
    void quasi_wing_maintenance::insert2(const shared_ptr<abstract_bipartite_graph> &G,
                                         const shared_ptr<unordered_map<shared_ptr<abstract_bipartite_edge>,shared_ptr<mutex>>>& edge_mutex_map,
                                         const shared_ptr<unordered_set<shared_ptr<abstract_bipartite_edge>>> &edge_set,
                                         const shared_ptr<unordered_map<shared_ptr<abstract_bipartite_edge>, uint32_t>>& edge_rank_map,
                                         const shared_ptr<unordered_map<shared_ptr<abstract_bipartite_edge>, uint32_t>>& candidate_edge_support_map,
                                         const shared_ptr<unordered_map<shared_ptr<abstract_bipartite_edge>, uint32_t>> &edge_wing_map,
                                         const shared_ptr<unordered_map<shared_ptr<abstract_bipartite_edge>, uint32_t>>& WS,
                                         const shared_ptr<uint32_t>& previous_k_max,
                                         const shared_ptr<thread_pool>& pool) {
        auto thread_number = pool->get_thread_number();
        auto global_mutex = make_shared<mutex>();

        G->insert_edge_collection(edge_set);
        for(const auto& e:*edge_set){
            edge_rank_map->insert({e,UINT32_MAX});
            edge_mutex_map->insert({e, make_shared<mutex>()});
            candidate_edge_support_map->insert({e,0});
            edge_wing_map->insert({e,0});
            WS->insert({e,0});
        }
        auto affected_edge_set = make_shared<unordered_set<shared_ptr<abstract_bipartite_edge>>>();
        {
            auto rank_id = make_shared<uint32_t>(0);
            set_edge_rank(edge_set, edge_rank_map, rank_id);
            auto location_vector = pool->split_task(edge_set);
            for(uint32_t i = 0; i < thread_number; ++i){
                pool->submit_task([=]{
                    auto sub_affected_set = make_shared<unordered_set<shared_ptr<abstract_bipartite_edge>>>();

                    for(auto iter = *location_vector->at(i); iter != *location_vector->at(i + 1); ++iter){
                        const auto &e1 = *iter;
                        auto l1 = e1->get_left_vertex_id();
                        auto r1 = e1->get_right_vertex_id();

                        for(const auto &[r2,e2]:*G->get_left_vertex(l1)->get_edge_map()){
                            if(r2 == r1 || edge_rank_map->at(e2) < edge_rank_map->at(e1)){
                                continue;
                            }
                            for(const auto&[l2,e3]:*G->get_right_vertex(r1)->get_edge_map()){
                                if(l2 == l1 || edge_rank_map->at(e3) < edge_rank_map->at(e1)){
                                    continue;
                                }

                                auto e4 = G->get_edge(l2,r2);
                                if(!e4 || edge_rank_map->at(e4) < edge_rank_map->at(e1)){
                                    continue;
                                }

                                edge_mutex_map->at(e1)->lock();
                                ++WS->at(e1);
                                edge_mutex_map->at(e1)->unlock();

                                if(edge_wing_map->at(e2) <= 1){
                                    edge_mutex_map->at(e2)->lock();
                                    ++WS->at(e2);
                                    edge_mutex_map->at(e2)->unlock();
                                }

                                if(edge_wing_map->at(e3) <= 1){
                                    edge_mutex_map->at(e3)->lock();
                                    ++WS->at(e3);
                                    edge_mutex_map->at(e3)->unlock();
                                }

                                if(edge_wing_map->at(e4) <= 1){
                                    edge_mutex_map->at(e4)->lock();
                                    ++WS->at(e4);
                                    edge_mutex_map->at(e4)->unlock();
                                }

                                sub_affected_set->insert(e1);
                                sub_affected_set->insert(e2);
                                sub_affected_set->insert(e3);
                                sub_affected_set->insert(e4);
                            }
                        }
                    }

                    global_mutex->lock();
                    affected_edge_set->merge(*sub_affected_set);
                    global_mutex->unlock();
                });
            }
            pool->barrier();
            for(const auto&e:*affected_edge_set){
                edge_rank_map->at(e) = UINT32_MAX;
                if (edge_wing_map->at(e) == 0) {
                    edge_wing_map->at(e) = 1;
                }
            }
        }




        /**
         * @brief update the remainder k-wings
         */
        uint32_t k = 2;

        while(!affected_edge_set->empty()){
            if (k <= *previous_k_max) {
                auto candidate_edge_set = make_shared<unordered_set<shared_ptr<abstract_bipartite_edge>>>();

                candidate_graph_finding2(G, edge_mutex_map, edge_rank_map, edge_wing_map, affected_edge_set,
                                         candidate_edge_set, candidate_edge_support_map, WS, k, pool);
                for (const auto &e: *candidate_edge_set) {
                    edge_wing_map->at(e) = k;
                }

                update_edge_wing_support(G, edge_mutex_map, edge_rank_map, edge_wing_map, WS, candidate_edge_set, k,
                                         pool);


                for(const auto &e:*candidate_edge_set){
                    affected_edge_set->insert(e);
                    WS->at(e) = candidate_edge_support_map->at(e);
                    candidate_edge_support_map->at(e) = 0;
                }
                /**
                  * @brief continue the loop
                  */
                ++k;
            } else {
                auto candidate_edge_set = make_shared<unordered_set<shared_ptr<abstract_bipartite_edge>>>();
                candidate_graph_finding2(G, edge_mutex_map, edge_rank_map, edge_wing_map, affected_edge_set, candidate_edge_set, candidate_edge_support_map, WS, k, pool);
                *previous_k_max = partial_wing_decomposition2(G, edge_mutex_map, edge_rank_map, edge_wing_map, candidate_edge_set, candidate_edge_support_map, WS, k, pool);
                break;
            }
        }
    }

    /**
   * @details a basic insert method
   * @param G
   * @param edge_set
   * @param edge_wing_map
   * @param previous_k_max
   * @param thread_number
   */
    void quasi_wing_maintenance::insert3(const shared_ptr<abstract_bipartite_graph> &G,
                                         const shared_ptr<unordered_map<shared_ptr<abstract_bipartite_edge>,shared_ptr<mutex>>>& edge_mutex_map,
                                         const shared_ptr<unordered_set<shared_ptr<abstract_bipartite_edge>>> &edge_set,
                                         const shared_ptr<unordered_map<shared_ptr<abstract_bipartite_edge>, uint32_t>>& edge_rank_map,
                                         const shared_ptr<unordered_map<shared_ptr<abstract_bipartite_edge>, uint32_t>>& candidate_edge_support_map,
                                         const shared_ptr<unordered_map<shared_ptr<abstract_bipartite_edge>, uint32_t>> &edge_wing_map,
                                         const shared_ptr<unordered_map<shared_ptr<abstract_bipartite_edge>, uint32_t>>& WS,
                                         const shared_ptr<uint32_t>& previous_k_max,
                                         const shared_ptr<thread_pool>& pool) {
        auto thread_number = pool->get_thread_number();
        auto global_mutex = make_shared<mutex>();

        G->insert_edge_collection(edge_set);
        for(const auto& e:*edge_set){
            edge_rank_map->insert({e,UINT32_MAX});
            edge_mutex_map->insert({e, make_shared<mutex>()});
            candidate_edge_support_map->insert({e,0});
            edge_wing_map->insert({e,0});
            WS->insert({e,0});
        }
        auto bloom_index_mutex = make_shared<mutex>();
        auto bloom_index = index_construction(G, edge_mutex_map, pool);

        auto affected_edge_set = make_shared<unordered_set<shared_ptr<abstract_bipartite_edge>>>();
        {
            auto rank_id = make_shared<uint32_t>(0);
            set_edge_rank(edge_set, edge_rank_map, rank_id);
            auto location_vector = pool->split_task(edge_set);
            for(uint32_t i = 0; i < thread_number; ++i){
                pool->submit_task([=]{
                    auto sub_affected_set = make_shared<unordered_set<shared_ptr<abstract_bipartite_edge>>>();

                    for(auto iter = *location_vector->at(i); iter != *location_vector->at(i + 1); ++iter){
                        const auto &e1 = *iter;
                        auto l1 = e1->get_left_vertex_id();
                        auto r1 = e1->get_right_vertex_id();
                        auto bloom_set = bloom_index->get_bloom_set(e1);
                        if(bloom_set){

                            for (const auto &B: *bloom_set) {
                                if (B->get_butterfly_count() > 0) {
                                    auto e2 = B->get_twin(e1);

                                    if (edge_rank_map->at(e2) < edge_rank_map->at(e1) ) {
                                        continue;
                                    }

                                    auto visited_set = make_shared<unordered_set<shared_ptr<abstract_bipartite_edge>>>();
                                    for (const auto &[e3, e4]: *B->get_edge_map()) {
                                        if (e3 == e1 || e3 == e2 || visited_set->count(e4)) {
                                            continue;
                                        }
                                        visited_set->insert(e3);
//                                        if ((e3->get_left_vertex_id() == e4->get_left_vertex_id() &&
//                                             e3->get_right_vertex_id() < e4->get_right_vertex_id())
//                                            || (e3->get_right_vertex_id() == e4->get_right_vertex_id() &&
//                                                e3->get_left_vertex_id() < e4->get_left_vertex_id())) {
//                                            continue;
//                                        }

                                        if ( edge_rank_map->at(e3) < edge_rank_map->at(e1)
                                             || edge_rank_map->at(e4) < edge_rank_map->at(e1) ) {
                                            continue;
                                        }

//                                        if(e3->get_left_vertex_id()!=l1 && e3->get_right_vertex_id()!=r1){
//                                            if(visited_set->count(e3)){
//                                                continue;
//                                            }
//                                            visited_set->insert(e3);
//                                        }
//                                        else{
//                                            if(visited_set->count(e4)){
//                                                continue;
//                                            }
//                                            visited_set->insert(e4);
//                                        }

                                        edge_mutex_map->at(e1)->lock();
                                        ++WS->at(e1);
                                        edge_mutex_map->at(e1)->unlock();

                                        if (edge_wing_map->at(e2) <= 1) {
                                            edge_mutex_map->at(e2)->lock();
                                            ++WS->at(e2);
                                            edge_mutex_map->at(e2)->unlock();
                                        }

                                        if (edge_wing_map->at(e3) <= 1) {
                                            edge_mutex_map->at(e3)->lock();
                                            ++WS->at(e3);
                                            edge_mutex_map->at(e3)->unlock();
                                        }

                                        if (edge_wing_map->at(e4) <= 1) {
                                            edge_mutex_map->at(e4)->lock();
                                            ++WS->at(e4);
                                            edge_mutex_map->at(e4)->unlock();
                                        }

                                        sub_affected_set->insert(e1);
                                        sub_affected_set->insert(e2);
                                        sub_affected_set->insert(e3);
                                        sub_affected_set->insert(e4);
                                    }
                                }
                            }
                        }
                        global_mutex->lock();
                        affected_edge_set->merge(*sub_affected_set);
                        global_mutex->unlock();
                    }

                });
            }
            pool->barrier();
            for(const auto&e:*affected_edge_set){
                edge_rank_map->at(e) = UINT32_MAX;
                if (edge_wing_map->at(e) == 0) {
                    edge_wing_map->at(e) = 1;
                }
            }
        }

        /**
         * @brief update the remainder k-wings
         */
        uint32_t k = 2;
        while(!affected_edge_set->empty()){
            if (k <= *previous_k_max) {
                auto candidate_edge_set = make_shared<unordered_set<shared_ptr<abstract_bipartite_edge>>>();
                simple_timer t;
                candidate_graph_finding3(G, bloom_index, edge_mutex_map, edge_rank_map, edge_wing_map, affected_edge_set, candidate_edge_set, candidate_edge_support_map, WS, k, pool);
                for(const auto &e:*candidate_edge_set){
                    edge_wing_map->at(e) = k;
                    WS->at(e) = candidate_edge_support_map->at(e);
                    candidate_edge_support_map->at(e) = 0;
                }
                update_edge_wing_support(G, bloom_index, edge_mutex_map, edge_rank_map, edge_wing_map, WS,
                                         candidate_edge_set, k, pool);
                affected_edge_set->merge(*candidate_edge_set);
                /**
                  * @brief continue the loop
                  */
                ++k;
            } else {
                auto candidate_edge_set = make_shared<unordered_set<shared_ptr<abstract_bipartite_edge>>>();
                candidate_graph_finding3(G, bloom_index, edge_mutex_map, edge_rank_map, edge_wing_map, affected_edge_set, candidate_edge_set, candidate_edge_support_map, WS, k, pool);
                *previous_k_max = partial_wing_decomposition3(G, bloom_index, edge_mutex_map, edge_rank_map, edge_wing_map, candidate_edge_set, candidate_edge_support_map, WS, k, pool);
                break;
            }
        }
    }

    void quasi_wing_maintenance::insert4(const shared_ptr<abstract_bipartite_graph> &G,
                                         const shared_ptr<unordered_map<shared_ptr<abstract_bipartite_edge>,shared_ptr<mutex>>>& edge_mutex_map,
                                         const shared_ptr<unordered_set<shared_ptr<abstract_bipartite_edge>>> &edge_set,
                                         const shared_ptr<unordered_map<shared_ptr<abstract_bipartite_edge>, uint32_t>>& edge_rank_map,
                                         const shared_ptr<unordered_map<shared_ptr<abstract_bipartite_edge>, uint32_t>>& candidate_edge_support_map,
                                         const shared_ptr<unordered_map<shared_ptr<abstract_bipartite_edge>, uint32_t>> &edge_wing_map,
                                         const shared_ptr<unordered_map<shared_ptr<abstract_bipartite_edge>, uint32_t>>& WS,
                                         const shared_ptr<unordered_map<uint32_t,shared_ptr<map<uint32_t,shared_ptr<unordered_set<uint32_t>>>>>> &WL,
                                         const shared_ptr<uint32_t>& previous_k_max,
                                         const shared_ptr<thread_pool>& pool) {
        auto thread_number = pool->get_thread_number();
        auto global_mutex = make_shared<mutex>();

        G->insert_edge_collection(edge_set);
        for(const auto& e:*edge_set){
            auto l = e->get_left_vertex_id();
            auto r = e->get_right_vertex_id();
            edge_rank_map->insert({e,UINT32_MAX});
            edge_mutex_map->insert({e, make_shared<mutex>()});
            candidate_edge_support_map->insert({e,0});
            edge_wing_map->insert({e,0});
            WS->insert({e,0});

            if(!WL->count(l)){
                WL->insert({l, make_shared<map<uint32_t ,shared_ptr<unordered_set<uint32_t>>>>()});
            }
            if(!WL->count(r)){
                WL->insert({r, make_shared<map<uint32_t ,shared_ptr<unordered_set<uint32_t>>>>()});
            }
            if(!WL->at(l)->count(0)){
                WL->at(l)->insert({0, make_shared<unordered_set<uint32_t>>()});
            }
            if(!WL->at(r)->count(0)){
                WL->at(r)->insert({0, make_shared<unordered_set<uint32_t>>()});
            }
            WL->at(l)->at(0)->insert(r);
            WL->at(r)->at(0)->insert(l);
        }
        auto affected_edge_set = make_shared<unordered_set<shared_ptr<abstract_bipartite_edge>>>();
        {
            auto rank_id = make_shared<uint32_t>(0);
            set_edge_rank(edge_set, edge_rank_map, rank_id);
            auto location_vector = pool->split_task(edge_set);
            for(uint32_t i = 0; i < thread_number; ++i){
                pool->submit_task([=]{
                    auto sub_affected_set = make_shared<unordered_set<shared_ptr<abstract_bipartite_edge>>>();

                    for(auto iter = *location_vector->at(i); iter != *location_vector->at(i + 1); ++iter){
                        const auto &e1 = *iter;
                        auto l1 = e1->get_left_vertex_id();
                        auto r1 = e1->get_right_vertex_id();

                        for(const auto &[r2,e2]:*G->get_left_vertex(l1)->get_edge_map()){
                            if(r2 == r1 || edge_rank_map->at(e2) < edge_rank_map->at(e1)){
                                continue;
                            }
                            for(const auto&[l2,e3]:*G->get_right_vertex(r1)->get_edge_map()){
                                if(l2 == l1 || edge_rank_map->at(e3) < edge_rank_map->at(e1)){
                                    continue;
                                }

                                auto e4 = G->get_edge(l2,r2);
                                if(!e4 || edge_rank_map->at(e4) < edge_rank_map->at(e1)){
                                    continue;
                                }

                                edge_mutex_map->at(e1)->lock();
                                ++WS->at(e1);
                                edge_mutex_map->at(e1)->unlock();

                                if(edge_wing_map->at(e2) <= 1){
                                    edge_mutex_map->at(e2)->lock();
                                    ++WS->at(e2);
                                    edge_mutex_map->at(e2)->unlock();
                                }

                                if(edge_wing_map->at(e3) <= 1){
                                    edge_mutex_map->at(e3)->lock();
                                    ++WS->at(e3);
                                    edge_mutex_map->at(e3)->unlock();
                                }

                                if(edge_wing_map->at(e4) <= 1){
                                    edge_mutex_map->at(e4)->lock();
                                    ++WS->at(e4);
                                    edge_mutex_map->at(e4)->unlock();
                                }

                                sub_affected_set->insert(e1);
                                sub_affected_set->insert(e2);
                                sub_affected_set->insert(e3);
                                sub_affected_set->insert(e4);
                            }
                        }
                    }

                    global_mutex->lock();
                    affected_edge_set->merge(*sub_affected_set);
                    global_mutex->unlock();
                });
            }
            pool->barrier();
            for(const auto&e:*affected_edge_set){
                edge_rank_map->at(e) = UINT32_MAX;
                if (edge_wing_map->at(e) == 0) {
                    edge_wing_map->at(e) = 1;
                    wl_move(e,WL,0,1);
                }
            }
        }

        /**
         * @brief update the remainder k-wings
         */
        uint32_t k = 2;

        while(!affected_edge_set->empty()){
            if (k <= *previous_k_max) {
                auto candidate_edge_set = make_shared<unordered_set<shared_ptr<abstract_bipartite_edge>>>();

                candidate_graph_finding4(G, edge_mutex_map, edge_rank_map, edge_wing_map, affected_edge_set, candidate_edge_set, candidate_edge_support_map, WS, WL, k, pool);
                for(const auto &e:*candidate_edge_set){
                    edge_wing_map->at(e) = k;
                    WS->at(e) = candidate_edge_support_map->at(e);
                    candidate_edge_support_map->at(e) = 0;
                    wl_move(e,WL,k-1,k);
                }
                update_edge_wing_support(G, edge_mutex_map, edge_rank_map, edge_wing_map, WS, WL, candidate_edge_set, k,
                                         pool);

                affected_edge_set->merge(*candidate_edge_set);
                /**
                  * @brief continue the loop
                  */
                ++k;
            } else {
                auto candidate_edge_set = make_shared<unordered_set<shared_ptr<abstract_bipartite_edge>>>();
                candidate_graph_finding4(G, edge_mutex_map, edge_rank_map, edge_wing_map, affected_edge_set, candidate_edge_set, candidate_edge_support_map, WS, WL, k, pool);
                *previous_k_max = partial_wing_decomposition4(G, edge_mutex_map, edge_rank_map, edge_wing_map, candidate_edge_set, candidate_edge_support_map, WS, WL, k, pool);
                break;
            }
        }
    }

    shared_ptr<BE_index> quasi_wing_maintenance::index_construction(const shared_ptr<abstract_bipartite_graph> &G,
                                                                    const shared_ptr<unordered_map<shared_ptr<abstract_bipartite_edge>, shared_ptr<mutex>>> &edge_mutex_map,
                                                                    const shared_ptr<thread_pool>& pool) {
        auto edge_set = G->get_edge_set();

        auto bloom_index = make_shared<BE_index>();
        auto vertex_priority_map = make_shared<unordered_map<uint32_t, uint32_t>>();

        vertex_priority_computation(G,vertex_priority_map, pool);

        auto edge_support_map = make_shared<unordered_map<shared_ptr<abstract_bipartite_edge>, uint32_t>>();
        for(const auto& e:*edge_set){
            edge_support_map->insert({e,0});
        }

        edge_support_computation(G, edge_set, edge_mutex_map,  edge_support_map, vertex_priority_map, pool);

        for(auto iter = edge_set->begin();iter!= edge_set->end();){
            auto &e = *iter;
            ++iter;
            auto support = edge_support_map->at(e);
            if(support > 0){
                bloom_index->insert_edge(e,support);
            }else
            {
                edge_set->erase(e);
                edge_support_map->erase(e);
            }
        }

        left_index_construction(G, vertex_priority_map, edge_support_map, edge_mutex_map, bloom_index, pool);
        right_index_construction(G, vertex_priority_map, edge_support_map, edge_mutex_map, bloom_index, pool);

        return bloom_index;
    }

    void quasi_wing_maintenance::left_index_construction(const shared_ptr<abstract_bipartite_graph> &G,
                                                         const shared_ptr<unordered_map<uint32_t, uint32_t>> &vertex_priority,
                                                         const shared_ptr<unordered_map<shared_ptr<abstract_bipartite_edge>, uint32_t>> &edge_support_map,
                                                         const shared_ptr<unordered_map<shared_ptr<abstract_bipartite_edge>, shared_ptr<mutex>>> &edge_mutex_map,
                                                         const shared_ptr<BE_index> &bloom_index,
                                                         const shared_ptr<thread_pool> &pool) {

        auto thread_number = pool->get_thread_number();
        auto bloom_index_mutex = make_shared<mutex>();

        auto left_vertex_map = G->get_left_vertex_map();
        auto location_vector = pool->split_task(left_vertex_map);
        for (uint32_t i = 0; i < thread_number; ++i) {
            pool->submit_task([=] {
                auto sub_bloom_map = make_shared<unordered_map<pair<uint32_t, uint32_t>,
                        shared_ptr<priority_obeyed_bloom>, hash_pair, equal_pair>>();

                for (auto iter = *location_vector->at(i); iter != *location_vector->at(i + 1); ++iter) {
                    auto &[l1, l1_vertex] = *iter;

                    auto wedge_count_map = make_shared<unordered_map<uint32_t, uint32_t>>();
                    wedge_count_map->reserve(l1_vertex->get_edge_map()->size());

                    for (const auto &[r1, l1r1_edge]: *l1_vertex->get_edge_map()) {
                        if (vertex_priority->at(r1) < vertex_priority->at(l1)) {
                            auto r1_vertex = G->get_right_vertex(r1);
                            for (const auto &[l2, l2r1_edge]: *r1_vertex->get_edge_map()) {
                                if (vertex_priority->at(l2) < vertex_priority->at(l1)) {
                                    if(!wedge_count_map->count(l2)){
                                        wedge_count_map->insert({l2, 0});
                                    }
                                    ++wedge_count_map->at(l2);
                                }
                            }
                        }
                    }

                    for (const auto &[r1, l1r1_edge]: *l1_vertex->get_edge_map()) {
                        if (vertex_priority->at(r1) < vertex_priority->at(l1)) {
                            auto r1_vertex = G->get_right_vertex(r1);
                            for (const auto &[l2, l2r1_edge]: *r1_vertex->get_edge_map()) {
                                if (vertex_priority->at(l2) < vertex_priority->at(l1)) {
                                    if (wedge_count_map->at(l2) > 1) {
                                        if (!sub_bloom_map->count({l1, l2})) {
                                            auto B = make_shared<priority_obeyed_bloom>(l1, l2);

                                            auto butterfly_count =
                                                    wedge_count_map->at(l2) * (wedge_count_map->at(l2) - 1) / 2;
                                            B->set_butterfly_count(butterfly_count);

                                            sub_bloom_map->insert({{l1,l2}, B});
                                        }
                                        auto  B = sub_bloom_map->at({l1, l2});
                                        B->link_twin(l1r1_edge, l2r1_edge);

                                        edge_mutex_map->at(l1r1_edge)->lock();
                                        bloom_index->link_bloom(l1r1_edge, B);
                                        edge_mutex_map->at(l1r1_edge)->unlock();

                                        edge_mutex_map->at(l2r1_edge)->lock();
                                        bloom_index->link_bloom(l2r1_edge, B);
                                        edge_mutex_map->at(l2r1_edge)->unlock();
                                    }
                                }
                            }
                        }
                    }
                }

                bloom_index_mutex->lock();
                bloom_index->get_bloom_map()->merge(*sub_bloom_map);
                bloom_index_mutex->unlock();
            });
        }
        pool->barrier();
    }


    uint32_t quasi_wing_maintenance::partial_wing_decomposition(const shared_ptr<abstract_bipartite_graph> &G,
                                                                const shared_ptr<unordered_map<shared_ptr<abstract_bipartite_edge>,shared_ptr<mutex>>>& edge_mutex_map,
                                                                const shared_ptr<unordered_map<shared_ptr<abstract_bipartite_edge>,uint32_t>>& edge_rank_map,
                                                                const shared_ptr<unordered_map<shared_ptr<abstract_bipartite_edge>, uint32_t>> &edge_wing_map,
                                                                const shared_ptr<unordered_set<shared_ptr<abstract_bipartite_edge>>> &candidate_edge_set,
                                                                const shared_ptr<unordered_map<shared_ptr<abstract_bipartite_edge>, uint32_t>> &candidate_edge_support_map,
                                                                uint32_t k,
                                                                const shared_ptr<thread_pool>& pool) {
        uint32_t k_max = 0;
        while(!candidate_edge_set->empty()){
            k_max = k;

            auto current_edge_set = make_shared<unordered_set<shared_ptr<abstract_bipartite_edge>>>();
            for (const auto &e:*candidate_edge_set) {
                edge_wing_map->at(e) = k;
                if (candidate_edge_support_map->at(e) < k + 1) {
                    current_edge_set->insert(e);
                }
            }
            auto evicted_edge_set = make_shared<unordered_set<shared_ptr<abstract_bipartite_edge>>>();
            remove_unsatisfied_edges(G, edge_mutex_map, edge_rank_map, edge_wing_map, current_edge_set, candidate_edge_set, candidate_edge_support_map, evicted_edge_set, k + 1, pool);
            ++k;
        }
        return k_max;
    }

    uint32_t quasi_wing_maintenance::partial_wing_decomposition2(const shared_ptr<abstract_bipartite_graph> &G,
                                                                 const shared_ptr<unordered_map<shared_ptr<abstract_bipartite_edge>,shared_ptr<mutex>>>& edge_mutex_map,
                                                                 const shared_ptr<unordered_map<shared_ptr<abstract_bipartite_edge>,uint32_t>>& edge_rank_map,
                                                                 const shared_ptr<unordered_map<shared_ptr<abstract_bipartite_edge>, uint32_t>> &edge_wing_map,
                                                                 const shared_ptr<unordered_set<shared_ptr<abstract_bipartite_edge>>> &candidate_edge_set,
                                                                 const shared_ptr<unordered_map<shared_ptr<abstract_bipartite_edge>, uint32_t>> &candidate_edge_support_map,
                                                                 const shared_ptr<unordered_map<shared_ptr<abstract_bipartite_edge>, uint32_t>> &WS,
                                                                 uint32_t k,
                                                                 const shared_ptr<thread_pool>& pool) {
        uint32_t k_max = 0;
        while(!candidate_edge_set->empty()){
            k_max = k;

            auto current_edge_set = make_shared<unordered_set<shared_ptr<abstract_bipartite_edge>>>();
            for (const auto &e:*candidate_edge_set) {
                edge_wing_map->at(e) = k;
                WS->at(e) = candidate_edge_support_map->at(e);

                if (candidate_edge_support_map->at(e) < k + 1) {
                    current_edge_set->insert(e);
                }
            }
            auto evicted_edge_set = make_shared<unordered_set<shared_ptr<abstract_bipartite_edge>>>();
            remove_unsatisfied_edges(G, edge_mutex_map, edge_rank_map, edge_wing_map, current_edge_set, candidate_edge_set, candidate_edge_support_map, evicted_edge_set, k + 1, pool);
            ++k;
        }
        return k_max;
    }

    uint32_t quasi_wing_maintenance::partial_wing_decomposition3(const shared_ptr<abstract_bipartite_graph> &G,
                                                                 const shared_ptr<BE_index>& bloom_index,
                                                                 const shared_ptr<unordered_map<shared_ptr<abstract_bipartite_edge>,shared_ptr<mutex>>>& edge_mutex_map,
                                                                 const shared_ptr<unordered_map<shared_ptr<abstract_bipartite_edge>,uint32_t>>& edge_rank_map,
                                                                 const shared_ptr<unordered_map<shared_ptr<abstract_bipartite_edge>, uint32_t>> &edge_wing_map,
                                                                 const shared_ptr<unordered_set<shared_ptr<abstract_bipartite_edge>>> &candidate_edge_set,
                                                                 const shared_ptr<unordered_map<shared_ptr<abstract_bipartite_edge>, uint32_t>> &candidate_edge_support_map,
                                                                 const shared_ptr<unordered_map<shared_ptr<abstract_bipartite_edge>, uint32_t>> &WS,
                                                                 uint32_t k,
                                                                 const shared_ptr<thread_pool>& pool) {
        uint32_t k_max = 0;
        while(!candidate_edge_set->empty()){
            k_max = k;

            auto current_edge_set = make_shared<unordered_set<shared_ptr<abstract_bipartite_edge>>>();
            for (const auto &e:*candidate_edge_set) {
                edge_wing_map->at(e) = k;
                WS->at(e) = candidate_edge_support_map->at(e);

                if (candidate_edge_support_map->at(e) < k + 1) {
                    current_edge_set->insert(e);
                }
            }
            auto evicted_edge_set = make_shared<unordered_set<shared_ptr<abstract_bipartite_edge>>>();
            remove_unsatisfied_edges(G, bloom_index, edge_mutex_map, edge_rank_map, edge_wing_map, current_edge_set, candidate_edge_set, candidate_edge_support_map, evicted_edge_set, k + 1, pool);
            ++k;
        }
        return k_max;
    }

    uint32_t quasi_wing_maintenance::partial_wing_decomposition4(const shared_ptr<abstract_bipartite_graph> &G,
                                                                 const shared_ptr<unordered_map<shared_ptr<abstract_bipartite_edge>,shared_ptr<mutex>>>& edge_mutex_map,
                                                                 const shared_ptr<unordered_map<shared_ptr<abstract_bipartite_edge>,uint32_t>>& edge_rank_map,
                                                                 const shared_ptr<unordered_map<shared_ptr<abstract_bipartite_edge>, uint32_t>> &edge_wing_map,
                                                                 const shared_ptr<unordered_set<shared_ptr<abstract_bipartite_edge>>> &candidate_edge_set,
                                                                 const shared_ptr<unordered_map<shared_ptr<abstract_bipartite_edge>, uint32_t>> &candidate_edge_support_map,
                                                                 const shared_ptr<unordered_map<shared_ptr<abstract_bipartite_edge>, uint32_t>> &WS,
                                                                 const shared_ptr<unordered_map<uint32_t,shared_ptr<map<uint32_t,shared_ptr<unordered_set<uint32_t>>>>>> &WL,
                                                                 uint32_t k,
                                                                 const shared_ptr<thread_pool>& pool) {
        uint32_t k_max = 0;
        while(!candidate_edge_set->empty()){
            k_max = k;

            auto current_edge_set = make_shared<unordered_set<shared_ptr<abstract_bipartite_edge>>>();
            for (const auto &e:*candidate_edge_set) {
                wl_move(e,WL,k-1,k);
                edge_wing_map->at(e) = k;
                WS->at(e) = candidate_edge_support_map->at(e);

                if (candidate_edge_support_map->at(e) < k + 1) {
                    current_edge_set->insert(e);
                }
            }
            auto evicted_edge_set = make_shared<unordered_set<shared_ptr<abstract_bipartite_edge>>>();
            remove_unsatisfied_edges(G, edge_mutex_map, edge_rank_map, edge_wing_map, current_edge_set, candidate_edge_set, candidate_edge_support_map, evicted_edge_set, WL, k + 1, pool);
            ++k;
        }
        return k_max;
    }

    /**
     * @details a top-down maintenance method
     * @param G
     * @param edge_set
     * @param edge_wing_map
     * @param previous_k_max
     * @param thread_number
     */
    void quasi_wing_maintenance::remove(const shared_ptr<abstract_bipartite_graph> &G,
                                        const shared_ptr<unordered_map<shared_ptr<abstract_bipartite_edge>,shared_ptr<mutex>>>& edge_mutex_map,
                                        const shared_ptr<unordered_set<shared_ptr<abstract_bipartite_edge>>> &edge_set,
                                        const shared_ptr<unordered_map<shared_ptr<abstract_bipartite_edge>, uint32_t>>& edge_rank_map,
                                        const shared_ptr<unordered_map<shared_ptr<abstract_bipartite_edge>, uint32_t>> &edge_wing_map,
                                        const shared_ptr<unordered_map<shared_ptr<abstract_bipartite_edge>, uint32_t>>& WS,
                                        const shared_ptr<uint32_t>& previous_k_max,
                                        const shared_ptr<thread_pool>& pool) {
        auto thread_number = pool->get_thread_number();
        auto global_mutex = make_shared<mutex>();

        auto affected_edge_set = make_shared<unordered_set<shared_ptr<abstract_bipartite_edge>>>();
        auto max_k = make_shared<uint32_t>(1);
        {
            auto rank_id = make_shared<uint32_t>(0);
            set_edge_rank(edge_set, edge_rank_map, rank_id);
            auto location_vector = pool->split_task(edge_set);
            for(uint32_t i = 0; i < thread_number; ++i){
                pool->submit_task([=]{
                    auto sub_affected_set = make_shared<unordered_set<shared_ptr<abstract_bipartite_edge>>>();
                    uint32_t sub_max_k = 1;

                    auto &sub_begin = *location_vector->at(i);
                    auto &sub_end = *location_vector->at(i + 1);

                    for(auto iter = sub_begin; iter!=sub_end; ++iter){
                        auto &e1 = *iter;
                        auto l1 = e1->get_left_vertex_id();
                        auto r1 = e1->get_right_vertex_id();

                        for (const auto &[r2, e2]: *G->get_left_vertex(l1)->get_edge_map()) {
                            if (r2 == r1 || edge_rank_map->at(e2) < edge_rank_map->at(e1)) {
                                continue;
                            }
                            for (const auto &[l2, e3]: *G->get_right_vertex(r1)->get_edge_map()) {
                                if (l2 == l1 || edge_rank_map->at(e3) < edge_rank_map->at(e1)) {
                                    continue;
                                }

                                auto e4 = G->get_edge(l2, r2);
                                if (!e4 || edge_rank_map->at(e4) < edge_rank_map->at(e1)) {
                                    continue;
                                }

                                if (!edge_set->count(e2)) {
                                    sub_affected_set->insert(e2);
                                    if (edge_wing_map->at(e2) > sub_max_k) {
                                        sub_max_k = edge_wing_map->at(e2);
                                    }
                                }

                                if (!edge_set->count(e3)) {
                                    sub_affected_set->insert(e3);
                                    if (edge_wing_map->at(e3) > sub_max_k) {
                                        sub_max_k = edge_wing_map->at(e3);
                                    }
                                }

                                if (!edge_set->count(e4)) {
                                    sub_affected_set->insert(e4);

                                    if (edge_wing_map->at(e4) > sub_max_k) {
                                        sub_max_k = edge_wing_map->at(e4);
                                    }
                                }
                            }
                        }
                    }

                    global_mutex->lock();
                    affected_edge_set->merge(*sub_affected_set);
                    if(sub_max_k > * max_k){
                        *max_k = sub_max_k;
                    }
                    global_mutex->unlock();
                });
            }
            pool->barrier();

            for(const auto&e:*edge_set){
                edge_rank_map->erase(e);
                G->remove_edge(e);
                edge_wing_map->erase(e);
            }
        }

        /**
         * @brief indicate k_max_flag is updated or not
         */
        auto update_flag = false;
        if(*max_k==*previous_k_max){
            update_flag = true;
        }

        auto k = static_cast<int64_t>(*max_k);
        while(k >= 1) {
            auto current_removed_edge_set = make_shared<unordered_set<shared_ptr<abstract_bipartite_edge>>>();
            update_single_wing(G, edge_mutex_map, edge_rank_map, edge_wing_map, affected_edge_set, WS,
                               current_removed_edge_set, k, pool);
            for (const auto &e: *current_removed_edge_set) {
                edge_wing_map->at(e) = k - 1;
                affected_edge_set->insert(e);
            }
            --k;
        }

        if (update_flag) {
            auto current_k_max = make_shared<uint32_t>(0);
            auto location_vector = pool->split_task(edge_wing_map);
            for(uint32_t i = 0; i < thread_number; ++i){
                pool->submit_task([=]{
                    auto &sub_begin = *location_vector->at(i);
                    auto &sub_end = *location_vector->at(i + 1);

                    uint32_t sub_k_max = 0;

                    for(auto iter = sub_begin; iter!=sub_end; ++iter){
                        auto &[e, e_wing_number] = *iter;
                        if (e_wing_number > sub_k_max) {
                            sub_k_max = e_wing_number;
                        }
                    }

                    global_mutex->lock();
                    if(sub_k_max > *current_k_max){
                        *current_k_max = sub_k_max;
                    }
                    global_mutex->unlock();
                });
            }
            pool->barrier();
            *previous_k_max = *current_k_max;
        }
    }

    void quasi_wing_maintenance::remove2(const shared_ptr<abstract_bipartite_graph> &G,
                                         const shared_ptr<unordered_map<shared_ptr<abstract_bipartite_edge>, shared_ptr<mutex>>> &edge_mutex_map,
                                         const shared_ptr<unordered_set<shared_ptr<abstract_bipartite_edge>>> &edge_set,
                                         const shared_ptr<unordered_map<shared_ptr<abstract_bipartite_edge>, uint32_t>>& edge_rank_map,
                                         const shared_ptr<unordered_map<shared_ptr<abstract_bipartite_edge>, uint32_t>> &edge_wing_map,
                                         const shared_ptr<unordered_map<shared_ptr<abstract_bipartite_edge>, uint32_t>> &WS,
                                         const shared_ptr<uint32_t> &previous_k_max,
                                         const shared_ptr<thread_pool> &pool) {
        auto thread_number = pool->get_thread_number();
        auto global_mutex = make_shared<mutex>();

        auto affected_edge_set = make_shared<unordered_set<shared_ptr<abstract_bipartite_edge>>>();
        auto max_k = make_shared<uint32_t>(1);
        {
            auto rank_id = make_shared<uint32_t>(0);
            set_edge_rank(edge_set, edge_rank_map, rank_id);
            auto location_vector = pool->split_task(edge_set);
            for(uint32_t i = 0; i < thread_number; ++i){
                pool->submit_task([=]{
                    auto sub_affected_set = make_shared<unordered_set<shared_ptr<abstract_bipartite_edge>>>();
                    uint32_t sub_max_k = 1;

                    auto &sub_begin = *location_vector->at(i);
                    auto &sub_end = *location_vector->at(i + 1);

                    for(auto iter = sub_begin; iter!=sub_end; ++iter){
                        auto &e1 = *iter;
                        auto l1 = e1->get_left_vertex_id();
                        auto r1 = e1->get_right_vertex_id();

                        for (const auto &[r2, e2]: *G->get_left_vertex(l1)->get_edge_map()) {
                            if (r2 == r1 || edge_rank_map->at(e2) < edge_rank_map->at(e1)) {
                                continue;
                            }
                            for (const auto &[l2, e3]: *G->get_right_vertex(r1)->get_edge_map()) {
                                if (l2 == l1 || edge_rank_map->at(e3) < edge_rank_map->at(e1)) {
                                    continue;
                                }

                                auto e4 = G->get_edge(l2, r2);
                                if (!e4 || edge_rank_map->at(e4) < edge_rank_map->at(e1)) {
                                    continue;
                                }

                                auto min_wing_number =std::min(std::min(edge_wing_map->at(e1),edge_wing_map->at(e2)), std::min(edge_wing_map->at(e3),edge_wing_map->at(e4)));

                                if (!edge_set->count(e2) && edge_wing_map->at(e2) == min_wing_number) {
                                    edge_mutex_map->at(e2)->lock();
                                    --WS->at(e2);
                                    edge_mutex_map->at(e2)->unlock();


                                    if (WS->at(e2) < edge_wing_map->at(e2)) {
                                        sub_affected_set->insert(e2);
                                    }

                                    if (edge_wing_map->at(e2) > sub_max_k) {
                                        sub_max_k = edge_wing_map->at(e2);
                                    }
                                }


                                if (!edge_set->count(e3) && edge_wing_map->at(e3) == min_wing_number) {
                                    edge_mutex_map->at(e3)->lock();
                                    --WS->at(e3);
                                    edge_mutex_map->at(e3)->unlock();

                                    if(WS->at(e3) < edge_wing_map->at(e3)){
                                        sub_affected_set->insert(e3);
                                    }

                                    if (edge_wing_map->at(e3) > sub_max_k) {
                                        sub_max_k = edge_wing_map->at(e3);
                                    }
                                }

                                if (!edge_set->count(e4) && edge_wing_map->at(e4) == min_wing_number) {
                                    edge_mutex_map->at(e4)->lock();
                                    --WS->at(e4);
                                    edge_mutex_map->at(e4)->unlock();

                                    if (WS->at(e4) < edge_wing_map->at(e4)) {
                                        sub_affected_set->insert(e4);
                                    }

                                    if (edge_wing_map->at(e4) > sub_max_k) {
                                        sub_max_k = edge_wing_map->at(e4);
                                    }
                                }
                            }
                        }
                    }

                    global_mutex->lock();
                    affected_edge_set->merge(*sub_affected_set);
                    if(sub_max_k > * max_k){
                        *max_k = sub_max_k;
                    }
                    global_mutex->unlock();
                });
            }
            pool->barrier();

            for(const auto&e:*edge_set){
                edge_rank_map->erase(e);
                G->remove_edge(e);
                edge_wing_map->erase(e);
                WS->at(e) = 0;
            }
        }

        /**
         * @brief indicate k_max_flag is updated or not
         */
        auto update_flag = false;
        if(*max_k==*previous_k_max){
            update_flag = true;
        }

        auto k = static_cast<int64_t>(*max_k);

        while(k >= 1) {
            auto current_removed_edge_set = make_shared<unordered_set<shared_ptr<abstract_bipartite_edge>>>();
            update_single_wing2(G, edge_mutex_map, edge_rank_map, edge_wing_map, affected_edge_set, WS,
                                current_removed_edge_set, k, pool);
            for (const auto &e: *current_removed_edge_set) {
                edge_wing_map->at(e) = k - 1;
                affected_edge_set->insert(e);
            }
            for(const auto&e1:*current_removed_edge_set){
                WS->at(e1) = 0;
                pool->submit_task([=]{
                    auto l1 = e1->get_left_vertex_id();
                    auto r1 = e1->get_right_vertex_id();
                    for (const auto &[r2, e2]:*G->get_left_vertex(l1)->get_edge_map()) {
                        if (r2 == r1 || edge_wing_map->at(e2) < edge_wing_map->at(e1)) {
                            continue;
                        }
                        for (const auto &[l2, e3]:*G->get_right_vertex(r1)->get_edge_map()) {
                            if (l2 == l1 || edge_wing_map->at(e3) < edge_wing_map->at(e1)) {
                                continue;
                            }

                            auto e4 = G->get_edge(l2, r2);
                            if (!e4 || edge_wing_map->at(e4) < edge_wing_map->at(e1)) {
                                continue;
                            }

                            ++WS->at(e1);
                        }
                    }
                });
            }
            pool->barrier();
            current_removed_edge_set->clear();
            --k;
        }


        if (update_flag) {
            auto current_k_max = make_shared<uint32_t>(0);
            auto location_vector = pool->split_task(edge_wing_map);
            for(uint32_t i = 0; i < thread_number; ++i){
                pool->submit_task([=]{
                    auto &sub_begin = *location_vector->at(i);
                    auto &sub_end = *location_vector->at(i + 1);

                    uint32_t sub_k_max = 0;

                    for(auto iter = sub_begin; iter!=sub_end; ++iter){
                        auto &[e, e_wing_number] = *iter;
                        if (e_wing_number > sub_k_max) {
                            sub_k_max = e_wing_number;
                        }
                    }

                    global_mutex->lock();
                    if(sub_k_max > *current_k_max){
                        *current_k_max = sub_k_max;
                    }
                    global_mutex->unlock();
                });
            }
            pool->barrier();
            *previous_k_max = *current_k_max;
        }
    }

    void quasi_wing_maintenance::remove3(const shared_ptr<abstract_bipartite_graph> &G,
                                         const shared_ptr<unordered_map<shared_ptr<abstract_bipartite_edge>, shared_ptr<mutex>>> &edge_mutex_map,
                                         const shared_ptr<unordered_set<shared_ptr<abstract_bipartite_edge>>> &edge_set,
                                         const shared_ptr<unordered_map<shared_ptr<abstract_bipartite_edge>, uint32_t>> &edge_rank_map,
                                         const shared_ptr<unordered_map<shared_ptr<abstract_bipartite_edge>, uint32_t>> &edge_wing_map,
                                         const shared_ptr<unordered_map<shared_ptr<abstract_bipartite_edge>, uint32_t>> &WS,
                                         const shared_ptr<uint32_t> &previous_k_max,
                                         const shared_ptr<thread_pool> &pool) {
        auto thread_number = pool->get_thread_number();
        auto global_mutex = make_shared<mutex>();

        auto affected_edge_set = make_shared<unordered_set<shared_ptr<abstract_bipartite_edge>>>();
        auto max_k = make_shared<uint32_t>(1);

        {
            auto rank_id = make_shared<uint32_t>(0);
            set_edge_rank(edge_set, edge_rank_map, rank_id);
            auto location_vector = pool->split_task(edge_set);
            for(uint32_t i = 0; i < thread_number; ++i){
                pool->submit_task([=]{
                    auto sub_affected_set = make_shared<unordered_set<shared_ptr<abstract_bipartite_edge>>>();
                    uint32_t sub_max_k = 1;

                    for(auto iter = *location_vector->at(i); iter != *location_vector->at(i + 1); ++iter){
                        auto &e1 = *iter;
                        auto l1 = e1->get_left_vertex_id();
                        auto r1 = e1->get_right_vertex_id();

                        for (const auto &[r2, e2]: *G->get_left_vertex(l1)->get_edge_map()) {
                            if (r2 == r1 || edge_rank_map->at(e2) < edge_rank_map->at(e1)) {
                                continue;
                            }
                            for (const auto &[l2, e3]: *G->get_right_vertex(r1)->get_edge_map()) {
                                if (l2 == l1 || edge_rank_map->at(e3) < edge_rank_map->at(e1)) {
                                    continue;
                                }

                                auto e4 = G->get_edge(l2, r2);
                                if (!e4 || edge_rank_map->at(e4) < edge_rank_map->at(e1)) {
                                    continue;
                                }

                                auto min_wing_number =std::min(std::min(edge_wing_map->at(e1),edge_wing_map->at(e2)), std::min(edge_wing_map->at(e3),edge_wing_map->at(e4)));

                                if (!edge_set->count(e2) && edge_wing_map->at(e2) == min_wing_number) {
                                    edge_mutex_map->at(e2)->lock();
                                    --WS->at(e2);
                                    edge_mutex_map->at(e2)->unlock();


                                    if (WS->at(e2) < edge_wing_map->at(e2)) {
                                        sub_affected_set->insert(e2);
                                    }

                                    if (edge_wing_map->at(e2) > sub_max_k) {
                                        sub_max_k = edge_wing_map->at(e2);
                                    }
                                }


                                if (!edge_set->count(e3) && edge_wing_map->at(e3) == min_wing_number) {
                                    edge_mutex_map->at(e3)->lock();
                                    --WS->at(e3);
                                    edge_mutex_map->at(e3)->unlock();

                                    if(WS->at(e3) < edge_wing_map->at(e3)){
                                        sub_affected_set->insert(e3);
                                    }

                                    if (edge_wing_map->at(e3) > sub_max_k) {
                                        sub_max_k = edge_wing_map->at(e3);
                                    }
                                }

                                if (!edge_set->count(e4) && edge_wing_map->at(e4) == min_wing_number) {
                                    edge_mutex_map->at(e4)->lock();
                                    --WS->at(e4);
                                    edge_mutex_map->at(e4)->unlock();

                                    if (WS->at(e4) < edge_wing_map->at(e4)) {
                                        sub_affected_set->insert(e4);
                                    }

                                    if (edge_wing_map->at(e4) > sub_max_k) {
                                        sub_max_k = edge_wing_map->at(e4);
                                    }
                                }
                            }
                        }
                    }

                    global_mutex->lock();
                    affected_edge_set->merge(*sub_affected_set);
                    if(sub_max_k > * max_k){
                        *max_k = sub_max_k;
                    }
                    global_mutex->unlock();
                });
            }
            pool->barrier();

            for(const auto&e:*edge_set){
                edge_rank_map->erase(e);
                G->remove_edge(e);
                edge_wing_map->erase(e);
                WS->at(e) = 0;
            }
        }

        /**
         * @brief indicate k_max_flag is updated or not
         */
        auto update_flag = false;
        if(*max_k==*previous_k_max){
            update_flag = true;
        }

        auto bloom_index = index_construction(G, edge_mutex_map, pool);

        auto k = static_cast<int64_t>(*max_k);

        while(k >= 1) {
            auto current_removed_edge_set = make_shared<unordered_set<shared_ptr<abstract_bipartite_edge>>>();
            update_single_wing3(G, bloom_index, edge_mutex_map, edge_rank_map, edge_wing_map, affected_edge_set, WS,
                                current_removed_edge_set, k, pool);
            for (const auto &e: *current_removed_edge_set) {
                edge_wing_map->at(e) = k - 1;
                affected_edge_set->insert(e);
            }
            {
                auto location_vector = pool->split_task(current_removed_edge_set);
                for(uint32_t i = 0; i < thread_number; ++i){
                    pool->submit_task([=]{
                        for(auto iter = *location_vector->at(i); iter!=*location_vector->at(i + 1); ++iter){
                            auto &e1 = *iter;
                            WS->at(e1) = 0;
                            edge_rank_map->at(e1) = UINT32_MAX;

                            auto bloom_set = bloom_index->get_bloom_set(e1);
                            if(bloom_set) {
                                for (const auto &B: *bloom_set) {
                                    if (B->get_butterfly_count() > 0) {
                                        auto e2 = B->get_twin(e1);

                                        if (edge_wing_map->at(e2) < edge_wing_map->at(e1)) {
                                            continue;
                                        }

                                        auto visited_set = make_shared<unordered_set<shared_ptr<abstract_bipartite_edge>>>();
                                        for (const auto &[e3, e4]: *B->get_edge_map()) {
                                            if (e3 == e1 || e3 == e2 || visited_set->count(e3)) {
                                                continue;
                                            }
                                            visited_set->insert(e4);

                                            if (edge_wing_map->at(e3) < edge_wing_map->at(e1)) {
                                                continue;
                                            }

                                            if (edge_wing_map->at(e4) < edge_wing_map->at(e1)) {
                                                continue;
                                            }

                                            ++WS->at(e1);
                                        }
                                    }
                                }
                            }
                        }
                    });
                }

            }

//            for(const auto&e1:*current_removed_edge_set){
//                pool->submit_task([=]{
//                    WS->at(e1) = 0;
//                    auto l1 = e1->get_left_vertex_id();
//                    auto r1 = e1->get_right_vertex_id();
//                    for (const auto &[r2, e2]:*G->get_left_vertex(l1)->get_edge_map()) {
//                        if (r2 == r1 || edge_wing_map->at(e2) < edge_wing_map->at(e1)) {
//                            continue;
//                        }
//                        for (const auto &[l2, e3]:*G->get_right_vertex(r1)->get_edge_map()) {
//                            if (l2 == l1 || edge_wing_map->at(e3) < edge_wing_map->at(e1)) {
//                                continue;
//                            }
//
//                            auto e4 = G->get_edge(l2, r2);
//                            if (!e4 || edge_wing_map->at(e4) < edge_wing_map->at(e1)) {
//                                continue;
//                            }
//
//                            ++WS->at(e1);
//                        }
//                    }
//                });
//            }
            pool->barrier();
            current_removed_edge_set->clear();
            --k;
        }

        if (update_flag) {
            auto current_k_max = make_shared<uint32_t>(0);
            auto location_vector = pool->split_task(edge_wing_map);
            for(uint32_t i = 0; i < thread_number; ++i){
                pool->submit_task([=]{
                    auto &sub_begin = *location_vector->at(i);
                    auto &sub_end = *location_vector->at(i + 1);

                    uint32_t sub_k_max = 0;

                    for(auto iter = sub_begin; iter!=sub_end; ++iter){
                        auto &[e, e_wing_number] = *iter;
                        if (e_wing_number > sub_k_max) {
                            sub_k_max = e_wing_number;
                        }
                    }

                    global_mutex->lock();
                    if(sub_k_max > *current_k_max){
                        *current_k_max = sub_k_max;
                    }
                    global_mutex->unlock();
                });
            }
            pool->barrier();
            *previous_k_max = *current_k_max;
        }
    }

    void quasi_wing_maintenance::remove4(const shared_ptr<abstract_bipartite_graph> &G,
                                         const shared_ptr<unordered_map<shared_ptr<abstract_bipartite_edge>, shared_ptr<mutex>>> &edge_mutex_map,
                                         const shared_ptr<unordered_set<shared_ptr<abstract_bipartite_edge>>> &edge_set,
                                         const shared_ptr<unordered_map<shared_ptr<abstract_bipartite_edge>, uint32_t>>& edge_rank_map,
                                         const shared_ptr<unordered_map<shared_ptr<abstract_bipartite_edge>, uint32_t>> &edge_wing_map,
                                         const shared_ptr<unordered_map<shared_ptr<abstract_bipartite_edge>, uint32_t>> &WS,
                                         const shared_ptr<unordered_map<uint32_t,shared_ptr<map<uint32_t,shared_ptr<unordered_set<uint32_t>>>>>> &WL,
                                         const shared_ptr<uint32_t> &previous_k_max,
                                         const shared_ptr<thread_pool> &pool) {
        auto thread_number = pool->get_thread_number();
        auto global_mutex = make_shared<mutex>();

        auto affected_edge_set = make_shared<unordered_set<shared_ptr<abstract_bipartite_edge>>>();
        auto max_k = make_shared<uint32_t>(1);
        {
            auto rank_id = make_shared<uint32_t>(0);
            set_edge_rank(edge_set, edge_rank_map, rank_id);
            auto location_vector = pool->split_task(edge_set);
            for(uint32_t i = 0; i < thread_number; ++i){
                pool->submit_task([=]{
                    auto sub_affected_set = make_shared<unordered_set<shared_ptr<abstract_bipartite_edge>>>();
                    uint32_t sub_max_k = 1;

                    auto &sub_begin = *location_vector->at(i);
                    auto &sub_end = *location_vector->at(i + 1);

                    for(auto iter = sub_begin; iter!=sub_end; ++iter){
                        auto &e1 = *iter;
                        auto l1 = e1->get_left_vertex_id();
                        auto r1 = e1->get_right_vertex_id();
                        for(auto iter2 = WL->at(l1)->begin();iter2!=WL->at(l1)->lower_bound(edge_wing_map->at(e1));iter2++){
                            for (const auto &r2: *iter2->second) {
                                auto e2 = G->get_edge(l1, r2);
                                if (r2 == r1 || edge_rank_map->at(e2) < edge_rank_map->at(e1)) {
                                    continue;
                                }

                                for(auto iter3 = WL->at(r1)->begin();iter3!= WL->at(r1)->lower_bound(edge_wing_map->at(e2));iter3++){
                                    for (const auto &l2: *iter3->second){
                                        auto e3 = G->get_edge(l2, r1);
                                        if (l2 == l1 || edge_rank_map->at(e3) < edge_rank_map->at(e1)) {
                                            continue;
                                        }

                                        auto e4 = G->get_edge(l2, r2);
                                        if (!e4 || edge_rank_map->at(e4) < edge_rank_map->at(e1)) {
                                            continue;
                                        }
                                        auto min_wing_number =std::min(edge_wing_map->at(e3), edge_wing_map->at(e4));
                                        if (!edge_set->count(e3) && edge_wing_map->at(e3) == min_wing_number) {
                                            edge_mutex_map->at(e3)->lock();
                                            --WS->at(e3);
                                            edge_mutex_map->at(e3)->unlock();

                                            if(WS->at(e3) < edge_wing_map->at(e3)){
                                                sub_affected_set->insert(e3);
                                            }

                                            if (edge_wing_map->at(e3) > sub_max_k) {
                                                sub_max_k = edge_wing_map->at(e3);
                                            }

                                        }

                                        if (!edge_set->count(e4) && edge_wing_map->at(e4) == min_wing_number) {
                                            edge_mutex_map->at(e4)->lock();
                                            --WS->at(e4);
                                            edge_mutex_map->at(e4)->unlock();

                                            if (WS->at(e4) < edge_wing_map->at(e4)) {
                                                sub_affected_set->insert(e4);
                                            }

                                            if (edge_wing_map->at(e4) > sub_max_k) {
                                                sub_max_k = edge_wing_map->at(e4);
                                            }
                                        }
                                    }
                                }

                                for(auto iter3 = WL->at(r1)->lower_bound(edge_wing_map->at(e2));iter3!= WL->at(r1)->end();iter3++){
                                    for (const auto &l2: *iter3->second){
                                        auto e3 = G->get_edge(l2, r1);
                                        if (l2 == l1 || edge_rank_map->at(e3) < edge_rank_map->at(e1)) {
                                            continue;
                                        }

                                        auto e4 = G->get_edge(l2, r2);
                                        if (!e4 || edge_rank_map->at(e4) < edge_rank_map->at(e1)) {
                                            continue;
                                        }
                                        auto min_wing_number =std::min(edge_wing_map->at(e2), edge_wing_map->at(e4));
                                        if (!edge_set->count(e2) && edge_wing_map->at(e2) == min_wing_number) {
                                            edge_mutex_map->at(e2)->lock();
                                            --WS->at(e2);
                                            edge_mutex_map->at(e2)->unlock();


                                            if (WS->at(e2) < edge_wing_map->at(e2)) {
                                                sub_affected_set->insert(e2);
                                            }

                                            if (edge_wing_map->at(e2) > sub_max_k) {
                                                sub_max_k = edge_wing_map->at(e2);
                                            }
                                        }

                                        if (!edge_set->count(e3) && edge_wing_map->at(e3) == min_wing_number) {
                                            edge_mutex_map->at(e3)->lock();
                                            --WS->at(e3);
                                            edge_mutex_map->at(e3)->unlock();

                                            if(WS->at(e3) < edge_wing_map->at(e3)){
                                                sub_affected_set->insert(e3);
                                            }

                                            if (edge_wing_map->at(e3) > sub_max_k) {
                                                sub_max_k = edge_wing_map->at(e3);
                                            }

                                        }

                                        if (!edge_set->count(e4) && edge_wing_map->at(e4) == min_wing_number) {
                                            edge_mutex_map->at(e4)->lock();
                                            --WS->at(e4);
                                            edge_mutex_map->at(e4)->unlock();

                                            if (WS->at(e4) < edge_wing_map->at(e4)) {
                                                sub_affected_set->insert(e4);
                                            }

                                            if (edge_wing_map->at(e4) > sub_max_k) {
                                                sub_max_k = edge_wing_map->at(e4);
                                            }
                                        }
                                    }
                                }


                            }
                        }



                        for(auto iter2 = WL->at(l1)->lower_bound(edge_wing_map->at(e1));iter2!=WL->at(l1)->end();iter2++){
                            for (const auto &r2: *iter2->second) {
                                auto e2 = G->get_edge(l1, r2);
                                if (r2 == r1 || edge_rank_map->at(e2) < edge_rank_map->at(e1)) {
                                    continue;
                                }

                                for(auto iter3 = WL->at(r1)->begin();iter3!= WL->at(r1)->lower_bound(edge_wing_map->at(e1));iter3++) {
                                    for (const auto &l2: *iter3->second) {
                                        auto e3 = G->get_edge(l2, r1);
                                        if (l2 == l1 || edge_rank_map->at(e3) < edge_rank_map->at(e1)) {
                                            continue;
                                        }
                                        auto e4 = G->get_edge(l2, r2);
                                        if (!e4 || edge_rank_map->at(e4) < edge_rank_map->at(e1)) {
                                            continue;
                                        }
                                        auto min_wing_number =std::min(edge_wing_map->at(e3), edge_wing_map->at(e4));
                                        if (!edge_set->count(e3) && edge_wing_map->at(e3) == min_wing_number) {
                                            edge_mutex_map->at(e3)->lock();
                                            --WS->at(e3);
                                            edge_mutex_map->at(e3)->unlock();

                                            if(WS->at(e3) < edge_wing_map->at(e3)){
                                                sub_affected_set->insert(e3);
                                            }

                                            if (edge_wing_map->at(e3) > sub_max_k) {
                                                sub_max_k = edge_wing_map->at(e3);
                                            }
                                        }

                                        if (!edge_set->count(e4) && edge_wing_map->at(e4) == min_wing_number) {
                                            edge_mutex_map->at(e4)->lock();
                                            --WS->at(e4);
                                            edge_mutex_map->at(e4)->unlock();

                                            if (WS->at(e4) < edge_wing_map->at(e4)) {
                                                sub_affected_set->insert(e4);
                                            }

                                            if (edge_wing_map->at(e4) > sub_max_k) {
                                                sub_max_k = edge_wing_map->at(e4);
                                            }
                                        }
                                    }
                                }

                                for(auto iter3 = WL->at(r1)->lower_bound(edge_wing_map->at(e1));iter3!= WL->at(r1)->end();iter3++) {
                                    for (const auto &l2: *iter3->second) {
                                        auto e3 = G->get_edge(l2, r1);
                                        if (l2 == l1 || edge_rank_map->at(e3) < edge_rank_map->at(e1)) {
                                            continue;
                                        }
                                        auto e4 = G->get_edge(l2, r2);
                                        if (!e4 || edge_rank_map->at(e4) < edge_rank_map->at(e1)) {
                                            continue;
                                        }
                                        auto min_wing_number =std::min(edge_wing_map->at(e1), edge_wing_map->at(e4));

                                        if (!edge_set->count(e2) && edge_wing_map->at(e2) == min_wing_number) {
                                            edge_mutex_map->at(e2)->lock();
                                            --WS->at(e2);
                                            edge_mutex_map->at(e2)->unlock();


                                            if (WS->at(e2) < edge_wing_map->at(e2)) {
                                                sub_affected_set->insert(e2);
                                            }

                                            if (edge_wing_map->at(e2) > sub_max_k) {
                                                sub_max_k = edge_wing_map->at(e2);
                                            }

                                        }

                                        if (!edge_set->count(e3) && edge_wing_map->at(e3) == min_wing_number) {
                                            edge_mutex_map->at(e3)->lock();
                                            --WS->at(e3);
                                            edge_mutex_map->at(e3)->unlock();

                                            if(WS->at(e3) < edge_wing_map->at(e3)){
                                                sub_affected_set->insert(e3);
                                            }

                                            if (edge_wing_map->at(e3) > sub_max_k) {
                                                sub_max_k = edge_wing_map->at(e3);
                                            }
                                        }


                                        if (!edge_set->count(e4) && edge_wing_map->at(e4) == min_wing_number) {
                                            edge_mutex_map->at(e4)->lock();
                                            --WS->at(e4);
                                            edge_mutex_map->at(e4)->unlock();

                                            if (WS->at(e4) < edge_wing_map->at(e4)) {
                                                sub_affected_set->insert(e4);
                                            }

                                            if (edge_wing_map->at(e4) > sub_max_k) {
                                                sub_max_k = edge_wing_map->at(e4);
                                            }
                                        }
                                    }
                                }

                            }

                        }


                        /*for (const auto &[r2, e2]: *G->get_left_vertex(l1)->get_edge_map()) {
                            if (r2 == r1 || edge_rank_map->at(e2) < edge_rank_map->at(e1)) {
                                continue;
                            }
                            for (const auto &[l2, e3]: *G->get_right_vertex(r1)->get_edge_map()) {
                                if (l2 == l1 || edge_rank_map->at(e3) < edge_rank_map->at(e1)) {
                                    continue;
                                }

                                auto e4 = G->get_edge(l2, r2);
                                if (!e4 || edge_rank_map->at(e4) < edge_rank_map->at(e1)) {
                                    continue;
                                }

                                auto min_wing_number =std::min(std::min(edge_wing_map->at(e1),edge_wing_map->at(e2)), std::min(edge_wing_map->at(e3),edge_wing_map->at(e4)));

                                if (!edge_set->count(e2) && edge_wing_map->at(e2) == min_wing_number) {
                                    edge_mutex_map->at(e2)->lock();
                                    --WS->at(e2);
                                    edge_mutex_map->at(e2)->unlock();


                                    if (WS->at(e2) < edge_wing_map->at(e2)) {
                                        sub_affected_set->insert(e2);
                                    }

                                    if (edge_wing_map->at(e2) > sub_max_k) {
                                        sub_max_k = edge_wing_map->at(e2);
                                    }
                                }


                                if (!edge_set->count(e3) && edge_wing_map->at(e3) == min_wing_number) {
                                    edge_mutex_map->at(e3)->lock();
                                    --WS->at(e3);
                                    edge_mutex_map->at(e3)->unlock();

                                    if(WS->at(e3) < edge_wing_map->at(e3)){
                                        sub_affected_set->insert(e3);
                                    }

                                    if (edge_wing_map->at(e3) > sub_max_k) {
                                        sub_max_k = edge_wing_map->at(e3);
                                    }
                                }

                                if (!edge_set->count(e4) && edge_wing_map->at(e4) == min_wing_number) {
                                    edge_mutex_map->at(e4)->lock();
                                    --WS->at(e4);
                                    edge_mutex_map->at(e4)->unlock();

                                    if (WS->at(e4) < edge_wing_map->at(e4)) {
                                        sub_affected_set->insert(e4);
                                    }

                                    if (edge_wing_map->at(e4) > sub_max_k) {
                                        sub_max_k = edge_wing_map->at(e4);
                                    }
                                }
                            }
                        }*/
                    }

                    global_mutex->lock();
                    affected_edge_set->merge(*sub_affected_set);
                    if(sub_max_k > * max_k){
                        *max_k = sub_max_k;
                    }
                    global_mutex->unlock();
                });
            }
            pool->barrier();

            for(const auto&e:*edge_set){
                auto l = e->get_left_vertex_id();
                auto r = e->get_right_vertex_id();
                WL->at(l)->at(edge_wing_map->at(e))->erase(r);
                WL->at(r)->at(edge_wing_map->at(e))->erase(l);
                edge_rank_map->erase(e);
                G->remove_edge(e);
                edge_wing_map->erase(e);
                WS->at(e) = 0;
            }
        }

        /**
         * @brief indicate k_max_flag is updated or not
         */
        auto update_flag = false;
        if(*max_k==*previous_k_max){
            update_flag = true;
        }

        auto k = static_cast<int64_t>(*max_k);

        while(k >= 1) {
            auto current_removed_edge_set = make_shared<unordered_set<shared_ptr<abstract_bipartite_edge>>>();
            update_single_wing4(G, edge_mutex_map, edge_rank_map, edge_wing_map, affected_edge_set, WS, WL,
                                current_removed_edge_set, k, pool);
            for (const auto &e: *current_removed_edge_set) {
                wl_move(e,WL,k,k-1);
                edge_wing_map->at(e) = k - 1;
                affected_edge_set->insert(e);
            }
            for(const auto&e1:*current_removed_edge_set){
                WS->at(e1) = 0;
                pool->submit_task([=]{
                    auto l1 = e1->get_left_vertex_id();
                    auto r1 = e1->get_right_vertex_id();
                    for(auto iter2 = WL->at(l1)->lower_bound(edge_wing_map->at(e1));iter2!=WL->at(l1)->end();++iter2) {
                        for (const auto &r2: *iter2->second) {
                            auto e2 = G->get_edge(l1, r2);
                            if (r2 == r1) {
                                continue;
                            }
                            for (auto iter3 = WL->at(r2)->lower_bound(edge_wing_map->at(e1)); iter3 != WL->at(r2)->end(); ++iter3) {
                                for (const auto &l2: *iter3->second) {
                                    auto e3 = G->get_edge(l2, r2);
                                    if (l2 == l1) {
                                        continue;
                                    }
                                    auto e4 = G->get_edge(l2, r1);
                                    if (!e4 || edge_wing_map->at(e4) < edge_wing_map->at(e1)) {
                                        continue;
                                    }

                                    ++WS->at(e1);
                                }
                            }
                        }
                    }
                });
            }
            pool->barrier();
            current_removed_edge_set->clear();
            --k;
        }


        if (update_flag) {
            auto current_k_max = make_shared<uint32_t>(0);
            auto location_vector = pool->split_task(edge_wing_map);
            for(uint32_t i = 0; i < thread_number; ++i){
                pool->submit_task([=]{
                    auto &sub_begin = *location_vector->at(i);
                    auto &sub_end = *location_vector->at(i + 1);

                    uint32_t sub_k_max = 0;

                    for(auto iter = sub_begin; iter!=sub_end; ++iter){
                        auto &[e, e_wing_number] = *iter;
                        if (e_wing_number > sub_k_max) {
                            sub_k_max = e_wing_number;
                        }
                    }

                    global_mutex->lock();
                    if(sub_k_max > *current_k_max){
                        *current_k_max = sub_k_max;
                    }
                    global_mutex->unlock();
                });
            }
            pool->barrier();
            *previous_k_max = *current_k_max;
        }
    }

    void quasi_wing_maintenance::remove_unsatisfied_edges(const shared_ptr<abstract_bipartite_graph> &G,
                                                          const shared_ptr<unordered_map<shared_ptr<abstract_bipartite_edge>, shared_ptr<mutex>>> &edge_mutex_map,
                                                          const shared_ptr<unordered_map<shared_ptr<abstract_bipartite_edge>, uint32_t>> &edge_rank_map,
                                                          const shared_ptr<unordered_map<shared_ptr<abstract_bipartite_edge>, uint32_t>> &edge_wing_map,
                                                          const shared_ptr<unordered_set<shared_ptr<abstract_bipartite_edge>>> &current_edge_set,
                                                          const shared_ptr<unordered_set<shared_ptr<abstract_bipartite_edge>>> &candidate_edge_set,
                                                          const shared_ptr<unordered_map<shared_ptr<abstract_bipartite_edge>,uint32_t>>& candidate_edge_support_map,
                                                          const shared_ptr<unordered_set<shared_ptr<abstract_bipartite_edge>>> &evicted_edge_set,
                                                          uint32_t k,
                                                          const shared_ptr<thread_pool> &pool) {
        auto thread_number = pool->get_thread_number();
        auto global_mutex = make_shared<mutex>();

        auto current_evicted_set =  make_shared<unordered_set<shared_ptr<abstract_bipartite_edge>>>();

        auto rank_id = make_shared<uint32_t>(0);

        while(!current_edge_set->empty()){
            auto next_edge_set = make_shared<unordered_set<shared_ptr<abstract_bipartite_edge>>>();
            set_edge_rank(current_edge_set, edge_rank_map, rank_id);
            auto location_vector = pool->split_task(current_edge_set);
            for(uint32_t i = 0; i < thread_number; ++i){
                pool->submit_task([=]{
                    auto sub_next_edge_set = make_shared<unordered_set<shared_ptr<abstract_bipartite_edge>>>();

                    auto &sub_begin =  *location_vector->at(i);
                    auto &sub_end =  *location_vector->at(i + 1);

                    for(auto iter = sub_begin; iter!=sub_end; ++iter){
                        auto& e1 = *iter;

                        auto l1 = e1->get_left_vertex_id();
                        auto r1 = e1->get_right_vertex_id();
                        for (const auto &[r2, e2]: *G->get_left_vertex(l1)->get_edge_map()) {
                            if (r2 == r1 || evicted_edge_set->count(e2)
                                || edge_wing_map->at(e2) < k - 1 || edge_rank_map->at(e2) < edge_rank_map->at(e1)||
                                (edge_wing_map->at(e2) == k - 1 && !candidate_edge_set->count(e2))) {
                                continue;
                            }
                            for (const auto &[l2, e3]: *G->get_right_vertex(r1)->get_edge_map()) {
                                if (l2 == l1 || evicted_edge_set->count(e3)
                                    || edge_wing_map->at(e3) < k - 1 || edge_rank_map->at(e3) < edge_rank_map->at(e1)||
                                    (edge_wing_map->at(e3) == k - 1 && !candidate_edge_set->count(e3))) {
                                    continue;
                                }

                                auto e4 = G->get_edge(l2, r2);
                                if (!e4 || evicted_edge_set->count(e4)
                                    || edge_wing_map->at(e4) < k - 1 || edge_rank_map->at(e4) < edge_rank_map->at(e1)||
                                    (edge_wing_map->at(e4) == k - 1 && !candidate_edge_set->count(e4))) {
                                    continue;
                                }

                                if (edge_wing_map->at(e2) == k - 1 && candidate_edge_support_map->at(e2) >= k) {
                                    edge_mutex_map->at(e2)->lock();
                                    --candidate_edge_support_map->at(e2);
                                    edge_mutex_map->at(e2)->unlock();
                                    if (candidate_edge_support_map->at(e2) < k && !current_edge_set->count(e2)) {
                                        sub_next_edge_set->insert(e2);
                                    }
                                }

                                if (edge_wing_map->at(e3) == k - 1 && candidate_edge_support_map->at(e3) >= k) {
                                    edge_mutex_map->at(e3)->lock();
                                    --candidate_edge_support_map->at(e3);
                                    edge_mutex_map->at(e3)->unlock();
                                    if (candidate_edge_support_map->at(e3) < k && !current_edge_set->count(e3)) {
                                        sub_next_edge_set->insert(e3);
                                    }
                                }

                                if (edge_wing_map->at(e4) == k - 1 && candidate_edge_support_map->at(e4) >= k) {
                                    edge_mutex_map->at(e4)->lock();
                                    --candidate_edge_support_map->at(e4);
                                    edge_mutex_map->at(e4)->unlock();
                                    if (candidate_edge_support_map->at(e4) < k && !current_edge_set->count(e4)) {
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
            for(const auto &e:*current_edge_set) {
                candidate_edge_set->erase(e);
                candidate_edge_support_map->at(e) = 0;
            }
            current_evicted_set->merge(*current_edge_set);
            swap(*current_edge_set, *next_edge_set);
        }


        for(const auto &e:*current_evicted_set){
            edge_rank_map->at(e) = UINT32_MAX;
        }
        evicted_edge_set->merge(*current_evicted_set);
    }

    void quasi_wing_maintenance::remove_unsatisfied_edges(const shared_ptr<abstract_bipartite_graph> &G,
                                                          const shared_ptr<unordered_map<shared_ptr<abstract_bipartite_edge>, shared_ptr<mutex>>> &edge_mutex_map,
                                                          const shared_ptr<unordered_map<shared_ptr<abstract_bipartite_edge>, uint32_t>> &edge_rank_map,
                                                          const shared_ptr<unordered_map<shared_ptr<abstract_bipartite_edge>, uint32_t>> &edge_wing_map,
                                                          const shared_ptr<unordered_set<shared_ptr<abstract_bipartite_edge>>> &current_edge_set,
                                                          const shared_ptr<unordered_set<shared_ptr<abstract_bipartite_edge>>> &candidate_edge_set,
                                                          const shared_ptr<unordered_map<shared_ptr<abstract_bipartite_edge>,uint32_t>>& candidate_edge_support_map,
                                                          const shared_ptr<unordered_set<shared_ptr<abstract_bipartite_edge>>> &evicted_edge_set,
                                                          const shared_ptr<unordered_map<uint32_t,shared_ptr<map<uint32_t,shared_ptr<unordered_set<uint32_t>>>>>> &WL,
                                                          uint32_t k,
                                                          const shared_ptr<thread_pool> &pool) {
        auto thread_number = pool->get_thread_number();
        auto global_mutex = make_shared<mutex>();

        auto current_evicted_set =  make_shared<unordered_set<shared_ptr<abstract_bipartite_edge>>>();

        auto rank_id = make_shared<uint32_t>(0);

        while(!current_edge_set->empty()){
            auto next_edge_set = make_shared<unordered_set<shared_ptr<abstract_bipartite_edge>>>();
            set_edge_rank(current_edge_set, edge_rank_map, rank_id);
            auto location_vector = pool->split_task(current_edge_set);
            for(uint32_t i = 0; i < thread_number; ++i){
                pool->submit_task([=]{
                    auto sub_next_edge_set = make_shared<unordered_set<shared_ptr<abstract_bipartite_edge>>>();

                    auto &sub_begin =  *location_vector->at(i);
                    auto &sub_end =  *location_vector->at(i + 1);

                    for(auto iter = sub_begin; iter!=sub_end; ++iter){
                        auto& e1 = *iter;

                        auto l1 = e1->get_left_vertex_id();
                        auto r1 = e1->get_right_vertex_id();

                        for (const auto &r2: *WL->at(l1)->at(k - 1)) {
                            auto e2 = G->get_edge(l1, r2);
                            if (r2 == r1 || evicted_edge_set->count(e2)
                                || edge_rank_map->at(e2) < edge_rank_map->at(e1)||
                                !candidate_edge_set->count(e2)) {
                                continue;
                            }
                            for (const auto &l2: *WL->at(r1)->at(k - 1)) {
                                auto e3 = G->get_edge(l2, r1);
                                if (l2 == l1 || evicted_edge_set->count(e3)
                                    ||  edge_rank_map->at(e3) < edge_rank_map->at(e1)||
                                    !candidate_edge_set->count(e3)) {
                                    continue;
                                }
                                auto e4 = G->get_edge(l2, r2);
                                if (!e4 || evicted_edge_set->count(e4)
                                    || edge_wing_map->at(e4) < k - 1 || edge_rank_map->at(e4) < edge_rank_map->at(e1)||
                                    (edge_wing_map->at(e4) == k - 1 && !candidate_edge_set->count(e4))) {
                                    continue;
                                }

                                if (candidate_edge_support_map->at(e2) >= k) {
                                    edge_mutex_map->at(e2)->lock();
                                    --candidate_edge_support_map->at(e2);
                                    edge_mutex_map->at(e2)->unlock();
                                    if (candidate_edge_support_map->at(e2) < k && !current_edge_set->count(e2)) {
                                        sub_next_edge_set->insert(e2);
                                    }
                                }

                                if (candidate_edge_support_map->at(e3) >= k) {
                                    edge_mutex_map->at(e3)->lock();
                                    --candidate_edge_support_map->at(e3);
                                    edge_mutex_map->at(e3)->unlock();
                                    if (candidate_edge_support_map->at(e3) < k && !current_edge_set->count(e3)) {
                                        sub_next_edge_set->insert(e3);
                                    }
                                }

                                if (edge_wing_map->at(e4) == k - 1 && candidate_edge_support_map->at(e4) >= k) {
                                    edge_mutex_map->at(e4)->lock();
                                    --candidate_edge_support_map->at(e4);
                                    edge_mutex_map->at(e4)->unlock();
                                    if (candidate_edge_support_map->at(e4) < k && !current_edge_set->count(e4)) {
                                        sub_next_edge_set->insert(e4);
                                    }
                                }
                            }
                            for (auto iter3 = WL->at(r1)->lower_bound(k); iter3 != WL->at(r1)->end(); ++iter3) {
                                for (const auto &l2: *iter3->second) {
                                    auto e3 = G->get_edge(l2, r1);
//                                    if (l2 == l1 || evicted_edge_set->count(e3)
//                                        ||  edge_rank_map->at(e3) < edge_rank_map->at(e1)||
//                                        !candidate_edge_set->count(e3)) {
//                                        continue;
//                                    }
                                    auto e4 = G->get_edge(l2, r2);
                                    if (!e4 || evicted_edge_set->count(e4)
                                        || edge_wing_map->at(e4) < k - 1 || edge_rank_map->at(e4) < edge_rank_map->at(e1)||
                                        (edge_wing_map->at(e4) == k - 1 && !candidate_edge_set->count(e4))) {
                                        continue;
                                    }

                                    if (candidate_edge_support_map->at(e2) >= k) {
                                        edge_mutex_map->at(e2)->lock();
                                        --candidate_edge_support_map->at(e2);
                                        edge_mutex_map->at(e2)->unlock();
                                        if (candidate_edge_support_map->at(e2) < k && !current_edge_set->count(e2)) {
                                            sub_next_edge_set->insert(e2);
                                        }
                                    }
//
//                                        if (candidate_edge_support_map->at(e3) >= k) {
//                                            edge_mutex_map->at(e3)->lock();
//                                            --candidate_edge_support_map->at(e3);
//                                            edge_mutex_map->at(e3)->unlock();
//                                            if (candidate_edge_support_map->at(e3) < k && !current_edge_set->count(e3)) {
//                                                sub_next_edge_set->insert(e3);
//                                            }
//                                        }

                                    if (edge_wing_map->at(e4) == k - 1 && candidate_edge_support_map->at(e4) >= k) {
                                        edge_mutex_map->at(e4)->lock();
                                        --candidate_edge_support_map->at(e4);
                                        edge_mutex_map->at(e4)->unlock();
                                        if (candidate_edge_support_map->at(e4) < k && !current_edge_set->count(e4)) {
                                            sub_next_edge_set->insert(e4);
                                        }
                                    }
                                }
                            }
                        }
                        for(auto iter2 = WL->at(l1)->lower_bound(k);iter2!=WL->at(l1)->end();++iter2) {
                            for (const auto &r2: *iter2->second) {
                                auto e2 = G->get_edge(l1, r2);
//                                if (r2 == r1 || evicted_edge_set->count(e2) ||
//                                    (current_edge_set->count(e2) && edge_rank_map->at(e2) < edge_rank_map->at(e1))||
//                                    (edge_wing_map->at(e2) == k - 1 && !rectangle_edge_set->count(e2))) {
//                                    continue;
//                                }
                                for (const auto &l2: *WL->at(r1)->at(k - 1)) {
                                    auto e3 = G->get_edge(l2, r1);
                                    if (l2 == l1 || evicted_edge_set->count(e3)
                                        ||  edge_rank_map->at(e3) < edge_rank_map->at(e1)||
                                        !candidate_edge_set->count(e3)) {
                                        continue;
                                    }
                                    auto e4 = G->get_edge(l2, r2);
                                    if (!e4 || evicted_edge_set->count(e4)
                                        || edge_wing_map->at(e4) < k - 1 || edge_rank_map->at(e4) < edge_rank_map->at(e1)||
                                        (edge_wing_map->at(e4) == k - 1 && !candidate_edge_set->count(e4))) {
                                        continue;
                                    }

//                                        if (candidate_edge_support_map->at(e2) >= k) {
//                                            edge_mutex_map->at(e2)->lock();
//                                            --candidate_edge_support_map->at(e2);
//                                            edge_mutex_map->at(e2)->unlock();
//                                            if (candidate_edge_support_map->at(e2) < k && !current_edge_set->count(e2)) {
//                                                sub_next_edge_set->insert(e2);
//                                            }
//                                        }

                                    if (candidate_edge_support_map->at(e3) >= k) {
                                        edge_mutex_map->at(e3)->lock();
                                        --candidate_edge_support_map->at(e3);
                                        edge_mutex_map->at(e3)->unlock();
                                        if (candidate_edge_support_map->at(e3) < k && !current_edge_set->count(e3)) {
                                            sub_next_edge_set->insert(e3);
                                        }
                                    }

                                    if (edge_wing_map->at(e4) == k - 1 && candidate_edge_support_map->at(e4) >= k) {
                                        edge_mutex_map->at(e4)->lock();
                                        --candidate_edge_support_map->at(e4);
                                        edge_mutex_map->at(e4)->unlock();
                                        if (candidate_edge_support_map->at(e4) < k && !current_edge_set->count(e4)) {
                                            sub_next_edge_set->insert(e4);
                                        }
                                    }
                                }
                                for (auto iter3 = WL->at(r1)->lower_bound(k); iter3 != WL->at(r1)->end(); ++iter3) {
                                    for (const auto &l2: *iter3->second) {
                                        auto e3 = G->get_edge(l2, r1);
//                                        if (l2 == l1 || evicted_edge_set->count(e3) ||
//                                            (current_edge_set->count(e3) && edge_rank_map->at(e3) < edge_rank_map->at(e1))||
//                                            (edge_wing_map->at(e3) == k - 1 && !rectangle_edge_set->count(e3))) {
//                                            continue;
//                                        }
                                        auto e4 = G->get_edge(l2, r2);
                                        if (!e4 || evicted_edge_set->count(e4)
                                            || edge_wing_map->at(e4) < k - 1 || edge_rank_map->at(e4) < edge_rank_map->at(e1)||
                                            (edge_wing_map->at(e4) == k - 1 && !candidate_edge_set->count(e4))) {
                                            continue;
                                        }

//                                        if (candidate_edge_support_map->at(e2) >= k) {
//                                            edge_mutex_map->at(e2)->lock();
//                                            --candidate_edge_support_map->at(e2);
//                                            edge_mutex_map->at(e2)->unlock();
//                                            if (candidate_edge_support_map->at(e2) < k && !current_edge_set->count(e2)) {
//                                                sub_next_edge_set->insert(e2);
//                                            }
//                                        }
//
//                                        if (candidate_edge_support_map->at(e3) >= k) {
//                                            edge_mutex_map->at(e3)->lock();
//                                            --candidate_edge_support_map->at(e3);
//                                            edge_mutex_map->at(e3)->unlock();
//                                            if (candidate_edge_support_map->at(e3) < k && !current_edge_set->count(e3)) {
//                                                sub_next_edge_set->insert(e3);
//                                            }
//                                        }

                                        if (edge_wing_map->at(e4) == k - 1 && candidate_edge_support_map->at(e4) >= k) {
                                            edge_mutex_map->at(e4)->lock();
                                            --candidate_edge_support_map->at(e4);
                                            edge_mutex_map->at(e4)->unlock();
                                            if (candidate_edge_support_map->at(e4) < k && !current_edge_set->count(e4)) {
                                                sub_next_edge_set->insert(e4);
                                            }
                                        }
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

            for(const auto &e:*current_edge_set) {
                candidate_edge_set->erase(e);
                candidate_edge_support_map->at(e) = 0;
            }
            current_evicted_set->merge(*current_edge_set);
            swap(*current_edge_set, *next_edge_set);
        }


        for(const auto &e:*current_evicted_set){
            edge_rank_map->at(e) = UINT32_MAX;
        }
        evicted_edge_set->merge(*current_evicted_set);
    }

    void quasi_wing_maintenance::remove_unsatisfied_edges(const shared_ptr<abstract_bipartite_graph> &G,
                                                          const shared_ptr<unordered_map<shared_ptr<abstract_bipartite_edge>, shared_ptr<mutex>>> &edge_mutex_map,
                                                          const shared_ptr<unordered_map<shared_ptr<abstract_bipartite_edge>, uint32_t>> &edge_rank_map,
                                                          const shared_ptr<unordered_map<shared_ptr<abstract_bipartite_edge>, uint32_t>> &edge_wing_map,
                                                          const shared_ptr<unordered_set<shared_ptr<abstract_bipartite_edge>>> &current_edge_set,
                                                          const shared_ptr<unordered_set<shared_ptr<abstract_bipartite_edge>>> &candidate_edge_set,
                                                          const shared_ptr<unordered_map<shared_ptr<abstract_bipartite_edge>,uint32_t>>& candidate_edge_support_map,
                                                          const shared_ptr<unordered_set<shared_ptr<abstract_bipartite_edge>>> &rectangle_edge_set,
                                                          const shared_ptr<unordered_map<shared_ptr<abstract_bipartite_edge>, shared_ptr<unordered_set<shared_ptr<abstract_bipartite_edge>>>>>& next_edge_map,
                                                          const shared_ptr<unordered_set<shared_ptr<abstract_bipartite_edge>>> &evicted_edge_set,
                                                          uint32_t k,
                                                          const shared_ptr<thread_pool> &pool) {
        auto thread_number = pool->get_thread_number();
        auto global_mutex = make_shared<mutex>();

        auto current_evicted_set =  make_shared<unordered_set<shared_ptr<abstract_bipartite_edge>>>();

        auto rank_id = make_shared<uint32_t>(0);

        while(!current_edge_set->empty()){
            auto next_edge_set = make_shared<unordered_set<shared_ptr<abstract_bipartite_edge>>>();
            set_edge_rank(current_edge_set, edge_rank_map, rank_id);
            auto location_vector = pool->split_task(current_edge_set);
            for(uint32_t i = 0; i < thread_number; ++i){
                pool->submit_task([=]{
                    auto sub_next_edge_set = make_shared<unordered_set<shared_ptr<abstract_bipartite_edge>>>();

                    for(auto iter = *location_vector->at(i); iter!=*location_vector->at(i + 1); ++iter){
                        auto& e1 = *iter;

                        auto l1 = e1->get_left_vertex_id();
                        auto r1 = e1->get_right_vertex_id();
                        for (const auto &[r2, e2]: *G->get_left_vertex(l1)->get_edge_map()) {
                            if (r2 == r1 || evicted_edge_set->count(e2) ||
                                edge_wing_map->at(e2) < k - 1 ||  ( current_edge_set->count(e2) && edge_rank_map->at(e2) < edge_rank_map->at(e1))||
                                (edge_wing_map->at(e2) == k - 1 && !rectangle_edge_set->count(e2))
                                    ) {
                                continue;
                            }
                            for (const auto &[l2, e3]: *G->get_right_vertex(r1)->get_edge_map()) {
                                if (l2 == l1 || evicted_edge_set->count(e3) ||
                                    edge_wing_map->at(e3) < k - 1 || (current_edge_set->count(e3) && edge_rank_map->at(e3) < edge_rank_map->at(e1))||
                                    (edge_wing_map->at(e3) == k - 1 && !rectangle_edge_set->count(e3))
                                        ) {
                                    continue;
                                }

                                auto e4 = G->get_edge(l2, r2);
                                if (!e4 || evicted_edge_set->count(e4) ||
                                    edge_wing_map->at(e4) < k - 1 || (current_edge_set->count(e4) && edge_rank_map->at(e4) < edge_rank_map->at(e1))||
                                    (edge_wing_map->at(e4) == k - 1 && !rectangle_edge_set->count(e4))
                                        ) {
                                    continue;
                                }

                                if (edge_wing_map->at(e2) == k - 1 && candidate_edge_support_map->at(e2) >= k) {
                                    edge_mutex_map->at(e2)->lock();
                                    --candidate_edge_support_map->at(e2);
                                    edge_mutex_map->at(e2)->unlock();
                                    if (candidate_edge_support_map->at(e2) < k && !current_edge_set->count(e2)) {
                                        sub_next_edge_set->insert(e2);
                                    }
                                }

                                if (edge_wing_map->at(e3) == k - 1 && candidate_edge_support_map->at(e3) >= k) {
                                    edge_mutex_map->at(e3)->lock();
                                    --candidate_edge_support_map->at(e3);
                                    edge_mutex_map->at(e3)->unlock();
                                    if (candidate_edge_support_map->at(e3) < k && !current_edge_set->count(e3)) {
                                        sub_next_edge_set->insert(e3);
                                    }
                                }

                                if (edge_wing_map->at(e4) == k - 1 && candidate_edge_support_map->at(e4) >= k) {
                                    edge_mutex_map->at(e4)->lock();
                                    --candidate_edge_support_map->at(e4);
                                    edge_mutex_map->at(e4)->unlock();
                                    if (candidate_edge_support_map->at(e4) < k && !current_edge_set->count(e4)) {
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
            for(const auto &e:*current_edge_set) {
                candidate_edge_set->erase(e);
                rectangle_edge_set->erase(e);
                next_edge_map->erase(e);
                candidate_edge_support_map->at(e) = 0;
            }
            current_evicted_set->merge(*current_edge_set);
            swap(*current_edge_set, *next_edge_set);
        }

        for(const auto &e:*current_evicted_set){
            edge_rank_map->at(e) = UINT32_MAX;
        }
        evicted_edge_set->merge(*current_evicted_set);
    }

    void quasi_wing_maintenance::remove_unsatisfied_edges(const shared_ptr<abstract_bipartite_graph> &G,
                                                          const shared_ptr<unordered_map<shared_ptr<abstract_bipartite_edge>, shared_ptr<mutex>>> &edge_mutex_map,
                                                          const shared_ptr<unordered_map<shared_ptr<abstract_bipartite_edge>, uint32_t>> &edge_rank_map,
                                                          const shared_ptr<unordered_map<shared_ptr<abstract_bipartite_edge>, uint32_t>> &edge_wing_map,
                                                          const shared_ptr<unordered_set<shared_ptr<abstract_bipartite_edge>>> &current_edge_set,
                                                          const shared_ptr<unordered_set<shared_ptr<abstract_bipartite_edge>>> &candidate_edge_set,
                                                          const shared_ptr<unordered_map<shared_ptr<abstract_bipartite_edge>,uint32_t>>& candidate_edge_support_map,
                                                          const shared_ptr<unordered_set<shared_ptr<abstract_bipartite_edge>>> &rectangle_edge_set,
                                                          const shared_ptr<unordered_map<shared_ptr<abstract_bipartite_edge>, shared_ptr<unordered_set<shared_ptr<abstract_bipartite_edge>>>>>& next_edge_map,
                                                          const shared_ptr<unordered_set<shared_ptr<abstract_bipartite_edge>>> &evicted_edge_set,
                                                          const shared_ptr<unordered_map<uint32_t,shared_ptr<map<uint32_t,shared_ptr<unordered_set<uint32_t>>>>>> &WL,
                                                          uint32_t k,
                                                          const shared_ptr<thread_pool> &pool) {
        auto thread_number = pool->get_thread_number();
        auto global_mutex = make_shared<mutex>();

        auto current_evicted_set =  make_shared<unordered_set<shared_ptr<abstract_bipartite_edge>>>();

        auto rank_id = make_shared<uint32_t>(0);

        while(!current_edge_set->empty()){
            auto next_edge_set = make_shared<unordered_set<shared_ptr<abstract_bipartite_edge>>>();
            set_edge_rank(current_edge_set, edge_rank_map, rank_id);
            auto location_vector = pool->split_task(current_edge_set);
            for(uint32_t i = 0; i < thread_number; ++i){
                pool->submit_task([=]{
                    auto sub_next_edge_set = make_shared<unordered_set<shared_ptr<abstract_bipartite_edge>>>();

                    for(auto iter = *location_vector->at(i); iter!=*location_vector->at(i + 1); ++iter){
                        auto& e1 = *iter;

                        auto l1 = e1->get_left_vertex_id();
                        auto r1 = e1->get_right_vertex_id();

                        for (const auto &r2: *WL->at(l1)->at(k - 1)) {
                            auto e2 = G->get_edge(l1, r2);
                            if (r2 == r1 || evicted_edge_set->count(e2) ||
                                (current_edge_set->count(e2) && edge_rank_map->at(e2) < edge_rank_map->at(e1)) ||
                                (edge_wing_map->at(e2) == k - 1 && !rectangle_edge_set->count(e2))) {
                                continue;
                            }
                            for (const auto &l2: *WL->at(r1)->at(k - 1)) {
                                auto e3 = G->get_edge(l2, r1);
                                if (l2 == l1 || evicted_edge_set->count(e3) ||
                                    (current_edge_set->count(e3) && edge_rank_map->at(e3) < edge_rank_map->at(e1)) ||
                                    (edge_wing_map->at(e3) == k - 1 && !rectangle_edge_set->count(e3))) {
                                    continue;
                                }
                                auto e4 = G->get_edge(l2, r2);
                                if (!e4 || evicted_edge_set->count(e4) ||
                                    edge_wing_map->at(e4) < k - 1 ||
                                    (current_edge_set->count(e4) && edge_rank_map->at(e4) < edge_rank_map->at(e1)) ||
                                    (edge_wing_map->at(e4) == k - 1 && !rectangle_edge_set->count(e4))
                                        ) {
                                    continue;
                                }

                                if (edge_wing_map->at(e2) == k - 1 && candidate_edge_support_map->at(e2) >= k) {
                                    edge_mutex_map->at(e2)->lock();
                                    --candidate_edge_support_map->at(e2);
                                    edge_mutex_map->at(e2)->unlock();
                                    if (candidate_edge_support_map->at(e2) < k && !current_edge_set->count(e2)) {
                                        sub_next_edge_set->insert(e2);
                                    }
                                }

                                if (edge_wing_map->at(e3) == k - 1 && candidate_edge_support_map->at(e3) >= k) {
                                    edge_mutex_map->at(e3)->lock();
                                    --candidate_edge_support_map->at(e3);
                                    edge_mutex_map->at(e3)->unlock();
                                    if (candidate_edge_support_map->at(e3) < k && !current_edge_set->count(e3)) {
                                        sub_next_edge_set->insert(e3);
                                    }
                                }

                                if (edge_wing_map->at(e4) == k - 1 && candidate_edge_support_map->at(e4) >= k) {
                                    edge_mutex_map->at(e4)->lock();
                                    --candidate_edge_support_map->at(e4);
                                    edge_mutex_map->at(e4)->unlock();
                                    if (candidate_edge_support_map->at(e4) < k && !current_edge_set->count(e4)) {
                                        sub_next_edge_set->insert(e4);
                                    }
                                }
                            }
                            for (auto iter3 = WL->at(r1)->lower_bound(k); iter3 != WL->at(r1)->end(); ++iter3) {
                                for (const auto &l2: *iter3->second) {
                                    auto e3 = G->get_edge(l2, r1);
//                                        if (l2 == l1 || evicted_edge_set->count(e3) ||
//                                            (current_edge_set->count(e3) && edge_rank_map->at(e3) < edge_rank_map->at(e1))||
//                                            (edge_wing_map->at(e3) == k - 1 && !rectangle_edge_set->count(e3))) {
//                                            continue;
//                                        }
                                    auto e4 = G->get_edge(l2, r2);
                                    if (!e4 || evicted_edge_set->count(e4) ||
                                        edge_wing_map->at(e4) < k - 1 || (current_edge_set->count(e4) &&
                                                                          edge_rank_map->at(e4) <
                                                                          edge_rank_map->at(e1)) ||
                                        (edge_wing_map->at(e4) == k - 1 && !rectangle_edge_set->count(e4))) {
                                        continue;
                                    }

                                    if (edge_wing_map->at(e2) == k - 1 && candidate_edge_support_map->at(e2) >= k) {
                                        edge_mutex_map->at(e2)->lock();
                                        --candidate_edge_support_map->at(e2);
                                        edge_mutex_map->at(e2)->unlock();
                                        if (candidate_edge_support_map->at(e2) < k && !current_edge_set->count(e2)) {
                                            sub_next_edge_set->insert(e2);
                                        }
                                    }
//
//                                        if (edge_wing_map->at(e3) == k - 1 && candidate_edge_support_map->at(e3) >= k) {
//                                            edge_mutex_map->at(e3)->lock();
//                                            --candidate_edge_support_map->at(e3);
//                                            edge_mutex_map->at(e3)->unlock();
//                                            if (candidate_edge_support_map->at(e3) < k && !current_edge_set->count(e3)) {
//                                                sub_next_edge_set->insert(e3);
//                                            }
//                                        }

                                    if (edge_wing_map->at(e4) == k - 1 && candidate_edge_support_map->at(e4) >= k) {
                                        edge_mutex_map->at(e4)->lock();
                                        --candidate_edge_support_map->at(e4);
                                        edge_mutex_map->at(e4)->unlock();
                                        if (candidate_edge_support_map->at(e4) < k && !current_edge_set->count(e4)) {
                                            sub_next_edge_set->insert(e4);
                                        }
                                    }
                                }
                            }
                        }
                        for(auto iter2 = WL->at(l1)->lower_bound(k);iter2!=WL->at(l1)->end();++iter2) {
                            for (const auto &r2: *iter2->second) {
                                auto e2 = G->get_edge(l1, r2);
//                                if (r2 == r1 || evicted_edge_set->count(e2) ||
//                                    (current_edge_set->count(e2) && edge_rank_map->at(e2) < edge_rank_map->at(e1))||
//                                    (edge_wing_map->at(e2) == k - 1 && !rectangle_edge_set->count(e2))) {
//                                    continue;
//                                }
                                for (const auto &l2: *WL->at(r1)->at(k - 1)) {
                                    auto e3 = G->get_edge(l2, r1);
                                    if (l2 == l1 || evicted_edge_set->count(e3) ||
                                        (current_edge_set->count(e3) && edge_rank_map->at(e3) < edge_rank_map->at(e1))||
                                        (edge_wing_map->at(e3) == k - 1 && !rectangle_edge_set->count(e3))) {
                                        continue;
                                    }
                                    auto e4 = G->get_edge(l2, r2);
                                    if (!e4 || evicted_edge_set->count(e4) ||
                                        edge_wing_map->at(e4) < k - 1 || (current_edge_set->count(e4) && edge_rank_map->at(e4) < edge_rank_map->at(e1))||
                                        (edge_wing_map->at(e4) == k - 1 && !rectangle_edge_set->count(e4))
                                            ) {
                                        continue;
                                    }

//                                        if (edge_wing_map->at(e2) == k - 1 && candidate_edge_support_map->at(e2) >= k) {
//                                            edge_mutex_map->at(e2)->lock();
//                                            --candidate_edge_support_map->at(e2);
//                                            edge_mutex_map->at(e2)->unlock();
//                                            if (candidate_edge_support_map->at(e2) < k && !current_edge_set->count(e2)) {
//                                                sub_next_edge_set->insert(e2);
//                                            }
//                                        }

                                    if (edge_wing_map->at(e3) == k - 1 && candidate_edge_support_map->at(e3) >= k) {
                                        edge_mutex_map->at(e3)->lock();
                                        --candidate_edge_support_map->at(e3);
                                        edge_mutex_map->at(e3)->unlock();
                                        if (candidate_edge_support_map->at(e3) < k && !current_edge_set->count(e3)) {
                                            sub_next_edge_set->insert(e3);
                                        }
                                    }

                                    if (edge_wing_map->at(e4) == k - 1 && candidate_edge_support_map->at(e4) >= k) {
                                        edge_mutex_map->at(e4)->lock();
                                        --candidate_edge_support_map->at(e4);
                                        edge_mutex_map->at(e4)->unlock();
                                        if (candidate_edge_support_map->at(e4) < k && !current_edge_set->count(e4)) {
                                            sub_next_edge_set->insert(e4);
                                        }
                                    }
                                }
                                for (auto iter3 = WL->at(r1)->lower_bound(k); iter3 != WL->at(r1)->end(); ++iter3) {
                                    for (const auto &l2: *iter3->second) {
                                        auto e3 = G->get_edge(l2, r1);
//                                        if (l2 == l1 || evicted_edge_set->count(e3) ||
//                                            (current_edge_set->count(e3) && edge_rank_map->at(e3) < edge_rank_map->at(e1))||
//                                            (edge_wing_map->at(e3) == k - 1 && !rectangle_edge_set->count(e3))) {
//                                            continue;
//                                        }
                                        auto e4 = G->get_edge(l2, r2);
                                        if (!e4 || evicted_edge_set->count(e4) ||
                                            edge_wing_map->at(e4) < k - 1 || (current_edge_set->count(e4) && edge_rank_map->at(e4) < edge_rank_map->at(e1))||
                                            (edge_wing_map->at(e4) == k - 1 && !rectangle_edge_set->count(e4))
                                                ) {
                                            continue;
                                        }

//                                        if (edge_wing_map->at(e2) == k - 1 && candidate_edge_support_map->at(e2) >= k) {
//                                            edge_mutex_map->at(e2)->lock();
//                                            --candidate_edge_support_map->at(e2);
//                                            edge_mutex_map->at(e2)->unlock();
//                                            if (candidate_edge_support_map->at(e2) < k && !current_edge_set->count(e2)) {
//                                                sub_next_edge_set->insert(e2);
//                                            }
//                                        }
//
//                                        if (edge_wing_map->at(e3) == k - 1 && candidate_edge_support_map->at(e3) >= k) {
//                                            edge_mutex_map->at(e3)->lock();
//                                            --candidate_edge_support_map->at(e3);
//                                            edge_mutex_map->at(e3)->unlock();
//                                            if (candidate_edge_support_map->at(e3) < k && !current_edge_set->count(e3)) {
//                                                sub_next_edge_set->insert(e3);
//                                            }
//                                        }

                                        if (edge_wing_map->at(e4) == k - 1 && candidate_edge_support_map->at(e4) >= k) {
                                            edge_mutex_map->at(e4)->lock();
                                            --candidate_edge_support_map->at(e4);
                                            edge_mutex_map->at(e4)->unlock();
                                            if (candidate_edge_support_map->at(e4) < k && !current_edge_set->count(e4)) {
                                                sub_next_edge_set->insert(e4);
                                            }
                                        }
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
            for(const auto &e:*current_edge_set) {
                candidate_edge_set->erase(e);
                rectangle_edge_set->erase(e);
                next_edge_map->erase(e);
                candidate_edge_support_map->at(e) = 0;
            }
            current_evicted_set->merge(*current_edge_set);
            swap(*current_edge_set, *next_edge_set);
        }

        for(const auto &e:*current_evicted_set){
            edge_rank_map->at(e) = UINT32_MAX;
        }
        evicted_edge_set->merge(*current_evicted_set);
    }

    void quasi_wing_maintenance::remove_unsatisfied_edges(const shared_ptr<abstract_bipartite_graph> &G,
                                                          const shared_ptr<BE_index>& bloom_index,
                                                          const shared_ptr<unordered_map<shared_ptr<abstract_bipartite_edge>, shared_ptr<mutex>>> &edge_mutex_map,
                                                          const shared_ptr<unordered_map<shared_ptr<abstract_bipartite_edge>, uint32_t>> &edge_rank_map,
                                                          const shared_ptr<unordered_map<shared_ptr<abstract_bipartite_edge>, uint32_t>> &edge_wing_map,
                                                          const shared_ptr<unordered_set<shared_ptr<abstract_bipartite_edge>>> &current_edge_set,
                                                          const shared_ptr<unordered_set<shared_ptr<abstract_bipartite_edge>>> &candidate_edge_set,
                                                          const shared_ptr<unordered_map<shared_ptr<abstract_bipartite_edge>,uint32_t>>& candidate_edge_support_map,
                                                          const shared_ptr<unordered_set<shared_ptr<abstract_bipartite_edge>>> &evicted_edge_set,
                                                          uint32_t k,
                                                          const shared_ptr<thread_pool> &pool) {
        auto thread_number = pool->get_thread_number();
        auto global_mutex = make_shared<mutex>();

        auto current_evicted_set =  make_shared<unordered_set<shared_ptr<abstract_bipartite_edge>>>();

        auto rank_id = make_shared<uint32_t>(0);
        while(!current_edge_set->empty()){
            auto next_edge_set = make_shared<unordered_set<shared_ptr<abstract_bipartite_edge>>>();
            set_edge_rank(current_edge_set, edge_rank_map, rank_id);
            auto location_vector = pool->split_task(current_edge_set);
            for(uint32_t i = 0; i < thread_number; ++i){
                pool->submit_task([=]{
                    auto sub_next_edge_set = make_shared<unordered_set<shared_ptr<abstract_bipartite_edge>>>();

                    auto &sub_begin =  *location_vector->at(i);
                    auto &sub_end =  *location_vector->at(i + 1);

                    for(auto iter = sub_begin; iter!=sub_end; ++iter) {
                        auto &e1 = *iter;
                        auto l1 = e1->get_left_vertex_id();
                        auto r1 = e1->get_right_vertex_id();
                        auto bloom_set = bloom_index->get_bloom_set(e1);
                        if(bloom_set){

                            for (const auto &B: *bloom_set) {
                                if (B->get_butterfly_count() > 0) {
                                    auto e2 = B->get_twin(e1);
                                    if (edge_wing_map->at(e2) < k - 1 ||
                                        edge_rank_map->at(e2) < edge_rank_map->at(e1) ||
                                        (edge_wing_map->at(e2) == k - 1 && !candidate_edge_set->count(e2))) {
                                        continue;
                                    }

                                    auto sub_visited_set = make_shared<unordered_set<shared_ptr<abstract_bipartite_edge>>>();
                                    for (const auto &[e3, e4]: *B->get_edge_map()) {
                                        if (e3 == e1 || e3 == e2 ||sub_visited_set->count(e3)) {
                                            continue;
                                        }

                                        sub_visited_set->insert(e4);

//                                        if ((e3->get_left_vertex_id() == e4->get_left_vertex_id() &&
//                                             e3->get_right_vertex_id() < e4->get_right_vertex_id())
//                                            || (e3->get_right_vertex_id() == e4->get_right_vertex_id() &&
//                                                e3->get_left_vertex_id() < e4->get_left_vertex_id())) {
//                                            continue;
//                                        }

                                        if (edge_wing_map->at(e3) < k - 1 ||
                                            edge_rank_map->at(e3) < edge_rank_map->at(e1) ||
                                            (edge_wing_map->at(e3) == k - 1 && !candidate_edge_set->count(e3))) {
                                            continue;
                                        }

                                        if (edge_wing_map->at(e4) < k - 1 ||
                                            edge_rank_map->at(e4) < edge_rank_map->at(e1) ||
                                            (edge_wing_map->at(e4) == k - 1 && !candidate_edge_set->count(e4))) {
                                            continue;
                                        }
//
//                                        if(e3->get_left_vertex_id()!=l1&&e3->get_right_vertex_id()!=r1){
//                                            if(sub_visited_set->count(e3)){
//                                                continue;
//                                            }
//                                            sub_visited_set->insert(e3);
//                                        }
//                                        else{
//                                            if(sub_visited_set->count(e4)){
//                                                continue;
//                                            }
//                                            sub_visited_set->insert(e4);
//                                        }

                                        if (edge_wing_map->at(e2) == k - 1 && candidate_edge_support_map->at(e2) >= k) {
                                            edge_mutex_map->at(e2)->lock();
                                            --candidate_edge_support_map->at(e2);
                                            edge_mutex_map->at(e2)->unlock();
                                            if (candidate_edge_support_map->at(e2) < k && !current_edge_set->count(e2)) {
                                                sub_next_edge_set->insert(e2);
                                            }
                                        }

                                        if (edge_wing_map->at(e3) == k - 1 && candidate_edge_support_map->at(e3) >= k) {
                                            edge_mutex_map->at(e3)->lock();
                                            --candidate_edge_support_map->at(e3);
                                            edge_mutex_map->at(e3)->unlock();
                                            if (candidate_edge_support_map->at(e3) < k && !current_edge_set->count(e3)) {
                                                sub_next_edge_set->insert(e3);
                                            }
                                        }

                                        if (edge_wing_map->at(e4) == k - 1 && candidate_edge_support_map->at(e4) >= k) {
                                            edge_mutex_map->at(e4)->lock();
                                            --candidate_edge_support_map->at(e4);
                                            edge_mutex_map->at(e4)->unlock();
                                            if (candidate_edge_support_map->at(e4) < k && !current_edge_set->count(e4)) {
                                                sub_next_edge_set->insert(e4);
                                            }
                                        }

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
            for(const auto &e:*current_edge_set) {
                candidate_edge_set->erase(e);
                candidate_edge_support_map->at(e) = 0;
            }
            current_evicted_set->merge(*current_edge_set);
            swap(*current_edge_set, *next_edge_set);
        }

        for(const auto &e:*current_evicted_set){
            edge_rank_map->at(e) = UINT32_MAX;
        }
        evicted_edge_set->merge(*current_evicted_set);
    }

    void quasi_wing_maintenance::remove_unsatisfied_edges(const shared_ptr<abstract_bipartite_graph> &G,
                                                          const shared_ptr<BE_index> &bloom_index,
                                                          const shared_ptr<unordered_map<shared_ptr<abstract_bipartite_edge>, shared_ptr<mutex>>> &edge_mutex_map,
                                                          const shared_ptr<unordered_map<shared_ptr<abstract_bipartite_edge>, uint32_t>> &edge_rank_map,
                                                          const shared_ptr<unordered_map<shared_ptr<abstract_bipartite_edge>, uint32_t>> &edge_wing_map,
                                                          const shared_ptr<unordered_set<shared_ptr<abstract_bipartite_edge>>> &current_edge_set,
                                                          const shared_ptr<unordered_set<shared_ptr<abstract_bipartite_edge>>> &candidate_edge_set,
                                                          const shared_ptr<unordered_map<shared_ptr<abstract_bipartite_edge>, uint32_t>> &candidate_edge_support_map,
                                                          const shared_ptr<unordered_set<shared_ptr<abstract_bipartite_edge>>> &rectangle_edge_set,
                                                          const shared_ptr<unordered_map<shared_ptr<abstract_bipartite_edge>, shared_ptr<unordered_set<shared_ptr<abstract_bipartite_edge>>>>> &next_edge_map,
                                                          const shared_ptr<unordered_set<shared_ptr<abstract_bipartite_edge>>> &evicted_edge_set,
                                                          uint32_t k, const shared_ptr<thread_pool> &pool) {
        auto thread_number = pool->get_thread_number();
        auto global_mutex = make_shared<mutex>();

        auto current_evicted_set =  make_shared<unordered_set<shared_ptr<abstract_bipartite_edge>>>();

        auto rank_id = make_shared<uint32_t>(0);
        while(!current_edge_set->empty()){
            auto next_edge_set = make_shared<unordered_set<shared_ptr<abstract_bipartite_edge>>>();
            set_edge_rank(current_edge_set, edge_rank_map, rank_id);
            auto location_vector = pool->split_task(current_edge_set);
            for(uint32_t i = 0; i < thread_number; ++i){
                pool->submit_task([=]{
                    auto sub_next_edge_set = make_shared<unordered_set<shared_ptr<abstract_bipartite_edge>>>();

                    for(auto iter =  *location_vector->at(i); iter!=*location_vector->at(i + 1); ++iter) {
                        auto &e1 = *iter;
                        auto l1 = e1->get_left_vertex_id();
                        auto r1 = e1->get_right_vertex_id();
                        auto bloom_set = bloom_index->get_bloom_set(e1);
                        if(bloom_set){

                            for (const auto &B: *bloom_set) {
                                if (B->get_butterfly_count() > 0) {
                                    auto e2 = B->get_twin(e1);
                                    if (edge_wing_map->at(e2) < k - 1 || evicted_edge_set->count(e2)
                                        || ( current_edge_set->count(e2) && edge_rank_map->at(e2) < edge_rank_map->at(e1))
                                        || (edge_wing_map->at(e2) == k - 1 && !rectangle_edge_set->count(e2))) {
                                        continue;
                                    }

                                    auto sub_visited_set = make_shared<unordered_set<shared_ptr<abstract_bipartite_edge>>>();
                                    for (const auto &[e3, e4]: *B->get_edge_map()) {
                                        if (e3 == e1 || e3 == e2 ||
                                            evicted_edge_set->count(e3) ||evicted_edge_set->count(e4) || sub_visited_set->count(e3)) {
                                            continue;
                                        }
                                        sub_visited_set->insert(e4);

//                                        if ((e3->get_left_vertex_id() == e4->get_left_vertex_id() &&
//                                             e3->get_right_vertex_id() < e4->get_right_vertex_id())
//                                            || (e3->get_right_vertex_id() == e4->get_right_vertex_id() &&
//                                                e3->get_left_vertex_id() < e4->get_left_vertex_id())) {
//                                            continue;
//                                        }

                                        if (edge_wing_map->at(e3) < k - 1
                                            || (current_edge_set->count(e3) && edge_rank_map->at(e3) < edge_rank_map->at(e1))
                                            || (edge_wing_map->at(e3) == k - 1 && !rectangle_edge_set->count(e3))) {
                                            continue;
                                        }

                                        if (edge_wing_map->at(e4) < k - 1
                                            || (current_edge_set->count(e4) && edge_rank_map->at(e4) < edge_rank_map->at(e1))
                                            || (edge_wing_map->at(e4) == k - 1 && !rectangle_edge_set->count(e4))) {
                                            continue;
                                        }

//                                        if(e3->get_left_vertex_id()!=l1 && e3->get_right_vertex_id()!=r1){
//                                            if(sub_visited_set->count(e3)){
//                                                continue;
//                                            }
//                                            sub_visited_set->insert(e3);
//                                        }
//                                        else{
//                                            if(sub_visited_set->count(e4)){
//                                                continue;
//                                            }
//                                            sub_visited_set->insert(e4);
//                                        }

                                        if (edge_wing_map->at(e2) == k - 1 && candidate_edge_support_map->at(e2) >= k) {
                                            edge_mutex_map->at(e2)->lock();
                                            --candidate_edge_support_map->at(e2);
                                            edge_mutex_map->at(e2)->unlock();
                                            if (candidate_edge_support_map->at(e2) < k && !current_edge_set->count(e2)) {
                                                sub_next_edge_set->insert(e2);
                                            }
                                        }

                                        if (edge_wing_map->at(e3) == k - 1 && candidate_edge_support_map->at(e3) >= k) {
                                            edge_mutex_map->at(e3)->lock();
                                            --candidate_edge_support_map->at(e3);
                                            edge_mutex_map->at(e3)->unlock();
                                            if (candidate_edge_support_map->at(e3) < k && !current_edge_set->count(e3)) {
                                                sub_next_edge_set->insert(e3);
                                            }
                                        }

                                        if (edge_wing_map->at(e4) == k - 1 && candidate_edge_support_map->at(e4) >= k) {
                                            edge_mutex_map->at(e4)->lock();
                                            --candidate_edge_support_map->at(e4);
                                            edge_mutex_map->at(e4)->unlock();
                                            if (candidate_edge_support_map->at(e4) < k && !current_edge_set->count(e4)) {
                                                sub_next_edge_set->insert(e4);
                                            }
                                        }

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
            for(const auto &e:*current_edge_set) {
                candidate_edge_set->erase(e);
                candidate_edge_support_map->at(e) = 0;
                rectangle_edge_set->erase(e);
                next_edge_map->erase(e);
            }
            current_evicted_set->merge(*current_edge_set);
            swap(*current_edge_set, *next_edge_set);
        }

        for(const auto &e:*current_evicted_set){
            edge_rank_map->at(e) = UINT32_MAX;
        }
        evicted_edge_set->merge(*current_evicted_set);

    }

    void quasi_wing_maintenance::right_index_construction(const shared_ptr<abstract_bipartite_graph> &G,
                                                          const shared_ptr<unordered_map<uint32_t, uint32_t>> &vertex_priority_map,
                                                          const shared_ptr<unordered_map<shared_ptr<abstract_bipartite_edge>, uint32_t>> &edge_support_map,
                                                          const shared_ptr<unordered_map<shared_ptr<abstract_bipartite_edge>, shared_ptr<mutex>>> &edge_mutex_map,
                                                          const shared_ptr<BE_index> &bloom_index,
                                                          const shared_ptr<thread_pool>& pool) {
        auto thread_number = pool->get_thread_number();
        auto global_mutex = make_shared<mutex>();

        auto right_vertex_map = G->get_right_vertex_map();
        auto location_vector = pool->split_task(right_vertex_map);
        for(uint32_t i = 0; i < thread_number; ++i){
            pool->submit_task([=]{
                auto sub_bloom_map = make_shared<unordered_map<pair<uint32_t, uint32_t>,
                        shared_ptr<priority_obeyed_bloom>, hash_pair, equal_pair>>();

                for(auto iter = *location_vector->at(i); iter!= *location_vector->at(i + 1); ++iter){
                    auto &[r1,r1_vertex] = *iter;

                    auto wedge_count_map = make_shared<unordered_map<uint32_t, uint32_t>>();
                    wedge_count_map->reserve(r1_vertex->get_edge_map()->size());

                    for(const auto&[l1,l1r1_edge]:*r1_vertex->get_edge_map()){
                        if(vertex_priority_map->at(l1) < vertex_priority_map->at(r1)){
                            auto l1_vertex = G->get_left_vertex(l1);
                            for(const auto&[r2,l1r2_edge]:*l1_vertex->get_edge_map()){
                                if(vertex_priority_map->at(r2) < vertex_priority_map->at(r1))
                                {
                                    if(!wedge_count_map->count(r2))
                                    {
                                        wedge_count_map->insert({r2,0});
                                    }
                                    ++wedge_count_map->at(r2);
                                }
                            }
                        }
                    }
                    for(const auto&[l1,e1]:*r1_vertex->get_edge_map()){
                        if(vertex_priority_map->at(l1) < vertex_priority_map->at(r1)){
                            auto l1_vertex = G->get_left_vertex(l1);
                            for(const auto&[r2,e2]:*l1_vertex->get_edge_map()){
                                if(vertex_priority_map->at(r2) < vertex_priority_map->at(r1)){
                                    if(wedge_count_map->at(r2) > 1){
                                        if(!sub_bloom_map->count({r1, r2})){
                                            auto B = make_shared<priority_obeyed_bloom>(r1,r2);
                                            sub_bloom_map->insert({{r1, r2},B});
                                        }
                                        auto B = sub_bloom_map->at({r1,r2});

                                        auto butterfly_count = wedge_count_map->at(r2)*(wedge_count_map->at(r2)-1)/2;
                                        B->set_butterfly_count(butterfly_count);
                                        B->link_twin(e1, e2);

                                        edge_mutex_map->at(e1)->lock();
                                        bloom_index->link_bloom(e1, B);
                                        edge_mutex_map->at(e1)->unlock();

                                        edge_mutex_map->at(e2)->lock();
                                        bloom_index->link_bloom(e2, B);
                                        edge_mutex_map->at(e2)->unlock();
                                    }
                                }
                            }
                        }
                    }
                }

                global_mutex->lock();
                bloom_index->get_bloom_map()->merge(*sub_bloom_map);
                global_mutex->unlock();
            });
        }
        pool->barrier();
    }

    void quasi_wing_maintenance::set_edge_rank(const shared_ptr<unordered_set<shared_ptr<scnu::abstract_bipartite_edge>>> &edge_set,
                                               const shared_ptr<unordered_map<shared_ptr<scnu::abstract_bipartite_edge>, uint32_t>> &edge_rank_map,
                                               const shared_ptr<uint32_t> &rank_id) {
        for (const auto &e: *edge_set) {
            *rank_id = *rank_id + 1;
            edge_rank_map->at(e) = *rank_id;
        }
    }

    void quasi_wing_maintenance::update_edge_wing_support(const shared_ptr<abstract_bipartite_graph> &G,
                                                          const shared_ptr<unordered_map<shared_ptr<abstract_bipartite_edge>, shared_ptr<mutex>>> &edge_mutex_map,
                                                          const shared_ptr<unordered_map<shared_ptr<abstract_bipartite_edge>, uint32_t>> &edge_rank_map,
                                                          const shared_ptr<unordered_map<shared_ptr<abstract_bipartite_edge>, uint32_t>> &edge_wing_map,
                                                          const shared_ptr<unordered_map<shared_ptr<abstract_bipartite_edge>, uint32_t>> &WS,
                                                          const shared_ptr<unordered_set<shared_ptr<abstract_bipartite_edge>>> &candidate_edge_set,
                                                          uint32_t k,
                                                          const shared_ptr<thread_pool> &pool)
    {
        auto rank_id = make_shared<uint32_t>(0);
        set_edge_rank(candidate_edge_set, edge_rank_map, rank_id);
        auto thread_number = pool->get_thread_number();
        auto global_mutex = make_shared<mutex>();
        auto location_vector = pool->split_task(candidate_edge_set);
        for(uint32_t i = 0; i < thread_number; ++i){
            pool->submit_task([=]{
                auto &sub_begin = *location_vector->at(i);
                auto &sub_end = *location_vector->at(i + 1);

                for(auto iter = sub_begin; iter!=sub_end; ++iter){
                    auto &e1 = *iter;
                    auto l1 = e1->get_left_vertex_id();
                    auto r1 = e1->get_right_vertex_id();

                    for(const auto &[r2,e2]:*G->get_left_vertex(l1)->get_edge_map()){
                        if(r2 == r1 || edge_wing_map->at(e2) < k ||edge_rank_map->at(e2) < edge_rank_map->at(e1)){
                            continue;
                        }
                        for(const auto &[l2,e3]:*G->get_right_vertex(r1)->get_edge_map()){
                            if(l2 == l1 || edge_wing_map->at(e3) < k || edge_rank_map->at(e3) < edge_rank_map->at(e1) )
                            {
                                continue;
                            }

                            auto e4 = G->get_edge(l2,r2);
                            if(!e4 || edge_wing_map->at(e4) < k || edge_rank_map->at(e4) < edge_rank_map->at(e1)){
                                continue;
                            }

                            if(edge_wing_map->at(e2) == k && !candidate_edge_set->count(e2))
                            {
                                edge_mutex_map->at(e2)->lock();
                                ++WS->at(e2);
                                edge_mutex_map->at(e2)->unlock();
                            }

                            if(edge_wing_map->at(e3) == k && !candidate_edge_set->count(e3))
                            {
                                edge_mutex_map->at(e3)->lock();
                                ++WS->at(e3);
                                edge_mutex_map->at(e3)->unlock();
                            }

                            if(edge_wing_map->at(e4) == k && !candidate_edge_set->count(e4))
                            {
                                edge_mutex_map->at(e4)->lock();
                                ++WS->at(e4);
                                edge_mutex_map->at(e4)->unlock();
                            }
                        }
                    }
                }
            });
        }
        pool->barrier();
        for (const auto &e: *candidate_edge_set) {
            edge_rank_map->at(e) = UINT32_MAX;
        }
    }

    void quasi_wing_maintenance::update_edge_wing_support(const shared_ptr<abstract_bipartite_graph> &G,
                                                          const shared_ptr<BE_index>& bloom_index,
                                                          const shared_ptr<unordered_map<shared_ptr<abstract_bipartite_edge>, shared_ptr<mutex>>> &edge_mutex_map,
                                                          const shared_ptr<unordered_map<shared_ptr<abstract_bipartite_edge>, uint32_t>> &edge_rank_map,
                                                          const shared_ptr<unordered_map<shared_ptr<abstract_bipartite_edge>, uint32_t>> &edge_wing_map,
                                                          const shared_ptr<unordered_map<shared_ptr<abstract_bipartite_edge>, uint32_t>> &WS,
                                                          const shared_ptr<unordered_set<shared_ptr<abstract_bipartite_edge>>> &candidate_edge_set,
                                                          uint32_t k,
                                                          const shared_ptr<thread_pool> &pool)
    {
        auto rank_id = make_shared<uint32_t>(0);
        set_edge_rank(candidate_edge_set, edge_rank_map, rank_id);
        auto thread_number = pool->get_thread_number();
        auto global_mutex = make_shared<mutex>();
        auto location_vector = pool->split_task(candidate_edge_set);
        for(uint32_t i = 0; i < thread_number; ++i){
            pool->submit_task([=]{
                for(auto iter = *location_vector->at(i); iter!=*location_vector->at(i + 1); ++iter) {
                    auto &e1 = *iter;
                    auto l1 = e1->get_left_vertex_id();
                    auto r1 = e1->get_right_vertex_id();
                    auto bloom_set = bloom_index->get_bloom_set(e1);
                    if(bloom_set){
                        for (const auto &B: *bloom_set) {
                            if (B->get_butterfly_count() > 0) {
                                auto e2 = B->get_twin(e1);

                                if (edge_wing_map->at(e2) < k || edge_rank_map->at(e2) < edge_rank_map->at(e1)) {
                                    continue;
                                }

                                auto visited_set = make_shared<unordered_set<shared_ptr<abstract_bipartite_edge>>>();
                                for (const auto &[e3, e4]: *B->get_edge_map()) {
                                    if (e3 == e1 || e3 == e2 || visited_set->count(e4)) {
                                        continue;
                                    }
                                    visited_set->insert(e3);
//                                    if ((e3->get_left_vertex_id() == e4->get_left_vertex_id() &&
//                                         e3->get_right_vertex_id() < e4->get_right_vertex_id())
//                                        || (e3->get_right_vertex_id() == e4->get_right_vertex_id() &&
//                                            e3->get_left_vertex_id() < e4->get_left_vertex_id())) {
//                                        continue;
//                                    }
                                    if (edge_wing_map->at(e3) < k || edge_rank_map->at(e3) < edge_rank_map->at(e1)) {
                                        continue;
                                    }

                                    if (edge_wing_map->at(e4) < k || edge_rank_map->at(e4) < edge_rank_map->at(e1)) {
                                        continue;
                                    }
//                                    if(e3->get_left_vertex_id()!=l1&&e3->get_right_vertex_id()!=r1){
//                                        if(visited_set->count(e3)){
//                                            continue;
//                                        }
//                                        visited_set->insert(e3);
//                                    }
//                                    else{
//                                        if(visited_set->count(e4)){
//                                            continue;
//                                        }
//                                        visited_set->insert(e4);
//                                    }

                                    if (edge_wing_map->at(e2) == k && !candidate_edge_set->count(e2)) {
                                        edge_mutex_map->at(e2)->lock();
                                        ++WS->at(e2);
                                        edge_mutex_map->at(e2)->unlock();
                                    }

                                    if (edge_wing_map->at(e3) == k && !candidate_edge_set->count(e3)) {
                                        edge_mutex_map->at(e3)->lock();
                                        ++WS->at(e3);
                                        edge_mutex_map->at(e3)->unlock();
                                    }

                                    if (edge_wing_map->at(e4) == k && !candidate_edge_set->count(e4)) {
                                        edge_mutex_map->at(e4)->lock();
                                        ++WS->at(e4);
                                        edge_mutex_map->at(e4)->unlock();
                                    }
                                }
                            }
                        }
                    }
                }
            });
        }
        pool->barrier();
        for(const auto&e:*candidate_edge_set){
            edge_rank_map->at(e) = UINT32_MAX;
        }
    }

    void quasi_wing_maintenance::update_edge_wing_support(const shared_ptr<abstract_bipartite_graph> &G,
                                                          const shared_ptr<unordered_map<shared_ptr<abstract_bipartite_edge>, shared_ptr<mutex>>> &edge_mutex_map,
                                                          const shared_ptr<unordered_map<shared_ptr<abstract_bipartite_edge>, uint32_t>> &edge_rank_map,
                                                          const shared_ptr<unordered_map<shared_ptr<abstract_bipartite_edge>, uint32_t>> &edge_wing_map,
                                                          const shared_ptr<unordered_map<shared_ptr<abstract_bipartite_edge>, uint32_t>> &WS,
                                                          const shared_ptr<unordered_map<uint32_t,shared_ptr<map<uint32_t,shared_ptr<unordered_set<uint32_t>>>>>> &WL,
                                                          const shared_ptr<unordered_set<shared_ptr<abstract_bipartite_edge>>> &candidate_edge_set,
                                                          uint32_t k,
                                                          const shared_ptr<thread_pool> &pool)
    {
        auto rank_id = make_shared<uint32_t>(0);
        set_edge_rank(candidate_edge_set, edge_rank_map, rank_id);
        auto thread_number = pool->get_thread_number();
        auto global_mutex = make_shared<mutex>();
        auto location_vector = pool->split_task(candidate_edge_set);
        for(uint32_t i = 0; i < thread_number; ++i){
            pool->submit_task([=]{
                for(auto iter = *location_vector->at(i); iter!=*location_vector->at(i + 1); ++iter){
                    auto &e1 = *iter;
                    auto l1 = e1->get_left_vertex_id();
                    auto r1 = e1->get_right_vertex_id();


                    for (const auto &r2: *WL->at(l1)->at(k)) {
                        auto e2 = G->get_edge(l1, r2);
                        if(r2 == r1 || edge_rank_map->at(e2) < edge_rank_map->at(e1)){
                            continue;
                        }
                        for (const auto &l2: *WL->at(r1)->at(k)) {
                            auto e3 = G->get_edge(l2, r1);
                            if(l2 == l1 || edge_rank_map->at(e3) < edge_rank_map->at(e1) )
                            {
                                continue;
                            }

                            auto e4 = G->get_edge(l2, r2);
                            if(!e4 || edge_wing_map->at(e4) < k || edge_rank_map->at(e4) < edge_rank_map->at(e1)){
                                continue;
                            }

                            if(!candidate_edge_set->count(e2))
                            {
                                edge_mutex_map->at(e2)->lock();
                                ++WS->at(e2);
                                edge_mutex_map->at(e2)->unlock();
                            }

                            if(!candidate_edge_set->count(e3))
                            {
                                edge_mutex_map->at(e3)->lock();
                                ++WS->at(e3);
                                edge_mutex_map->at(e3)->unlock();
                            }

                            if(edge_wing_map->at(e4) == k && !candidate_edge_set->count(e4))
                            {
                                edge_mutex_map->at(e4)->lock();
                                ++WS->at(e4);
                                edge_mutex_map->at(e4)->unlock();
                            }
                        }
                        for (auto iter3 = WL->at(r1)->lower_bound(k + 1); iter3 != WL->at(r1)->end(); ++iter3) {
                            for (const auto &l2: *iter3->second) {
                                auto e3 = G->get_edge(l2, r1);
//                                    if(l2 == l1 || edge_rank_map->at(e3) < edge_rank_map->at(e1) )
//                                    {
//                                        continue;
//                                    }
                                auto e4 = G->get_edge(l2, r2);
                                if(!e4 || edge_wing_map->at(e4) < k || edge_rank_map->at(e4) < edge_rank_map->at(e1)){
                                    continue;
                                }
                                if(edge_wing_map->at(e2) == k && !candidate_edge_set->count(e2))
                                {
                                    edge_mutex_map->at(e2)->lock();
                                    ++WS->at(e2);
                                    edge_mutex_map->at(e2)->unlock();
                                }

//                                    if(edge_wing_map->at(e3) == k && !candidate_edge_set->count(e3))
//                                    {
//                                        edge_mutex_map->at(e3)->lock();
//                                        ++WS->at(e3);
//                                        edge_mutex_map->at(e3)->unlock();
//                                    }

                                if(edge_wing_map->at(e4) == k && !candidate_edge_set->count(e4))
                                {
                                    edge_mutex_map->at(e4)->lock();
                                    ++WS->at(e4);
                                    edge_mutex_map->at(e4)->unlock();
                                }
                            }
                        }
                    }
                    for(auto iter2 = WL->at(l1)->lower_bound(k + 1);iter2!=WL->at(l1)->end();++iter2) {
                        for (const auto &r2: *iter2->second) {
                            auto e2 = G->get_edge(l1, r2);
//                            if(r2 == r1 || edge_rank_map->at(e2) < edge_rank_map->at(e1)){
//                                continue;
//                            }
                            for (const auto &l2: *WL->at(r1)->at(k)) {
                                auto e3 = G->get_edge(l2, r1);
                                if(l2 == l1 || edge_rank_map->at(e3) < edge_rank_map->at(e1))
                                {
                                    continue;
                                }
                                auto e4 = G->get_edge(l2, r2);
                                if(!e4 || edge_wing_map->at(e4) < k || edge_rank_map->at(e4) < edge_rank_map->at(e1)){
                                    continue;
                                }
//                                    if(edge_wing_map->at(e2) == k && !candidate_edge_set->count(e2))
//                                    {
//                                        edge_mutex_map->at(e2)->lock();
//                                        ++WS->at(e2);
//                                        edge_mutex_map->at(e2)->unlock();
//                                    }

                                if(edge_wing_map->at(e3) == k && !candidate_edge_set->count(e3))
                                {
                                    edge_mutex_map->at(e3)->lock();
                                    ++WS->at(e3);
                                    edge_mutex_map->at(e3)->unlock();
                                }

                                if(edge_wing_map->at(e4) == k && !candidate_edge_set->count(e4))
                                {
                                    edge_mutex_map->at(e4)->lock();
                                    ++WS->at(e4);
                                    edge_mutex_map->at(e4)->unlock();
                                }
                            }
                            for (auto iter3 = WL->at(r1)->lower_bound(k + 1); iter3 != WL->at(r1)->end(); ++iter3) {
                                for (const auto &l2: *iter3->second) {
                                    auto e3 = G->get_edge(l2, r1);
//                                    if(l2 == l1 || edge_rank_map->at(e3) < edge_rank_map->at(e1))
//                                    {
//                                        continue;
//                                    }
                                    auto e4 = G->get_edge(l2, r2);
                                    if(!e4 || edge_wing_map->at(e4) < k || edge_rank_map->at(e4) < edge_rank_map->at(e1)){
                                        continue;
                                    }
//                                    if(edge_wing_map->at(e2) == k && !candidate_edge_set->count(e2))
//                                    {
//                                        edge_mutex_map->at(e2)->lock();
//                                        ++WS->at(e2);
//                                        edge_mutex_map->at(e2)->unlock();
//                                    }

//                                    if(edge_wing_map->at(e3) == k && !candidate_edge_set->count(e3))
//                                    {
//                                        edge_mutex_map->at(e3)->lock();
//                                        ++WS->at(e3);
//                                        edge_mutex_map->at(e3)->unlock();
//                                    }

                                    if(edge_wing_map->at(e4) == k && !candidate_edge_set->count(e4))
                                    {
                                        edge_mutex_map->at(e4)->lock();
                                        ++WS->at(e4);
                                        edge_mutex_map->at(e4)->unlock();
                                    }
                                }
                            }
                        }
                    }
                }
            });
        }
        pool->barrier();
        for(const auto &e:*candidate_edge_set){
            edge_rank_map->at(e) = UINT32_MAX;
        }
    }


    void quasi_wing_maintenance::update_single_wing(const shared_ptr<abstract_bipartite_graph> &G,
                                                    const shared_ptr<unordered_map<shared_ptr<abstract_bipartite_edge>, shared_ptr<mutex>>> &edge_mutex_map,
                                                    const shared_ptr<unordered_map<shared_ptr<abstract_bipartite_edge>, uint32_t>> &edge_rank_map,
                                                    const shared_ptr<unordered_map<shared_ptr<abstract_bipartite_edge>, uint32_t>> &edge_wing_map,
                                                    const shared_ptr<unordered_set<shared_ptr<abstract_bipartite_edge>>> &edge_set,
                                                    const shared_ptr<unordered_map<shared_ptr<abstract_bipartite_edge>, uint32_t>> &WS,
                                                    const shared_ptr<unordered_set<shared_ptr<abstract_bipartite_edge>>>& evicted_edge_set,
                                                    uint32_t k,
                                                    const shared_ptr<thread_pool> &pool) {
        auto thread_number = pool->get_thread_number();
        auto global_mutex = make_shared<mutex>();

        auto current_edge_set = make_shared<unordered_set<shared_ptr<abstract_bipartite_edge>>>();
        {
            auto removed_edge_set = make_shared<unordered_set<shared_ptr<abstract_bipartite_edge>>>();
            auto location_vector = pool->split_task(edge_set);
            for(uint32_t i = 0; i < thread_number; ++i){
                pool->submit_task([=]{
                    auto sub_current_edge_set = make_shared<unordered_set<shared_ptr<abstract_bipartite_edge>>>();
                    auto sub_removed_edge_set = make_shared<unordered_set<shared_ptr<abstract_bipartite_edge>>>();

                    auto &sub_begin = *location_vector->at(i);
                    auto &sub_end = *location_vector->at(i + 1);

                    for(auto iter = sub_begin; iter != sub_end; ++iter){
                        auto &e = *iter;
                        if(edge_wing_map->at(e) > k){
                            sub_removed_edge_set->insert(e);
                        }else{
                            if(edge_wing_map->at(e) == k){
                                WS->at(e) = edge_support_computation(G, e, edge_wing_map, k);

                                if (WS->at(e) < k) {
                                    sub_current_edge_set->insert(e);
                                }
                            }
                        }
                    }
                    global_mutex->lock();
                    current_edge_set->merge(*sub_current_edge_set);
                    removed_edge_set->merge(*sub_removed_edge_set);
                    global_mutex->unlock();
                });
            }
            pool->barrier();
            for(const auto &e:*removed_edge_set){
                edge_set->erase(e);
            }
        }

        auto clear_set = make_shared<unordered_set<shared_ptr<abstract_bipartite_edge>>>();
        /**
         * @brief a set of edges revising their rank
         */
        auto rank_id = make_shared<uint32_t>(0);
        while(!current_edge_set->empty()){
            set_edge_rank(current_edge_set, edge_rank_map, rank_id);
            auto next_edge_set = make_shared<unordered_set<shared_ptr<abstract_bipartite_edge>>>();
            auto location_vector = pool->split_task(current_edge_set);
            for(uint32_t i = 0; i < thread_number; ++i){
                pool->submit_task([=]{
                    auto sub_next_edge_set = make_shared<unordered_set<shared_ptr<abstract_bipartite_edge>>>();
                    auto sub_clear_set = make_shared<unordered_set<shared_ptr<abstract_bipartite_edge>>>();

                    auto &sub_begin = *location_vector->at(i);
                    auto &sub_end = *location_vector->at(i + 1);

                    for(auto iter = sub_begin; iter != sub_end; ++iter){
                        auto &e1 = *iter;

                        auto l1 = e1->get_left_vertex_id();
                        auto r1 = e1->get_right_vertex_id();
                        for (const auto &[r2, e2]:*G->get_left_vertex(l1)->get_edge_map()) {
                            if (r2 == r1 || edge_wing_map->at(e2) < k ||
                                edge_rank_map->at(e2) < edge_rank_map->at(e1)) {
                                continue;
                            }
                            for (const auto &[l2, e3]: *G->get_right_vertex(r1)->get_edge_map()) {
                                if (l2 == l1 || edge_wing_map->at(e3) < k ||
                                    edge_rank_map->at(e3) < edge_rank_map->at(e1)) {
                                    continue;
                                }

                                auto e4 = G->get_edge(l2, r2);
                                if (!e4 || edge_wing_map->at(e4) < k || edge_rank_map->at(e4) < edge_rank_map->at(e1)) {
                                    continue;
                                }


                                if (!current_edge_set->count(e2) && edge_wing_map->at(e2) == k) {
                                    sub_clear_set->insert(e2);
                                    edge_mutex_map->at(e2)->lock();
                                    if (WS->at(e2) == 0) {
                                        WS->at(e2) = edge_support_computation(G, e2, edge_wing_map, k);
                                    }
                                    --WS->at(e2);
                                    edge_mutex_map->at(e2)->unlock();

                                    if (WS->at(e2) < k) {
                                        sub_next_edge_set->insert(e2);
                                    }
                                }


                                if (!current_edge_set->count(e3) && edge_wing_map->at(e3) == k) {
                                    sub_clear_set->insert(e3);
                                    edge_mutex_map->at(e3)->lock();
                                    if (WS->at(e3) == 0) {
                                        WS->at(e3) = edge_support_computation(G, e3, edge_wing_map, k);
                                    }
                                    --WS->at(e3);
                                    edge_mutex_map->at(e3)->unlock();

                                    if (WS->at(e3) < k) {
                                        sub_next_edge_set->insert(e3);
                                    }
                                }

                                if (!current_edge_set->count(e4) && edge_wing_map->at(e4) == k) {
                                    sub_clear_set->insert(e4);
                                    edge_mutex_map->at(e4)->lock();
                                    if (WS->at(e4) == 0) {
                                        WS->at(e4) = edge_support_computation(G, e4, edge_wing_map, k);
                                    }
                                    --WS->at(e4);
                                    edge_mutex_map->at(e4)->unlock();

                                    if (WS->at(e4) < k) {
                                        sub_next_edge_set->insert(e4);
                                    }
                                }
                            }
                        }
                    }

                    global_mutex->lock();
                    next_edge_set->merge(*sub_next_edge_set);
                    clear_set->merge(*sub_clear_set);
                    global_mutex->unlock();
                });
            }
            pool->barrier();
            evicted_edge_set->merge(*current_edge_set);
            swap(*current_edge_set, *next_edge_set);
        }

        for(const auto &e:*clear_set){
            WS->at(e) = 0;
        }
        for(const auto &e:*evicted_edge_set){
            edge_rank_map->at(e) = UINT32_MAX;
        }
    }

    void quasi_wing_maintenance::update_single_wing2(const shared_ptr<abstract_bipartite_graph> &G,
                                                     const shared_ptr<unordered_map<shared_ptr<abstract_bipartite_edge>, shared_ptr<mutex>>> &edge_mutex_map,
                                                     const shared_ptr<unordered_map<shared_ptr<abstract_bipartite_edge>, uint32_t>> &edge_rank_map,
                                                     const shared_ptr<unordered_map<shared_ptr<abstract_bipartite_edge>, uint32_t>> &edge_wing_map,
                                                     const shared_ptr<unordered_set<shared_ptr<abstract_bipartite_edge>>> &edge_set,
                                                     const shared_ptr<unordered_map<shared_ptr<abstract_bipartite_edge>, uint32_t>> &WS,
                                                     const shared_ptr<unordered_set<shared_ptr<abstract_bipartite_edge>>>& evicted_edge_set,
                                                     uint32_t k,
                                                     const shared_ptr<thread_pool> &pool) {
        auto thread_number = pool->get_thread_number();
        auto global_mutex = make_shared<mutex>();

        auto current_edge_set = make_shared<unordered_set<shared_ptr<abstract_bipartite_edge>>>();
        {
            auto removed_edge_set = make_shared<unordered_set<shared_ptr<abstract_bipartite_edge>>>();
            auto location_vector = pool->split_task(edge_set);
            for(uint32_t i = 0; i < thread_number; ++i){
                pool->submit_task([=]{
                    auto sub_current_edge_set = make_shared<unordered_set<shared_ptr<abstract_bipartite_edge>>>();
                    auto sub_removed_edge_set = make_shared<unordered_set<shared_ptr<abstract_bipartite_edge>>>();

                    auto &sub_begin = *location_vector->at(i);
                    auto &sub_end = *location_vector->at(i + 1);

                    for(auto iter = sub_begin; iter != sub_end; ++iter){
                        auto &e = *iter;
                        if(edge_wing_map->at(e) > k){
                            sub_removed_edge_set->insert(e);
                        }else{
                            if(edge_wing_map->at(e) == k){
                                //WS->at(e) = edge_support_computation(G, e, edge_wing_map, k);

                                if (WS->at(e) < k) {
                                    sub_current_edge_set->insert(e);
                                }
                            }
                        }
                    }
                    global_mutex->lock();
                    current_edge_set->merge(*sub_current_edge_set);
                    removed_edge_set->merge(*sub_removed_edge_set);
                    global_mutex->unlock();
                });
            }
            pool->barrier();
            for(const auto &e:*removed_edge_set){
                edge_set->erase(e);
            }
        }
        /**
         * @brief a set of edges revising their rank
         */
        auto rank_id = make_shared<uint32_t>(0);
        auto current_task_vector = make_shared<vector<shared_ptr<unordered_set<shared_ptr<abstract_bipartite_edge>>>>>();
        for (uint32_t i = 0; i < thread_number; i++) {
            current_task_vector->emplace_back(make_shared<unordered_set<shared_ptr<abstract_bipartite_edge>>>());
        }
        while(!current_edge_set->empty()){
            set_edge_rank(current_edge_set, edge_rank_map, rank_id);
            auto next_edge_set = make_shared<unordered_set<shared_ptr<abstract_bipartite_edge>>>();
            auto location_vector = pool->split_task(current_edge_set);
            for(uint32_t i = 0; i < thread_number; ++i){
                pool->submit_task([=]{
                    auto sub_next_edge_set = make_shared<unordered_set<shared_ptr<abstract_bipartite_edge>>>();


                    for(auto iter = *location_vector->at(i); iter!=*location_vector->at(i + 1); ++iter){
                        auto& e1 = *iter;
                        auto l1 = e1->get_left_vertex_id();
                        auto r1 = e1->get_right_vertex_id();
                        for (const auto &[r2, e2]:*G->get_left_vertex(l1)->get_edge_map()) {
                            if (r2 == r1 || edge_wing_map->at(e2) < k || edge_rank_map->at(e2) < edge_rank_map->at(e1)) {
                                continue;
                            }
                            for (const auto &[l2, e3]:*G->get_right_vertex(r1)->get_edge_map()) {
                                if (l2 == l1 || edge_wing_map->at(e3) < k || edge_rank_map->at(e3) < edge_rank_map->at(e1)) {
                                    continue;
                                }

                                auto e4 = G->get_edge(l2, r2);
                                if (!e4 || edge_wing_map->at(e4) < k || edge_rank_map->at(e4) < edge_rank_map->at(e1)) {
                                    continue;
                                }


                                if(!current_edge_set->count(e2)){
                                    if(edge_wing_map->at(e2) == k){
                                        edge_mutex_map->at(e2)->lock();
                                        --WS->at(e2);
                                        edge_mutex_map->at(e2)->unlock();

                                        if (WS->at(e2) < k) {
                                            sub_next_edge_set->insert(e2);
                                        }
                                    }
                                }


                                if(!current_edge_set->count(e3)){
                                    if(edge_wing_map->at(e3) == k){
                                        edge_mutex_map->at(e3)->lock();
                                        --WS->at(e3);
                                        edge_mutex_map->at(e3)->unlock();

                                        if (WS->at(e3) < k) {
                                            sub_next_edge_set->insert(e3);
                                        }
                                    }
                                }

                                if(!current_edge_set->count(e4)){
                                    if(edge_wing_map->at(e4) == k){
                                        edge_mutex_map->at(e4)->lock();
                                        --WS->at(e4);
                                        edge_mutex_map->at(e4)->unlock();

                                        if (WS->at(e4) < k) {
                                            sub_next_edge_set->insert(e4);
                                        }
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
            evicted_edge_set->merge(*current_edge_set);
            swap(*current_edge_set, *next_edge_set);
        }

        for(const auto &e:*evicted_edge_set){
            edge_rank_map->at(e) = UINT32_MAX;
        }
    }

    void quasi_wing_maintenance::update_single_wing3(const shared_ptr<abstract_bipartite_graph> &G,
                                                     const shared_ptr<BE_index>& bloom_index,
                                                     const shared_ptr<unordered_map<shared_ptr<abstract_bipartite_edge>, shared_ptr<mutex>>> &edge_mutex_map,
                                                     const shared_ptr<unordered_map<shared_ptr<abstract_bipartite_edge>, uint32_t>> &edge_rank_map,
                                                     const shared_ptr<unordered_map<shared_ptr<abstract_bipartite_edge>, uint32_t>> &edge_wing_map,
                                                     const shared_ptr<unordered_set<shared_ptr<abstract_bipartite_edge>>> &edge_set,
                                                     const shared_ptr<unordered_map<shared_ptr<abstract_bipartite_edge>, uint32_t>> &WS,
                                                     const shared_ptr<unordered_set<shared_ptr<abstract_bipartite_edge>>>& evicted_edge_set,
                                                     uint32_t k,
                                                     const shared_ptr<thread_pool> &pool) {
        auto thread_number = pool->get_thread_number();
        auto global_mutex = make_shared<mutex>();

        auto current_edge_set = make_shared<unordered_set<shared_ptr<abstract_bipartite_edge>>>();
        {
            auto removed_edge_set = make_shared<unordered_set<shared_ptr<abstract_bipartite_edge>>>();
            auto location_vector = pool->split_task(edge_set);
            for(uint32_t i = 0; i < thread_number; ++i){
                pool->submit_task([=]{
                    auto sub_current_edge_set = make_shared<unordered_set<shared_ptr<abstract_bipartite_edge>>>();
                    auto sub_removed_edge_set = make_shared<unordered_set<shared_ptr<abstract_bipartite_edge>>>();

                    for(auto iter = *location_vector->at(i); iter != *location_vector->at(i + 1); ++iter){
                        auto &e = *iter;
                        if(edge_wing_map->at(e) > k){
                            sub_removed_edge_set->insert(e);
                        }else{
                            if(edge_wing_map->at(e) == k){
                                //WS->at(e) = edge_support_computation(G, e, edge_wing_map, k);

                                if (WS->at(e) < k) {
                                    sub_current_edge_set->insert(e);
                                }
                            }
                        }
                    }
                    global_mutex->lock();
                    current_edge_set->merge(*sub_current_edge_set);
                    removed_edge_set->merge(*sub_removed_edge_set);
                    global_mutex->unlock();
                });
            }
            pool->barrier();
            for(const auto &e:*removed_edge_set){
                edge_set->erase(e);
            }
        }
        /**
         * @brief a set of edges revising their rank
         */
        auto rank_id = make_shared<uint32_t>(0);
        while(!current_edge_set->empty()){
            set_edge_rank(current_edge_set, edge_rank_map, rank_id);
            auto next_edge_set = make_shared<unordered_set<shared_ptr<abstract_bipartite_edge>>>();
            auto location_vector = pool->split_task(current_edge_set);
            for(uint32_t i = 0; i < thread_number; ++i){
                pool->submit_task([=] {
                    auto sub_next_edge_set = make_shared<unordered_set<shared_ptr<abstract_bipartite_edge>>>();

                    for(auto iter = *location_vector->at(i); iter != *location_vector->at(i + 1); ++iter){
                        auto &e1 = *iter;
                        auto l1 = e1->get_left_vertex_id();
                        auto r1 = e1->get_right_vertex_id();

                        auto bloom_set = bloom_index->get_bloom_set(e1);
                        if(bloom_set){
                            for (const auto &B: *bloom_set) {
                                if (B->get_butterfly_count() > 0) {
                                    auto e2 = B->get_twin(e1);

                                    if (edge_wing_map->at(e2) < k || edge_rank_map->at(e2) < edge_rank_map->at(e1)) {
                                        continue;
                                    }

                                    auto visited_set = make_shared<unordered_set<shared_ptr<abstract_bipartite_edge>>>();

                                    for (const auto &[e3, e4]: *B->get_edge_map()) {
                                        if (e3 == e1 || e3 == e2 || visited_set->count(e3)) {
                                            continue;
                                        }
                                        visited_set->insert(e4);

//                                        if ((e3->get_left_vertex_id() == e4->get_left_vertex_id() &&
//                                             e3->get_right_vertex_id() < e4->get_right_vertex_id())
//                                            || (e3->get_right_vertex_id() == e4->get_right_vertex_id() &&
//                                                e3->get_left_vertex_id() < e4->get_left_vertex_id())) {
//                                            continue;
//                                        }


                                        if (edge_wing_map->at(e3) < k ||
                                            edge_rank_map->at(e3) < edge_rank_map->at(e1)) {
                                            continue;
                                        }

                                        if (edge_wing_map->at(e4) < k ||
                                            edge_rank_map->at(e4) < edge_rank_map->at(e1)) {
                                            continue;
                                        }



//                                        if(e3->get_left_vertex_id()!= l1 && e3->get_right_vertex_id() != r1){
//                                            if(visited_set->count(e3)){
//                                                continue;
//                                            }
//                                            visited_set->insert(e3);
//                                        }
//                                        else{
//                                            if(visited_set->count(e4)){
//                                                continue;
//                                            }
//                                            visited_set->insert(e4);
//                                        }

                                        if (!current_edge_set->count(e2) && edge_wing_map->at(e2) == k) {
                                            edge_mutex_map->at(e2)->lock();
                                            --WS->at(e2);
                                            edge_mutex_map->at(e2)->unlock();

                                            if (WS->at(e2) < k) {
                                                sub_next_edge_set->insert(e2);
                                            }
                                        }


                                        if (!current_edge_set->count(e3) && edge_wing_map->at(e3) == k) {
                                            edge_mutex_map->at(e3)->lock();
                                            --WS->at(e3);
                                            edge_mutex_map->at(e3)->unlock();

                                            if (WS->at(e3) < k) {
                                                sub_next_edge_set->insert(e3);
                                            }
                                        }

                                        if (!current_edge_set->count(e4) && edge_wing_map->at(e4) == k) {

                                            edge_mutex_map->at(e4)->lock();
                                            --WS->at(e4);
                                            edge_mutex_map->at(e4)->unlock();

                                            if (WS->at(e4) < k) {
                                                sub_next_edge_set->insert(e4);
                                            }
                                        }
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
            evicted_edge_set->merge(*current_edge_set);
            swap(*current_edge_set, *next_edge_set);
        }
    }

    void quasi_wing_maintenance::update_single_wing4(const shared_ptr<abstract_bipartite_graph> &G,
                                                     const shared_ptr<unordered_map<shared_ptr<abstract_bipartite_edge>, shared_ptr<mutex>>> &edge_mutex_map,
                                                     const shared_ptr<unordered_map<shared_ptr<abstract_bipartite_edge>, uint32_t>> &edge_rank_map,
                                                     const shared_ptr<unordered_map<shared_ptr<abstract_bipartite_edge>, uint32_t>> &edge_wing_map,
                                                     const shared_ptr<unordered_set<shared_ptr<abstract_bipartite_edge>>> &edge_set,
                                                     const shared_ptr<unordered_map<shared_ptr<abstract_bipartite_edge>, uint32_t>> &WS,
                                                     const shared_ptr<unordered_map<uint32_t,shared_ptr<map<uint32_t,shared_ptr<unordered_set<uint32_t>>>>>> &WL,
                                                     const shared_ptr<unordered_set<shared_ptr<abstract_bipartite_edge>>>& evicted_edge_set,
                                                     uint32_t k,
                                                     const shared_ptr<thread_pool> &pool) {
        auto thread_number = pool->get_thread_number();
        auto global_mutex = make_shared<mutex>();

        auto current_edge_set = make_shared<unordered_set<shared_ptr<abstract_bipartite_edge>>>();
        {
            auto removed_edge_set = make_shared<unordered_set<shared_ptr<abstract_bipartite_edge>>>();
            auto location_vector = pool->split_task(edge_set);
            for(uint32_t i = 0; i < thread_number; ++i){
                pool->submit_task([=]{
                    auto sub_current_edge_set = make_shared<unordered_set<shared_ptr<abstract_bipartite_edge>>>();
                    auto sub_removed_edge_set = make_shared<unordered_set<shared_ptr<abstract_bipartite_edge>>>();

                    auto &sub_begin = *location_vector->at(i);
                    auto &sub_end = *location_vector->at(i + 1);

                    for(auto iter = sub_begin; iter != sub_end; ++iter){
                        auto &e = *iter;
                        if(edge_wing_map->at(e) > k){
                            sub_removed_edge_set->insert(e);
                        }else{
                            if(edge_wing_map->at(e) == k){
                                //WS->at(e) = edge_support_computation(G, e, edge_wing_map, k);

                                if (WS->at(e) < k) {
                                    sub_current_edge_set->insert(e);
                                }
                            }
                        }
                    }
                    global_mutex->lock();
                    current_edge_set->merge(*sub_current_edge_set);
                    removed_edge_set->merge(*sub_removed_edge_set);
                    global_mutex->unlock();
                });
            }
            pool->barrier();
            for(const auto &e:*removed_edge_set){
                edge_set->erase(e);
            }
        }
        /**
         * @brief a set of edges revising their rank
         */
        auto rank_id = make_shared<uint32_t>(0);
        auto current_task_vector = make_shared<vector<shared_ptr<unordered_set<shared_ptr<abstract_bipartite_edge>>>>>();
        for (uint32_t i = 0; i < thread_number; i++) {
            current_task_vector->emplace_back(make_shared<unordered_set<shared_ptr<abstract_bipartite_edge>>>());
        }
        while(!current_edge_set->empty()){
            set_edge_rank(current_edge_set, edge_rank_map, rank_id);
            auto next_edge_set = make_shared<unordered_set<shared_ptr<abstract_bipartite_edge>>>();
            auto location_vector = pool->split_task(current_edge_set);
            for(uint32_t i = 0; i < thread_number; ++i){
                pool->submit_task([=] {
                    auto sub_next_edge_set = make_shared<unordered_set<shared_ptr<abstract_bipartite_edge>>>();

                    for (auto iter = *location_vector->at(i); iter != *location_vector->at(i + 1); ++iter) {
                        auto &e1 = *iter;
                        auto l1 = e1->get_left_vertex_id();
                        auto r1 = e1->get_right_vertex_id();

                        for (const auto &r2: *WL->at(l1)->at(k)) {
                            auto e2 = G->get_edge(l1, r2);
                            if (r2 == r1 || edge_rank_map->at(e2) < edge_rank_map->at(e1)) {
                                continue;
                            }
                            for (const auto &l2: *WL->at(r1)->at(k)) {
                                auto e3 = G->get_edge(l2, r1);
                                if (l2 == l1 || edge_rank_map->at(e3) < edge_rank_map->at(e1)) {
                                    continue;
                                }

                                auto e4 = G->get_edge(l2, r2);
                                if (!e4 || edge_wing_map->at(e4) < k || edge_rank_map->at(e4) < edge_rank_map->at(e1)) {
                                    continue;
                                }

                                if (!current_edge_set->count(e2)) {
                                    edge_mutex_map->at(e2)->lock();
                                    --WS->at(e2);
                                    edge_mutex_map->at(e2)->unlock();

                                    if (WS->at(e2) < k) {
                                        sub_next_edge_set->insert(e2);
                                    }
                                }

                                if (!current_edge_set->count(e3)) {
                                    edge_mutex_map->at(e3)->lock();
                                    --WS->at(e3);
                                    edge_mutex_map->at(e3)->unlock();

                                    if (WS->at(e3) < k) {
                                        sub_next_edge_set->insert(e3);
                                    }
                                }

                                if (edge_wing_map->at(e4) == k && !current_edge_set->count(e4)) {
                                    edge_mutex_map->at(e4)->lock();
                                    --WS->at(e4);
                                    edge_mutex_map->at(e4)->unlock();

                                    if (WS->at(e4) < k) {
                                        sub_next_edge_set->insert(e4);
                                    }
                                }
                            }
                            for (auto iter3 = WL->at(r1)->lower_bound(k + 1); iter3 != WL->at(r1)->end(); ++iter3) {
                                for (const auto &l2: *iter3->second) {
                                    auto e3 = G->get_edge(l2, r1);
//                                    if(l2 == l1 || edge_rank_map->at(e3) < edge_rank_map->at(e1) )
//                                    {
//                                        continue;
//                                    }
                                    auto e4 = G->get_edge(l2, r2);
                                    if (!e4 || edge_wing_map->at(e4) < k ||
                                        edge_rank_map->at(e4) < edge_rank_map->at(e1)) {
                                        continue;
                                    }
                                    if (!current_edge_set->count(e2)) {
                                        edge_mutex_map->at(e2)->lock();
                                        --WS->at(e2);
                                        edge_mutex_map->at(e2)->unlock();

                                        if (WS->at(e2) < k) {
                                            sub_next_edge_set->insert(e2);
                                        }
                                    }

//                                    if (!current_edge_set->count(e3)) {
//                                        edge_mutex_map->at(e3)->lock();
//                                        --WS->at(e3);
//                                        edge_mutex_map->at(e3)->unlock();
//
//                                        if (WS->at(e3) < k) {
//                                            sub_next_edge_set->insert(e3);
//                                        }
//                                    }

                                    if (edge_wing_map->at(e4) == k && !current_edge_set->count(e4)) {
                                        edge_mutex_map->at(e4)->lock();
                                        --WS->at(e4);
                                        edge_mutex_map->at(e4)->unlock();

                                        if (WS->at(e4) < k) {
                                            sub_next_edge_set->insert(e4);
                                        }
                                    }
                                }
                            }
                        }
                        for (auto iter2 = WL->at(l1)->lower_bound(k + 1); iter2 != WL->at(l1)->end(); ++iter2) {
                            for (const auto &r2: *iter2->second) {
                                auto e2 = G->get_edge(l1, r2);
//                            if(r2 == r1 || edge_rank_map->at(e2) < edge_rank_map->at(e1)){
//                                continue;
//                            }
                                for (const auto &l2: *WL->at(r1)->at(k)) {
                                    auto e3 = G->get_edge(l2, r1);
                                    if (l2 == l1 || edge_rank_map->at(e3) < edge_rank_map->at(e1)) {
                                        continue;
                                    }
                                    auto e4 = G->get_edge(l2, r2);
                                    if (!e4 || edge_wing_map->at(e4) < k ||
                                        edge_rank_map->at(e4) < edge_rank_map->at(e1)) {
                                        continue;
                                    }
//                                    if (!current_edge_set->count(e2)) {
//                                        edge_mutex_map->at(e2)->lock();
//                                        --WS->at(e2);
//                                        edge_mutex_map->at(e2)->unlock();
//
//                                        if (WS->at(e2) < k) {
//                                            sub_next_edge_set->insert(e2);
//                                        }
//                                    }


                                    if (!current_edge_set->count(e3)) {
                                        edge_mutex_map->at(e3)->lock();
                                        --WS->at(e3);
                                        edge_mutex_map->at(e3)->unlock();

                                        if (WS->at(e3) < k) {
                                            sub_next_edge_set->insert(e3);
                                        }
                                    }

                                    if (edge_wing_map->at(e4) == k && !current_edge_set->count(e4)) {
                                        edge_mutex_map->at(e4)->lock();
                                        --WS->at(e4);
                                        edge_mutex_map->at(e4)->unlock();

                                        if (WS->at(e4) < k) {
                                            sub_next_edge_set->insert(e4);
                                        }
                                    }
                                }
                                for (auto iter3 = WL->at(r1)->lower_bound(k + 1); iter3 != WL->at(r1)->end(); ++iter3) {
                                    for (const auto &l2: *iter3->second) {
                                        auto e3 = G->get_edge(l2, r1);
//                                    if(l2 == l1 || edge_rank_map->at(e3) < edge_rank_map->at(e1))
//                                    {
//                                        continue;
//                                    }
                                        auto e4 = G->get_edge(l2, r2);
                                        if (!e4 || edge_wing_map->at(e4) < k ||
                                            edge_rank_map->at(e4) < edge_rank_map->at(e1)) {
                                            continue;
                                        }
//                                    if (!current_edge_set->count(e2)) {
//                                        edge_mutex_map->at(e2)->lock();
//                                        --WS->at(e2);
//                                        edge_mutex_map->at(e2)->unlock();
//
//                                        if (WS->at(e2) < k) {
//                                            sub_next_edge_set->insert(e2);
//                                        }
//                                    }

//                                    if (!current_edge_set->count(e3)) {
//                                        edge_mutex_map->at(e3)->lock();
//                                        --WS->at(e3);
//                                        edge_mutex_map->at(e3)->unlock();
//
//                                        if (WS->at(e3) < k) {
//                                            sub_next_edge_set->insert(e3);
//                                        }
//                                    }

                                        if (edge_wing_map->at(e4) == k && !current_edge_set->count(e4)) {
                                            edge_mutex_map->at(e4)->lock();
                                            --WS->at(e4);
                                            edge_mutex_map->at(e4)->unlock();

                                            if (WS->at(e4) < k) {
                                                sub_next_edge_set->insert(e4);
                                            }
                                        }
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
            evicted_edge_set->merge(*current_edge_set);
            swap(*current_edge_set, *next_edge_set);
        }

        for(const auto &e:*evicted_edge_set){
            edge_rank_map->at(e) = UINT32_MAX;
        }
    }

    void quasi_wing_maintenance::vertex_priority_computation(const shared_ptr<abstract_bipartite_graph> &G,
                                                             const shared_ptr<unordered_map<uint32_t,uint32_t>>& vertex_priority_map,
                                                             const shared_ptr<thread_pool> &pool) {
        auto global_mutex = make_shared<mutex>();
        auto max_degree = make_shared<uint32_t>(0);
        {
            {
                auto left_vertex_map = G->get_left_vertex_map();
                auto location_vector = pool->split_task(left_vertex_map);
                for (uint32_t i = 0; i < pool->get_thread_number(); ++i) {
                    pool->submit_task([=] {
                        uint32_t sub_max = 0;
                        for (auto iter = *location_vector->at(i); iter != *location_vector->at(i + 1); ++iter) {
                            auto [l, l_vertex] = *iter;
                            if (l_vertex->get_degree() > sub_max) {
                                sub_max = l_vertex->get_degree();
                            }
                        }

                        global_mutex->lock();
                        if (sub_max > *max_degree) {
                            *max_degree = sub_max;
                        }
                        global_mutex->unlock();
                    });
                }
            }

            {
                auto right_vertex_map = G->get_right_vertex_map();
                auto location_vector = pool->split_task(right_vertex_map);
                for (uint32_t i = 0; i < pool->get_thread_number(); ++i) {
                    pool->submit_task([=] {
                        uint32_t sub_max = 0;
                        for (auto iter = *location_vector->at(i); iter != *location_vector->at(i + 1); ++iter) {
                            auto [r, r_vertex] = *iter;
                            if (r_vertex->get_degree() > sub_max) {
                                sub_max = r_vertex->get_degree();
                            }
                        }

                        global_mutex->lock();
                        if (sub_max > *max_degree) {
                            *max_degree = sub_max;
                        }
                        global_mutex->unlock();
                    });
                }
            }
            pool->barrier();
        }

        auto degree_vector = make_shared<vector<shared_ptr<set<uint32_t>>>>(*max_degree + 1);
        auto degree_mutex_vector = make_shared<vector<shared_ptr<mutex>>>(*max_degree + 1);
        {
            auto location_vector = pool->split_task(degree_mutex_vector);
            for (uint32_t i = 0; i < pool->get_thread_number(); ++i) {
                pool->submit_task([=] {
                    for (auto iter = *location_vector->at(i); iter != *location_vector->at(i + 1); ++iter) {
                        *iter = make_shared<mutex>();
                    }
                });
            }
            pool->barrier();
        }
        {
            {
                auto left_vertex_map = G->get_left_vertex_map();
                auto location_vector = pool->split_task(left_vertex_map);
                for (uint32_t i = 0; i < pool->get_thread_number(); ++i) {
                    pool->submit_task([=] {
                        for (auto iter = *location_vector->at(i); iter != *location_vector->at(i + 1); ++iter) {
                            auto [l, l_vertex] = *iter;
                            auto degree = l_vertex->get_degree();
                            degree_mutex_vector->at(degree)->lock();
                            if (!degree_vector->at(degree)) {
                                degree_vector->at(degree) = make_shared<set<uint32_t>>();
                            }
                            degree_vector->at(degree)->insert(l);
                            degree_mutex_vector->at(degree)->unlock();
                        }
                    });
                }
            }

            {
                auto right_vertex_map = G->get_right_vertex_map();
                auto location_vector = pool->split_task(right_vertex_map);
                for (uint32_t i = 0; i < pool->get_thread_number(); ++i) {
                    pool->submit_task([=] {
                        for (auto iter = *location_vector->at(i); iter != *location_vector->at(i + 1); ++iter) {
                            auto [r, r_vertex] = *iter;
                            auto degree = r_vertex->get_degree();
                            degree_mutex_vector->at(degree)->lock();
                            if (!degree_vector->at(degree)) {
                                degree_vector->at(degree) = make_shared<set<uint32_t>>();
                            }
                            degree_vector->at(degree)->insert(r);
                            degree_mutex_vector->at(degree)->unlock();
                        }
                    });
                }
            }
            pool->barrier();
        }
        uint32_t i = 0;
        vertex_priority_map->reserve(G->get_left_vertex_number() + G->get_right_vertex_number());
        for (auto iter1 = degree_vector->rbegin(); iter1 != degree_vector->rend(); ++iter1) {
            if (*iter1) {
                auto v_set = *iter1;
                for (auto iter2 = v_set->rbegin(); iter2 != v_set->rend(); ++iter2) {
                    vertex_priority_map->insert({*iter2, ++i});
                }
            }
        }
    }

    void quasi_wing_maintenance::previous_k_max(
            const shared_ptr<unordered_map<shared_ptr<abstract_bipartite_edge>, uint32_t>> &edge_wing_map,
            uint32_t &previous_k_max,
            const shared_ptr<thread_pool>& pool) {
        auto thread_number = pool->get_thread_number();
        auto global_mutex = make_shared<mutex>();
        auto k_max = make_shared<uint32_t>(edge_wing_map->begin()->second);
        auto location_vector = pool->split_task(edge_wing_map);
        for(uint32_t i = 0; i < thread_number; ++i){
            pool->submit_task([=]{
                auto &sub_begin = *location_vector->at(i);
                auto &sub_end = *location_vector->at(i + 1);

                uint32_t sub_max = 0 ;

                for(auto iter = sub_begin; iter!=sub_end; ++iter){
                    auto &e = iter->first;
                    if(edge_wing_map->at(e) > sub_max)
                    {
                        sub_max = edge_wing_map->at(e);
                    }
                }

                global_mutex->lock();
                if(sub_max > *k_max){
                    *k_max = sub_max;
                }
                global_mutex->unlock();
            });
        }
        pool->barrier();
        previous_k_max = *k_max;
    }

    void quasi_wing_maintenance::wl_move(const shared_ptr<abstract_bipartite_edge> &e,
                                         const shared_ptr<unordered_map<uint32_t,shared_ptr<map<uint32_t,shared_ptr<unordered_set<uint32_t>>>>>> &WL,
                                         uint32_t source,
                                         uint32_t destination){

        auto l = e->get_left_vertex_id();
        auto r = e->get_right_vertex_id();
        WL->at(l)->at(source)->erase(r);
        if(WL->at(l)->at(source)->empty()){
            WL->at(l)->erase(source);
        }
        WL->at(r)->at(source)->erase(l);
        if(WL->at(r)->at(source)->empty()){
            WL->at(r)->erase(source);
        }

        if(!WL->at(l)->count(destination)){
            WL->at(l)->insert({destination, make_shared<unordered_set<uint32_t>>()});
        }
        if(!WL->at(r)->count(destination)){
            WL->at(r)->insert({destination, make_shared<unordered_set<uint32_t>>()});
        }
        WL->at(l)->at(destination)->insert(r);
        WL->at(r)->at(destination)->insert(l);
    }
}





