
#include "bipartite_core/branch_bipartite_core_decomposition.h"

namespace scnu {

    void branch_bipartite_core_decomposition::init(const shared_ptr<scnu::abstract_bipartite_graph> &B,
                                                   const shared_ptr<unordered_map<uint32_t, shared_ptr<std::mutex>>> &left_mutex_map,
                                                   const shared_ptr<unordered_map<uint32_t, shared_ptr<std::mutex>>> &right_mutex_map,
                                                   const shared_ptr<unordered_map<uint32_t, shared_ptr<scnu::bipartite_core_left_store_index>>> &left_index_map,
                                                   const shared_ptr<unordered_map<uint32_t, shared_ptr<scnu::bipartite_core_right_store_index>>> &right_index_map,
                                                   const shared_ptr<::thread_pool> &pool) {
        for (const auto &[l, l_vertex]: *B->get_left_vertex_map()) {
            left_mutex_map->insert({l, shared_ptr<mutex>()});
            left_index_map->insert({l, shared_ptr<bipartite_core_left_store_index>()});
        }

        for (const auto &[r, r_vertex]: *B->get_right_vertex_map()) {
            right_mutex_map->insert({r, shared_ptr<mutex>()});
            right_index_map->insert({r, shared_ptr<bipartite_core_right_store_index>()});
        }

        {
            auto location_vector = pool->split_task(left_mutex_map);
            for (uint32_t i = 0; i < pool->get_thread_number(); ++i) {
                pool->submit_task([=] {
                    for (auto iter = *location_vector->at(i); iter != *location_vector->at(i + 1); ++iter) {
                        auto &l = iter->first;
                        left_mutex_map->at(l) = make_shared<mutex>();
                        left_index_map->at(l) = make_shared<bipartite_core_left_store_index>();
                    }
                });
            }
        }

        {
            auto location_vector = pool->split_task(right_mutex_map);
            for (uint32_t i = 0; i < pool->get_thread_number(); ++i) {
                pool->submit_task([=] {
                    for (auto iter = *location_vector->at(i); iter != *location_vector->at(i + 1); ++iter) {
                        auto &r = iter->first;
                        right_mutex_map->at(r) = make_shared<mutex>();
                        right_index_map->at(r) = make_shared<bipartite_core_right_store_index>();
                    }
                });
            }
        }
        pool->barrier();
    }

    uint32_t branch_bipartite_core_decomposition::decompose(const shared_ptr<abstract_bipartite_graph> &B,
                                                            const shared_ptr<unordered_map<uint32_t, shared_ptr<bipartite_core_left_store_index>>> &left_index_map,
                                                            const shared_ptr<unordered_map<uint32_t, shared_ptr<bipartite_core_right_store_index>>> &right_index_map) {
        auto left_map = make_shared<unordered_map<uint32_t, uint32_t>>();
        auto right_map = make_shared<unordered_map<uint32_t, uint32_t>>();

        for (const auto &[l, l_vertex]: *B->get_left_vertex_map()) {
            auto degree = l_vertex->get_degree();
            left_map->insert({l, degree});

            left_index_map->insert({l, make_shared<bipartite_core_left_store_index>()});
            for(uint32_t i = 1; i<= degree; ++i){
                left_index_map->at(l)->insert(i, 1);
            }
        }

        for (const auto&[r, r_vertex]: *B->get_right_vertex_map()) {
            auto degree = r_vertex->get_degree();
            right_map->insert({r, degree});

            right_index_map->insert({r, make_shared<bipartite_core_right_store_index>()});
            for(uint32_t j = 1; j <= degree; ++j){
                right_index_map->at(r)->insert(j, 1);
            }
        }

        for (const auto&[l, l_vertex]: *B->get_left_vertex_map()) {
            auto degree = l_vertex->get_degree();
            for(const auto &[r, e]:*l_vertex->get_edge_map()){
                right_index_map->at(r)->insert(1, degree);
            }
        }

        for (const auto&[r, r_vertex]: *B->get_right_vertex_map()) {
            auto degree = r_vertex->get_degree();
            for(const auto &[l, e]:*r_vertex->get_edge_map()){
                left_index_map->at(l)->insert(1, degree);
            }
        }


        uint32_t max_k = 1;

        for (uint32_t k = 2; true; ++k) {
            {
                find_core(B, left_map, right_map, k);
                if(!left_map->empty() && !right_map->empty()) {
                    max_k = k;
                }else{
                    break;
                }
            }

            {
                auto sub_left_degree_map = container_copy::to_unordered_map<uint32_t, uint32_t>(left_map);
                auto sub_right_degree_map = container_copy::to_unordered_map<uint32_t, uint32_t>(right_map);
                find_left_core(B, sub_left_degree_map, sub_right_degree_map,
                                   left_index_map,
                                   right_index_map, k);
            }

            {
                auto sub_left_degree_map = container_copy::to_unordered_map<uint32_t, uint32_t>(left_map);
                auto sub_right_degree_map = container_copy::to_unordered_map<uint32_t, uint32_t>(right_map);
                find_right_core(B, sub_left_degree_map, sub_right_degree_map, left_index_map,
                                    right_index_map,
                                    k);
            }
        }

        return max_k;
    }

    uint32_t branch_bipartite_core_decomposition::decompose(const shared_ptr<abstract_bipartite_graph> &B,
                                                            const shared_ptr<unordered_map<uint32_t, shared_ptr<mutex>>> &left_mutex_map,
                                                            const shared_ptr<unordered_map<uint32_t, shared_ptr<mutex>>> &right_mutex_map,
                                                            const shared_ptr<unordered_map<uint32_t, shared_ptr<bipartite_core_left_store_index>>> &left_index_map,
                                                            const shared_ptr<unordered_map<uint32_t, shared_ptr<bipartite_core_right_store_index>>> &right_index_map,
                                                            const shared_ptr<thread_pool> &pool) {
        auto left_map = make_shared<unordered_map<uint32_t, uint32_t>>();
        auto right_map = make_shared<unordered_map<uint32_t, uint32_t>>();

        for (const auto&[l, l_vertex]: *B->get_left_vertex_map()) {
            auto degree = l_vertex->get_degree();
            left_map->insert({l, degree});
        }

        for (const auto&[r, r_vertex]: *B->get_right_vertex_map()) {
            auto degree = r_vertex->get_degree();
            right_map->insert({r, degree});
        }

        for (const auto&[l, l_vertex]: *B->get_left_vertex_map()) {
            pool->submit_task([=]{
                auto degree = l_vertex->get_degree();
                for(uint32_t i = 1; i<= degree; ++i){
                    left_index_map->at(l)->insert(i, 1);
                }
            });
        }

        for (const auto&[r, r_vertex]: *B->get_right_vertex_map()) {
            pool->submit_task([=]{
                auto degree =  r_vertex->get_degree();
                for(uint32_t j = 1; j <= degree; ++j){
                    right_index_map->at(r)->insert(j, 1);
                }
            });
        }
        pool->barrier();
        for (const auto&[l, l_vertex]: *B->get_left_vertex_map()) {
            pool->submit_task([=]{
                auto degree = l_vertex->get_degree();
                for(const auto &[r, e]:*l_vertex->get_edge_map()){
                    right_mutex_map->at(r)->lock();
                    right_index_map->at(r)->insert(1, degree);
                    right_mutex_map->at(r)->unlock();
                }
            });
        }
        for (const auto&[r, r_vertex]: *B->get_right_vertex_map()) {
            pool->submit_task([=]{
                auto degree =  r_vertex->get_degree();
                for(const auto &[l, e]:*r_vertex->get_edge_map()){
                    left_mutex_map->at(l)->lock();
                    left_index_map->at(l)->insert(1, degree);
                    left_mutex_map->at(l)->unlock();
                }
            });
        }

        uint32_t max_k = 1;

        for (uint32_t k = 2; true; ++k) {
            {
                find_core(B, left_map, right_map, k);
                if(!left_map->empty() && !right_map->empty()) {
                    max_k = k;
                }else{
                    break;
                }
            }

            {
                auto sub_left_degree_map = container_copy::to_unordered_map<uint32_t, uint32_t>(left_map);
                auto sub_right_degree_map = container_copy::to_unordered_map<uint32_t, uint32_t>(right_map);
                pool->submit_task([=] {
                    find_left_core(B, left_mutex_map, right_mutex_map, sub_left_degree_map, sub_right_degree_map,
                                   left_index_map,
                                   right_index_map, k);
                });

            }

            {
                auto sub_left_degree_map = container_copy::to_unordered_map<uint32_t, uint32_t>(left_map);
                auto sub_right_degree_map = container_copy::to_unordered_map<uint32_t, uint32_t>(right_map);
                pool->submit_task([=] {
                    find_right_core(B, left_mutex_map, right_mutex_map, sub_left_degree_map, sub_right_degree_map, left_index_map,
                                    right_index_map,
                                    k);
                });
            }

        }
        pool->barrier();

        return max_k;
    }

    uint32_t branch_bipartite_core_decomposition::decompose(const shared_ptr<abstract_bipartite_graph> &B,
                                                            const shared_ptr<unordered_map<uint32_t, shared_ptr<mutex>>> &left_mutex_map,
                                                            const shared_ptr<unordered_map<uint32_t, shared_ptr<mutex>>> &right_mutex_map,
                                                            const shared_ptr<unordered_map<uint32_t, shared_ptr<bipartite_core_left_store_index>>> &left_index_map,
                                                            const shared_ptr<unordered_map<uint32_t, shared_ptr<bipartite_core_right_store_index>>> &right_index_map,
                                                            const shared_ptr<bipartite_core_order_index> &core_order_index,
                                                            const shared_ptr<bipartite_core_rem_degree_index> &core_rem_degree_index,
                                                            const shared_ptr<bipartite_core_degree_index> &core_degree_index,
                                                            const shared_ptr<thread_pool> &pool) {
        auto l_set = make_shared<unordered_set<uint32_t>>();
        auto r_set = make_shared<unordered_set<uint32_t>>();


        core_degree_index->get_middle_map()->reserve(B->get_left_vertex_number()+B->get_right_vertex_number());
        core_rem_degree_index->get_middle_map()->reserve(B->get_left_vertex_number()+B->get_right_vertex_number());

        auto vertex_degree_map = make_shared<unordered_map<uint32_t, uint32_t>>(B->get_left_vertex_number()+B->get_right_vertex_number());

        for (const auto&[l, l_vertex]: *B->get_left_vertex_map()) {
            l_set->insert(l);

            auto degree = l_vertex->get_degree();
            vertex_degree_map->insert({l, degree});
            core_degree_index->get_middle_map()->insert({l, 0});
            core_rem_degree_index->get_middle_map()->insert({l, 0});
        }

        for (const auto&[r, r_vertex]: *B->get_right_vertex_map()) {
            r_set->insert(r);

            auto degree = r_vertex->get_degree();
            vertex_degree_map->insert({r, degree});

            core_degree_index->get_middle_map()->insert({r, 0});
            core_rem_degree_index->get_middle_map()->insert({r, 0});
        }

        for (const auto&[l, l_vertex]: *B->get_left_vertex_map()) {
            pool->submit_task([=]{
                left_index_map->at(l) =  make_shared<bipartite_core_left_store_index>();

                auto degree = l_vertex->get_degree();
                for(uint32_t i = 1; i <= degree; ++i){
                    left_index_map->at(l)->insert(i, 1);
                }

                uint32_t neighbor_degree = 0;
                for(const auto &[r, e]:*l_vertex->get_edge_map()){
                    auto r_degree = B->get_right_vertex(r)->get_degree();
                    if(r_degree > neighbor_degree){
                        neighbor_degree = r_degree;
                    }
                }
                left_index_map->at(l)->insert(1, neighbor_degree);
            });
        }

        for (const auto&[r, r_vertex]: *B->get_right_vertex_map()) {
            pool->submit_task([=]{
                right_index_map->at(r) = make_shared<bipartite_core_right_store_index>();

                auto degree =  r_vertex->get_degree();
                for(uint32_t j = 1; j <= degree; ++j){
                    right_index_map->at(r)->insert(j, 1);
                }

                uint32_t neighbor_degree = 0;
                for (const auto &[l, e]: *r_vertex->get_edge_map()) {
                    auto l_degree = B->get_left_vertex(l)->get_degree();
                    if (l_degree > neighbor_degree) {
                        neighbor_degree = l_degree;
                    }
                }
                right_index_map->at(r)->insert(1, neighbor_degree);
            });
        }
        pool->barrier();

        uint32_t max_k = B->get_left_vertex_number() > 0 ? 1 : 0;
        for (uint32_t k = 2; !vertex_degree_map->empty(); ++k) {
            /**
             * @brief middle path
             */
            {
                auto k_order_map = core_order_index->get_middle_map();
                auto k_rem_degree_map = core_rem_degree_index->get_middle_map();
                auto k_core_degree_map = core_degree_index->get_middle_map();

                find_core(B, l_set, r_set, vertex_degree_map, k_order_map, k_rem_degree_map,
                          k_core_degree_map, k);
            }

            if (vertex_degree_map->empty()) {
                break;
            } else {
                max_k = k;

                core_order_index->insert_left(k,
                                              make_shared<unordered_map<uint32_t, shared_ptr<extend_list<int, uint32_t>>>>());
                core_rem_degree_index->insert_left(k, make_shared<unordered_map<uint32_t, uint32_t>>(
                        l_set->size() + r_set->size()));
                core_degree_index->insert_left(k, make_shared<unordered_map<uint32_t, uint32_t>>(
                        l_set->size() + r_set->size()));

                core_order_index->insert_right(k,
                                               make_shared<unordered_map<uint32_t, shared_ptr<extend_list<int, uint32_t>>>>());
                core_rem_degree_index->insert_right(k, make_shared<unordered_map<uint32_t, uint32_t>>(
                        l_set->size() + r_set->size()));
                core_degree_index->insert_right(k, make_shared<unordered_map<uint32_t, uint32_t>>(
                        l_set->size() + r_set->size()));

                {
                    auto left_k_order_map = core_order_index->get_left_map(k);
                    auto left_k_rem_degree_map = core_rem_degree_index->get_left_map(k);
                    auto left_k_core_degree_map = core_degree_index->get_left_map(k);

                    auto right_k_order_map = core_order_index->get_right_map(k);
                    auto right_k_rem_degree_map = core_rem_degree_index->get_right_map(k);
                    auto right_k_core_degree_map = core_degree_index->get_right_map(k);


                    auto left_sub_vertex_degree_map = container_copy::to_unordered_map<uint32_t, uint32_t>(
                            vertex_degree_map);

                    auto right_sub_vertex_degree_map = container_copy::to_unordered_map<uint32_t, uint32_t>(
                            vertex_degree_map);

                    pool->submit_task([=] {
                        /**
                         * @brief left path
                         */
                        {
                            auto sub_l_set = make_shared<unordered_set<uint32_t>>();
                            auto sub_r_set = make_shared<unordered_set<uint32_t>>();

                            for (const auto &[v, degree]: *left_sub_vertex_degree_map) {
                                left_k_core_degree_map->insert({v, degree});
                                left_k_rem_degree_map->insert({v, 0});
                                if (B->get_left_vertex(v)) {
                                    sub_l_set->insert(v);
                                } else {
                                    sub_r_set->insert(v);
                                }
                            }

                            find_left_core(B, left_mutex_map, right_mutex_map, sub_l_set, sub_r_set,
                                           left_sub_vertex_degree_map,
                                           left_index_map,
                                           right_index_map, left_k_order_map, left_k_rem_degree_map,
                                           left_k_core_degree_map, k);
                        }


                        /**
                         * @brief right path
                         */
                        {
                            auto sub_l_set = make_shared<unordered_set<uint32_t>>();
                            auto sub_r_set = make_shared<unordered_set<uint32_t>>();

                            for (const auto &[v, degree]: *right_sub_vertex_degree_map) {
                                right_k_core_degree_map->insert({v, degree});
                                right_k_rem_degree_map->insert({v, 0});
                                if (B->get_left_vertex(v)) {
                                    sub_l_set->insert(v);
                                } else {
                                    sub_r_set->insert(v);
                                }
                            }

                            find_right_core(B, left_mutex_map, right_mutex_map, sub_l_set, sub_r_set,
                                            right_sub_vertex_degree_map,
                                            left_index_map,
                                            right_index_map, right_k_order_map, right_k_rem_degree_map,
                                            right_k_core_degree_map,
                                            k);
                        }

                    });
                }
            }
        }

        pool->barrier();
        return max_k;
    }

    uint32_t branch_bipartite_core_decomposition::decompose2(const shared_ptr<abstract_bipartite_graph> &B,
                                                             const shared_ptr<unordered_map<uint32_t, shared_ptr<mutex>>> &left_mutex_map,
                                                             const shared_ptr<unordered_map<uint32_t, shared_ptr<mutex>>> &right_mutex_map,
                                                             const shared_ptr<unordered_map<uint32_t, shared_ptr<bipartite_core_left_store_index>>> &left_index_map,
                                                             const shared_ptr<unordered_map<uint32_t, shared_ptr<bipartite_core_right_store_index>>> &right_index_map,
                                                             const shared_ptr<bipartite_core_order_index> &core_order_index,
                                                             const shared_ptr<bipartite_core_rem_degree_index> &core_rem_degree_index,
                                                             const shared_ptr<bipartite_core_degree_index> &core_degree_index,
                                                             const shared_ptr<thread_pool> &pool) {
        auto l_set = make_shared<unordered_set<uint32_t>>();
        auto r_set = make_shared<unordered_set<uint32_t>>();

        core_degree_index->get_middle_map()->reserve(B->get_left_vertex_number() + B->get_right_vertex_number());
        core_rem_degree_index->get_middle_map()->reserve(B->get_left_vertex_number() + B->get_right_vertex_number());

        auto vertex_degree_map = make_shared<unordered_map<uint32_t, uint32_t>>(
                B->get_left_vertex_number() + B->get_right_vertex_number());

        for (const auto &[l, l_vertex]: *B->get_left_vertex_map()) {
            l_set->insert(l);

            auto degree = l_vertex->get_degree();
            vertex_degree_map->insert({l, degree});
            core_degree_index->get_middle_map()->insert({l, 0});
            core_rem_degree_index->get_middle_map()->insert({l, 0});
        }

        for (const auto &[r, r_vertex]: *B->get_right_vertex_map()) {
            r_set->insert(r);

            auto degree = r_vertex->get_degree();
            vertex_degree_map->insert({r, degree});

            core_degree_index->get_middle_map()->insert({r, 0});
            core_rem_degree_index->get_middle_map()->insert({r, 0});
        }

        for (const auto &[l, l_vertex]: *B->get_left_vertex_map()) {
            pool->submit_task([=] {
                left_index_map->at(l) = make_shared<bipartite_core_left_store_index>();

                auto degree = l_vertex->get_degree();
                for (uint32_t i = 1; i <= degree; ++i) {
                    left_index_map->at(l)->insert(i, 1);
                }

                uint32_t neighbor_degree = 0;
                for (const auto &[r, e]: *l_vertex->get_edge_map()) {
                    auto r_degree = B->get_right_vertex(r)->get_degree();
                    if (r_degree > neighbor_degree) {
                        neighbor_degree = r_degree;
                    }
                }
                left_index_map->at(l)->insert(1, neighbor_degree);
            });
        }

        for (const auto &[r, r_vertex]: *B->get_right_vertex_map()) {
            pool->submit_task([=] {
                right_index_map->at(r) = make_shared<bipartite_core_right_store_index>();

                auto degree = r_vertex->get_degree();
                for (uint32_t j = 1; j <= degree; ++j) {
                    right_index_map->at(r)->insert(j, 1);
                }

                uint32_t neighbor_degree = 0;
                for (const auto &[l, e]: *r_vertex->get_edge_map()) {
                    auto l_degree = B->get_left_vertex(l)->get_degree();
                    if (l_degree > neighbor_degree) {
                        neighbor_degree = l_degree;
                    }
                }
                right_index_map->at(r)->insert(1, neighbor_degree);
            });
        }
        pool->barrier();

        auto k_map = make_shared<map<uint32_t, shared_ptr<unordered_map<uint32_t, uint32_t>>>>();

        for (uint32_t k = 2; !vertex_degree_map->empty(); ++k) {
            {
                auto k_order_map = core_order_index->get_middle_map();
                auto k_rem_degree_map = core_rem_degree_index->get_middle_map();
                auto k_core_degree_map = core_degree_index->get_middle_map();

                find_core(B, l_set, r_set, vertex_degree_map, k_order_map, k_rem_degree_map,
                          k_core_degree_map, k);

                if (vertex_degree_map->empty()) {
                    break;
                }

                k_map->insert({k, make_shared<unordered_map<uint32_t, uint32_t>>(*vertex_degree_map)});

                core_order_index->insert_left(k,
                                              make_shared<unordered_map<uint32_t, shared_ptr<extend_list<int, uint32_t>>>>());
                core_rem_degree_index->insert_left(k, make_shared<unordered_map<uint32_t, uint32_t>>(
                        l_set->size() + r_set->size()));
                core_degree_index->insert_left(k, make_shared<unordered_map<uint32_t, uint32_t>>(
                        l_set->size() + r_set->size()));

                core_order_index->insert_right(k,
                                               make_shared<unordered_map<uint32_t, shared_ptr<extend_list<int, uint32_t>>>>());
                core_rem_degree_index->insert_right(k, make_shared<unordered_map<uint32_t, uint32_t>>(
                        l_set->size() + r_set->size()));
                core_degree_index->insert_right(k, make_shared<unordered_map<uint32_t, uint32_t>>(
                        l_set->size() + r_set->size()));
            }
        }

        /**
         * @brief left path
         */
        for (const auto &[k, k_vertex_degree_map]: *k_map) {
            pool->submit_task([=] {
                auto k_order_map = core_order_index->get_left_map(k);
                auto k_rem_degree_map = core_rem_degree_index->get_left_map(k);
                auto k_core_degree_map = core_degree_index->get_left_map(k);


                auto sub_vertex_degree_map = container_copy::to_unordered_map<uint32_t, uint32_t>(
                        k_vertex_degree_map);

                auto sub_l_set = make_shared<unordered_set<uint32_t>>();
                auto sub_r_set = make_shared<unordered_set<uint32_t>>();

                for (const auto &[v, degree]: *sub_vertex_degree_map) {
                    k_core_degree_map->insert({v, degree});
                    k_rem_degree_map->insert({v, 0});
                    if (B->get_left_vertex(v)) {
                        sub_l_set->insert(v);
                    } else {
                        sub_r_set->insert(v);
                    }
                }

                find_left_core(B, left_mutex_map, right_mutex_map, sub_l_set, sub_r_set, sub_vertex_degree_map,
                               left_index_map,
                               right_index_map, k_order_map, k_rem_degree_map, k_core_degree_map, k);
            });
        }

        /**
         * @brief right path
         */
        for (const auto &[k, k_vertex_degree_map]: *k_map) {
            pool->submit_task([=] {
                auto k_order_map = core_order_index->get_right_map(k);
                auto k_rem_degree_map = core_rem_degree_index->get_right_map(k);
                auto k_core_degree_map = core_degree_index->get_right_map(k);

                auto sub_vertex_degree_map = container_copy::to_unordered_map<uint32_t, uint32_t>(k_vertex_degree_map);

                auto sub_l_set = make_shared<unordered_set<uint32_t>>();
                auto sub_r_set = make_shared<unordered_set<uint32_t>>();

                for (const auto &[v, degree]: *sub_vertex_degree_map) {
                    k_core_degree_map->insert({v, degree});
                    k_rem_degree_map->insert({v, 0});
                    if (B->get_left_vertex(v)) {
                        sub_l_set->insert(v);
                    } else {
                        sub_r_set->insert(v);
                    }
                }

                find_right_core(B, left_mutex_map, right_mutex_map, sub_l_set, sub_r_set, sub_vertex_degree_map,
                                left_index_map,
                                right_index_map, k_order_map, k_rem_degree_map, k_core_degree_map,
                                k);
            });
        }
        pool->barrier();

        return k_map->empty() ? 0 : k_map->rbegin()->first;
    }

    uint32_t branch_bipartite_core_decomposition::decompose(const shared_ptr<abstract_bipartite_graph> &B,
                                                            const shared_ptr<unordered_map<uint32_t, shared_ptr<mutex>>> &left_mutex_map,
                                                            const shared_ptr<unordered_map<uint32_t, shared_ptr<mutex>>> &right_mutex_map,
                                                            const shared_ptr<unordered_map<uint32_t, shared_ptr<bipartite_core_branch_store_index>>> &left_index_map,
                                                            const shared_ptr<unordered_map<uint32_t, shared_ptr<bipartite_core_branch_store_index>>> &right_index_map,
                                                            const shared_ptr<thread_pool> &pool) {
        auto left_map = make_shared<unordered_map<uint32_t, uint32_t>>();
        auto right_map = make_shared<unordered_map<uint32_t, uint32_t>>();

        for (const auto &[l, l_vertex]: *B->get_left_vertex_map()) {
            auto degree = l_vertex->get_degree();
            left_map->insert({l, degree});
        }

        for (const auto &[r, r_vertex]: *B->get_right_vertex_map()) {
            auto degree = r_vertex->get_degree();
            right_map->insert({r, degree});
        }

        for (const auto &[l, l_vertex]: *B->get_left_vertex_map()) {
            pool->submit_task([=] {
                auto degree = l_vertex->get_degree();
                for (uint32_t i = 1; i <= degree; ++i) {
                    left_index_map->at(l)->insert(i, 1);
                }
            });
        }

        for (const auto &[r, r_vertex]: *B->get_right_vertex_map()) {
            pool->submit_task([=] {
                auto degree = r_vertex->get_degree();
                for (uint32_t j = 1; j <= degree; ++j) {
                    right_index_map->at(r)->insert(j, 1);
                }
            });
        }
        pool->barrier();
        for (const auto &[l, l_vertex]: *B->get_left_vertex_map()) {
            pool->submit_task([=] {
                auto degree = l_vertex->get_degree();
                for (const auto &[r, e]: *l_vertex->get_edge_map()) {
                    right_mutex_map->at(r)->lock();
                    right_index_map->at(r)->insert(1, degree);
                    right_mutex_map->at(r)->unlock();
                }
            });
        }
        for (const auto &[r, r_vertex]: *B->get_right_vertex_map()) {
            pool->submit_task([=] {
                auto degree = r_vertex->get_degree();
                for (const auto &[l, e]: *r_vertex->get_edge_map()) {
                    left_mutex_map->at(l)->lock();
                    left_index_map->at(l)->insert(1, degree);
                    left_mutex_map->at(l)->unlock();
                }
            });
        }

        uint32_t max_k = 1;

        for (uint32_t k = 2; true; ++k) {
            {
                find_core(B, left_map, right_map, k);
                if (!left_map->empty() && !right_map->empty()) {
                    max_k = k;
                } else {
                    break;
                }
            }

            {
                auto sub_left_degree_map = container_copy::to_unordered_map<uint32_t, uint32_t>(left_map);
                auto sub_right_degree_map = container_copy::to_unordered_map<uint32_t, uint32_t>(right_map);
                pool->submit_task([=] {
                    find_left_core(B, left_mutex_map, right_mutex_map, sub_left_degree_map, sub_right_degree_map,
                                   left_index_map,
                                   right_index_map, k);
                });

            }

            {
                auto sub_left_degree_map = container_copy::to_unordered_map<uint32_t, uint32_t>(left_map);
                auto sub_right_degree_map = container_copy::to_unordered_map<uint32_t, uint32_t>(right_map);
                pool->submit_task([=] {
                    find_right_core(B, left_mutex_map, right_mutex_map, sub_left_degree_map, sub_right_degree_map,
                                    left_index_map,
                                    right_index_map,
                                    k);
                });
            }

        }
        pool->barrier();

        return max_k;
    }

    void branch_bipartite_core_decomposition::find_core(const shared_ptr<abstract_bipartite_graph> &B,
                                                        const shared_ptr<unordered_map<uint32_t, uint32_t>> &left_degree_map,
                                                        const shared_ptr<unordered_map<uint32_t, uint32_t>> &right_degree_map,
                                                        uint32_t k) {
        auto evicted_l_set = make_shared<unordered_set<uint32_t>>();
        auto evicted_r_set = make_shared<unordered_set<uint32_t>>();


        for (const auto &[l, l_degree]: *left_degree_map) {
            if (l_degree < k) {
                evicted_l_set->emplace(l);
            }
        }

        for (const auto &[r, r_degree]: *right_degree_map) {
            if (r_degree < k) {
                evicted_r_set->emplace(r);
            }
        }


        while (!evicted_l_set->empty() || !evicted_r_set->empty()) {
            while(!evicted_l_set->empty()){
                auto l = *evicted_l_set->begin();
                evicted_l_set->erase(l);

                left_degree_map->erase(l);

                for (const auto&[r, e]: *B->get_left_vertex(l)->get_edge_map()) {
                    if (right_degree_map->count(r) && right_degree_map->at(r) >= k) {
                        --right_degree_map->at(r);
                        if (right_degree_map->at(r) < k) {
                            evicted_r_set->insert(r);
                        }
                    }
                }
            }


            while (!evicted_r_set->empty()) {
                auto r = *evicted_r_set->begin();
                evicted_r_set->erase(r);

                right_degree_map->erase(r);

                for (const auto &[l, e]: *B->get_right_vertex(r)->get_edge_map()) {
                    if (left_degree_map->count(l) && left_degree_map->at(l) >= k) {
                        --left_degree_map->at(l);
                        if (left_degree_map->at(l) < k) {
                            evicted_l_set->insert(l);
                        }
                    }
                }
            }
        }
    }

    void branch_bipartite_core_decomposition::find_core(const shared_ptr<abstract_bipartite_graph> &B,
                                                        const shared_ptr<unordered_set<uint32_t>> &l_set,
                                                        const shared_ptr<unordered_set<uint32_t>> &r_set,
                                                        const shared_ptr<unordered_map<uint32_t, uint32_t>> &vertex_degree_map,
                                                        const shared_ptr<unordered_map<uint32_t, shared_ptr<extend_list<int, uint32_t>>>>& middle_order_map,
                                                        const shared_ptr<unordered_map<uint32_t, uint32_t>>& middle_rem_degree_map,
                                                        const shared_ptr<unordered_map<uint32_t, uint32_t>>& middle_core_degree_map,
                                                        uint32_t k) {
        auto evicted_l_set = make_shared<unordered_set<uint32_t>>();
        auto evicted_r_set = make_shared<unordered_set<uint32_t>>();

        for (const auto &l: *l_set) {
            middle_core_degree_map->at(l) = vertex_degree_map->at(l);
            if (vertex_degree_map->at(l) < k) {
                evicted_l_set->insert(l);
            }
        }

        for (const auto &r: *r_set) {
            middle_core_degree_map->at(r) = vertex_degree_map->at(r);
            if (vertex_degree_map->at(r) < k) {
                evicted_r_set->insert(r);
            }
        }

        if(!middle_order_map->count(k - 1)){
            middle_order_map->insert({k - 1, make_shared<extend_list<int, uint32_t>>()});
        }

        while (!evicted_l_set->empty() || !evicted_r_set->empty()) {
            while(!evicted_l_set->empty()){
                auto l = *evicted_l_set->begin();
                evicted_l_set->erase(l);

                middle_rem_degree_map->at(l) = vertex_degree_map->at(l);
                middle_order_map->at(k - 1)->push_back(l);

                l_set->erase(l);
                vertex_degree_map->erase(l);

                for (const auto&[r, e]: *B->get_left_vertex(l)->get_edge_map()) {
                    if (vertex_degree_map->count(r)) {
                        --vertex_degree_map->at(r);
                        if (vertex_degree_map->at(r) < k) {
                            evicted_r_set->insert(r);
                        }
                    }
                }
            }


            while (!evicted_r_set->empty()) {
                auto r = *evicted_r_set->begin();
                evicted_r_set->erase(r);

                middle_order_map->at(k - 1)->push_back(r);
                middle_rem_degree_map->at(r) = vertex_degree_map->at(r);

                r_set->erase(r);
                vertex_degree_map->erase(r);

                for (const auto &[l, e]: *B->get_right_vertex(r)->get_edge_map()) {
                    if (vertex_degree_map->count(l)) {
                        --vertex_degree_map->at(l);
                        if (vertex_degree_map->at(l) < k) {
                            evicted_l_set->insert(l);
                        }
                    }
                }
            }
        }
    }

    uint32_t branch_bipartite_core_decomposition::find_left_core(const shared_ptr<abstract_bipartite_graph> &B,
                                                                 const shared_ptr<unordered_map<uint32_t, uint32_t>> &left_map,
                                                                 const shared_ptr<unordered_map<uint32_t, uint32_t>> &right_map,
                                                                 const shared_ptr<unordered_map<uint32_t, shared_ptr<bipartite_core_left_store_index>>>& left_index_map,
                                                                 const shared_ptr<unordered_map<uint32_t, shared_ptr<bipartite_core_right_store_index>>>& right_index_map,
                                                                 uint32_t k) {
        auto evicted_l_set = make_shared<unordered_set<uint32_t>>();

        uint32_t max_i = k;
        while (!left_map->empty()) {
            uint32_t  i = left_map->begin()->second;

            for (const auto&[l, l_degree]: *left_map) {
                if (l_degree < i) {
                    i = l_degree;

                    evicted_l_set->clear();
                    evicted_l_set->insert(l);
                } else if(l_degree == i)
                {
                    evicted_l_set->insert(l);
                }
            }

            max_i = i;

            while (!evicted_l_set->empty()) {
                auto l = *evicted_l_set->begin();
                evicted_l_set->erase(l);

                left_map->erase(l);

                for (uint32_t index = k; index <= i; ++index) {
                    left_index_map->at(l)->insert(index, k);
                }

                for (const auto &[r, e]: *B->get_left_vertex(l)->get_edge_map()) {
                    if (right_map->count(r) && right_map->at(r) >= k) {
                        --right_map->at(r);
                        if (right_map->at(r) < k) {
                            right_map->erase(r);
                            right_index_map->at(r)->insert(k, i);

                            for (const auto &[l2, e2]: *B->get_right_vertex(r)->get_edge_map()) {
                                if (left_map->count(l2) && left_map->at(l2) > i) {
                                    --left_map->at(l2);
                                    if (left_map->at(l2) == i) {
                                        evicted_l_set->emplace(l2);
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
        return max_i;
    }


    uint32_t branch_bipartite_core_decomposition::find_left_core(const shared_ptr<abstract_bipartite_graph> &B,
                                                                 const shared_ptr<unordered_map<uint32_t, shared_ptr<mutex>>>& left_mutex_map,
                                                                 const shared_ptr<unordered_map<uint32_t, shared_ptr<mutex>>>& right_mutex_map,
                                                                 const shared_ptr<unordered_map<uint32_t, uint32_t>> &left_map,
                                                                 const shared_ptr<unordered_map<uint32_t, uint32_t>> &right_map,
                                                                 const shared_ptr<unordered_map<uint32_t, shared_ptr<bipartite_core_left_store_index>>>& left_index_map,
                                                                 const shared_ptr<unordered_map<uint32_t, shared_ptr<bipartite_core_right_store_index>>>& right_index_map,
                                                                 uint32_t k) {
        auto evicted_l_set = make_shared<unordered_set<uint32_t>>();

        uint32_t max_i = k;
        while (!left_map->empty()) {
            uint32_t  i = left_map->begin()->second;

            for (const auto&[l, l_degree]: *left_map) {
                if (l_degree < i) {
                    i = l_degree;

                    evicted_l_set->clear();
                    evicted_l_set->insert(l);
                } else if(l_degree == i)
                {
                    evicted_l_set->insert(l);
                }
            }

            max_i = i;

            while (!evicted_l_set->empty()) {
                auto l = *evicted_l_set->begin();
                evicted_l_set->erase(l);

                left_map->erase(l);

                left_mutex_map->at(l)->lock();
                for (uint32_t index = k; index <= i; ++index) {
                    left_index_map->at(l)->insert(index, k);
                }
                left_mutex_map->at(l)->unlock();

                for (const auto &[r, e]: *B->get_left_vertex(l)->get_edge_map()) {
                    if (right_map->count(r) && right_map->at(r) >= k) {
                        --right_map->at(r);
                        if (right_map->at(r) < k) {
                            right_map->erase(r);

                            right_mutex_map->at(r)->lock();
                            right_index_map->at(r)->insert(k, i);
                            right_mutex_map->at(r)->unlock();

                            for (const auto &[l2, e2]: *B->get_right_vertex(r)->get_edge_map()) {
                                if (left_map->count(l2) && left_map->at(l2) > i) {
                                    --left_map->at(l2);
                                    if (left_map->at(l2) == i) {
                                        evicted_l_set->emplace(l2);
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
        return max_i;
    }

    uint32_t branch_bipartite_core_decomposition::find_left_core(const shared_ptr<abstract_bipartite_graph> &B,
                                                                 const shared_ptr<unordered_map<uint32_t, shared_ptr<mutex>>> &left_mutex_map,
                                                                 const shared_ptr<unordered_map<uint32_t, shared_ptr<mutex>>> &right_mutex_map,
                                                                 const shared_ptr<unordered_set<uint32_t>> &l_set,
                                                                 const shared_ptr<unordered_set<uint32_t>> &r_set,
                                                                 const shared_ptr<unordered_map<uint32_t, uint32_t>> &vertex_degree_map,
                                                                 const shared_ptr<unordered_map<uint32_t, shared_ptr<bipartite_core_left_store_index>>> &left_index_map,
                                                                 const shared_ptr<unordered_map<uint32_t, shared_ptr<bipartite_core_right_store_index>>> &right_index_map,
                                                                 const shared_ptr<unordered_map<uint32_t, shared_ptr<extend_list<int, uint32_t>>>> &k_oder_map,
                                                                 const shared_ptr<unordered_map<uint32_t, uint32_t>> &k_rem_degree_map,
                                                                 const shared_ptr<unordered_map<uint32_t, uint32_t>> &k_core_degree_map,
                                                                 uint32_t k) {
        auto evicted_l_set = make_shared<unordered_set<uint32_t>>();
        auto evicted_r_set = make_shared<unordered_set<uint32_t>>();

        uint32_t max_i = k;
        while (!l_set->empty()) {
            uint32_t i = UINT32_MAX;

            for (const auto &l: *l_set) {
                k_core_degree_map->at(l) = vertex_degree_map->at(l);

                if (vertex_degree_map->at(l) < i) {
                    i = vertex_degree_map->at(l);

                    evicted_l_set->clear();
                    evicted_l_set->insert(l);
                } else if (vertex_degree_map->at(l) == i) {
                    evicted_l_set->insert(l);
                }
            }

            for (const auto &r: *r_set) {
                k_core_degree_map->at(r) = vertex_degree_map->at(r);
            }

            if (!k_oder_map->count(i)) {
                k_oder_map->insert({i, make_shared<extend_list<int, uint32_t>>()});
            }

            max_i = i;
            while (!evicted_l_set->empty() || !evicted_r_set->empty()) {
                while (!evicted_l_set->empty()) {
                    auto l = *evicted_l_set->begin();
                    evicted_l_set->erase(l);

                    left_mutex_map->at(l)->lock();
                    for (uint32_t index = k; index <= i; ++index) {
                        left_index_map->at(l)->insert(index, k);
                    }
                    left_mutex_map->at(l)->unlock();

                    k_oder_map->at(i)->push_back(l);
                    k_rem_degree_map->at(l) = vertex_degree_map->at(l);
                    l_set->erase(l);
                    vertex_degree_map->erase(l);

                    for (const auto &[r, e]: *B->get_left_vertex(l)->get_edge_map()) {
                        if (vertex_degree_map->count(r)) {
                            --vertex_degree_map->at(r);
                            if (vertex_degree_map->at(r) < k) {
                                evicted_r_set->insert(r);
                            }
                        }
                    }
                }

                while (!evicted_r_set->empty()) {
                    auto r = *evicted_r_set->begin();
                    evicted_r_set->erase(r);

                    right_mutex_map->at(r)->lock();
                    right_index_map->at(r)->insert(k, i);
                    right_mutex_map->at(r)->unlock();

                    k_oder_map->at(i)->push_back(r);
                    k_rem_degree_map->at(r) = vertex_degree_map->at(r);

                    r_set->erase(r);
                    vertex_degree_map->erase(r);

                    for (const auto &[l, e]: *B->get_right_vertex(r)->get_edge_map()) {
                        if (vertex_degree_map->count(l)) {
                            --vertex_degree_map->at(l);
                            if (vertex_degree_map->at(l) == i) {
                                evicted_l_set->insert(l);
                            }
                        }
                    }
                }
            }
        }
        return max_i;
    }

    uint32_t branch_bipartite_core_decomposition::find_left_core(const shared_ptr<abstract_bipartite_graph> &B,
                                                                 const shared_ptr<unordered_map<uint32_t, shared_ptr<mutex>>> &left_mutex_map,
                                                                 const shared_ptr<unordered_map<uint32_t, shared_ptr<mutex>>> &right_mutex_map,
                                                                 const shared_ptr<unordered_map<uint32_t, uint32_t>> &left_map,
                                                                 const shared_ptr<unordered_map<uint32_t, uint32_t>> &right_map,
                                                                 const shared_ptr<unordered_map<uint32_t, shared_ptr<bipartite_core_branch_store_index>>> &left_index_map,
                                                                 const shared_ptr<unordered_map<uint32_t, shared_ptr<bipartite_core_branch_store_index>>> &right_index_map,
                                                                 uint32_t k) {
        auto evicted_l_set = make_shared<unordered_set<uint32_t>>();

        uint32_t max_i = k;
        while (!left_map->empty()) {
            uint32_t i = left_map->begin()->second;

            for (const auto &[l, l_degree]: *left_map) {
                if (l_degree < i) {
                    i = l_degree;

                    evicted_l_set->clear();
                    evicted_l_set->insert(l);
                } else if (l_degree == i) {
                    evicted_l_set->insert(l);
                }
            }

            max_i = i;

            while (!evicted_l_set->empty()) {
                auto l = *evicted_l_set->begin();
                evicted_l_set->erase(l);

                left_map->erase(l);

                left_mutex_map->at(l)->lock();
                for (uint32_t index = k; index <= i; ++index) {
                    left_index_map->at(l)->insert(index, k);
                }
                left_mutex_map->at(l)->unlock();

                for (const auto &[r, e]: *B->get_left_vertex(l)->get_edge_map()) {
                    if (right_map->count(r) && right_map->at(r) >= k) {
                        --right_map->at(r);
                        if (right_map->at(r) < k) {
                            right_map->erase(r);

                            right_mutex_map->at(r)->lock();
                            right_index_map->at(r)->insert(k, i);
                            right_mutex_map->at(r)->unlock();

                            for (const auto &[l2, e2]: *B->get_right_vertex(r)->get_edge_map()) {
                                if (left_map->count(l2) && left_map->at(l2) > i) {
                                    --left_map->at(l2);
                                    if (left_map->at(l2) == i) {
                                        evicted_l_set->emplace(l2);
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
        return max_i;
    }

    uint32_t branch_bipartite_core_decomposition::find_right_core(const shared_ptr<abstract_bipartite_graph> &B,
                                                                  const shared_ptr<unordered_map<uint32_t, uint32_t>> &left_map,
                                                                  const shared_ptr<unordered_map<uint32_t, uint32_t>> &right_map,
                                                                  const shared_ptr<unordered_map<uint32_t, shared_ptr<bipartite_core_left_store_index>>> &left_index_map,
                                                                  const shared_ptr<unordered_map<uint32_t, shared_ptr<bipartite_core_right_store_index>>> &right_index_map,
                                                                  uint32_t k) {
        auto evicted_r_set = make_shared<unordered_set<uint32_t>>();

        uint32_t max_j = k;
        while (!right_map->empty()) {
            uint32_t  j = right_map->begin()->second;

            for (const auto&[r, r_degree]: *right_map) {
                if (r_degree < j) {
                    j = r_degree;

                    evicted_r_set->clear();
                    evicted_r_set->insert(r);
                } else if(r_degree == j)
                {
                    evicted_r_set->insert(r);
                }
            }
            max_j = j;

            while (!evicted_r_set->empty()) {
                auto r = *evicted_r_set->begin();
                evicted_r_set->erase(r);

                right_map->erase(r);

                for (uint32_t index = k; index <= j; ++index) {
                    right_index_map->at(r)->insert(index, k);
                }

                for (const auto &[l, e]: *B->get_right_vertex(r)->get_edge_map()) {
                    if (left_map->count(l) && left_map->at(l) >= k) {
                        --left_map->at(l);
                        if (left_map->at(l) < k) {
                            left_map->erase(l);
                            left_index_map->at(l)->insert(k, j);

                            for (const auto &[r2, e2]: *B->get_left_vertex(l)->get_edge_map()) {
                                if (right_map->count(r2) && right_map->at(r2) > j) {
                                    --right_map->at(r2);
                                    if (right_map->at(r2) == j) {
                                        evicted_r_set->insert(r2);
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
        return max_j;
    }


    uint32_t branch_bipartite_core_decomposition::find_right_core(const shared_ptr<abstract_bipartite_graph> &B,
                                                                  const shared_ptr<unordered_map<uint32_t, shared_ptr<mutex>>>& left_mutex_map,
                                                                  const shared_ptr<unordered_map<uint32_t, shared_ptr<mutex>>>& right_mutex_map,
                                                                 const shared_ptr<unordered_map<uint32_t, uint32_t>> &left_map,
                                                                 const shared_ptr<unordered_map<uint32_t, uint32_t>> &right_map,
                                                                 const shared_ptr<unordered_map<uint32_t, shared_ptr<bipartite_core_left_store_index>>>& left_index_map,
                                                                 const shared_ptr<unordered_map<uint32_t, shared_ptr<bipartite_core_right_store_index>>>& right_index_map,
                                                                 uint32_t k) {
        auto evicted_r_set = make_shared<unordered_set<uint32_t>>();

        uint32_t max_j = k;
        while (!right_map->empty()) {
            uint32_t  j = right_map->begin()->second;

            for (const auto&[r, r_degree]: *right_map) {
                if (r_degree < j) {
                    j = r_degree;

                    evicted_r_set->clear();
                    evicted_r_set->insert(r);
                } else if(r_degree == j)
                {
                    evicted_r_set->insert(r);
                }
            }
            max_j = j;

            while (!evicted_r_set->empty()) {
                auto r = *evicted_r_set->begin();
                evicted_r_set->erase(r);

                right_map->erase(r);

                right_mutex_map->at(r)->lock();
                for (uint32_t index = k; index <= j; ++index) {
                    right_index_map->at(r)->insert(index, k);
                }
                right_mutex_map->at(r)->unlock();

                for (const auto &[l, e]: *B->get_right_vertex(r)->get_edge_map()) {
                    if (left_map->count(l) && left_map->at(l) >= k) {
                        --left_map->at(l);
                        if (left_map->at(l) < k) {
                            left_mutex_map->at(l)->lock();
                            left_index_map->at(l)->insert(k, j);
                            left_mutex_map->at(l)->unlock();

                            for (const auto &[r2, e2]: *B->get_left_vertex(l)->get_edge_map()) {
                                if (right_map->count(r2) && right_map->at(r2) > j) {
                                    --right_map->at(r2);
                                    if (right_map->at(r2) == j) {
                                        evicted_r_set->emplace(r2);
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
        return max_j;
    }

    uint32_t branch_bipartite_core_decomposition::find_right_core(const shared_ptr<abstract_bipartite_graph> &B,
                                                                  const shared_ptr<unordered_map<uint32_t, shared_ptr<mutex>>> &left_mutex_map,
                                                                  const shared_ptr<unordered_map<uint32_t, shared_ptr<mutex>>> &right_mutex_map,
                                                                  const shared_ptr<unordered_set<uint32_t>> &l_set,
                                                                  const shared_ptr<unordered_set<uint32_t>> &r_set,
                                                                  const shared_ptr<unordered_map<uint32_t, uint32_t>> &vertex_degree_map,
                                                                  const shared_ptr<unordered_map<uint32_t, shared_ptr<bipartite_core_left_store_index>>> &left_index_map,
                                                                  const shared_ptr<unordered_map<uint32_t, shared_ptr<bipartite_core_right_store_index>>> &right_index_map,
                                                                  const shared_ptr<unordered_map<uint32_t, shared_ptr<extend_list<int, uint32_t>>>> &k_oder_map,
                                                                  const shared_ptr<unordered_map<uint32_t, uint32_t>> &k_rem_degree_map,
                                                                  const shared_ptr<unordered_map<uint32_t, uint32_t>> &k_core_degree_map,
                                                                  uint32_t k) {
        auto evicted_l_set = make_shared<unordered_set<uint32_t>>();
        auto evicted_r_set = make_shared<unordered_set<uint32_t>>();

        uint32_t max_j = k;
        while (!r_set->empty()) {
            uint32_t j = UINT32_MAX;

            for (const auto &l: *l_set) {
                k_core_degree_map->at(l) = vertex_degree_map->at(l);
            }

            for (const auto &r: *r_set) {
                k_core_degree_map->at(r) = vertex_degree_map->at(r);

                if (vertex_degree_map->at(r) < j) {
                    j = vertex_degree_map->at(r);

                    evicted_r_set->clear();
                    evicted_r_set->insert(r);
                } else if (vertex_degree_map->at(r) == j) {
                    evicted_r_set->insert(r);
                }
            }

            if (!k_oder_map->count(j)) {
                k_oder_map->insert({j, make_shared<extend_list<int, uint32_t>>()});
            }

            max_j = j;

            while (!evicted_l_set->empty() || !evicted_r_set->empty()) {

                while (!evicted_r_set->empty()) {
                    auto r = *evicted_r_set->begin();
                    evicted_r_set->erase(r);

                    right_mutex_map->at(r)->lock();
                    for (uint32_t index = k; index <= j; ++index) {
                        right_index_map->at(r)->insert(index, k);
                    }
                    right_mutex_map->at(r)->unlock();;

                    k_oder_map->at(j)->push_back(r);
                    k_rem_degree_map->at(r) = vertex_degree_map->at(r);
                    r_set->erase(r);
                    vertex_degree_map->erase(r);

                    for (const auto &[l, e]: *B->get_right_vertex(r)->get_edge_map()) {
                        if (vertex_degree_map->count(l)) {
                            --vertex_degree_map->at(l);
                            if (vertex_degree_map->at(l) < k) {
                                evicted_l_set->insert(l);
                            }
                        }
                    }
                }



                while (!evicted_l_set->empty()) {
                    auto l = *evicted_l_set->begin();
                    evicted_l_set->erase(l);

                    left_mutex_map->at(l)->lock();
                    left_index_map->at(l)->insert(k, j);
                    left_mutex_map->at(l)->unlock();

                    k_oder_map->at(j)->push_back(l);
                    k_rem_degree_map->at(l) = vertex_degree_map->at(l);
                    l_set->erase(l);
                    vertex_degree_map->erase(l);

                    for (const auto &[r, e]: *B->get_left_vertex(l)->get_edge_map()) {
                        if (r_set->count(r)) {
                            --vertex_degree_map->at(r);
                            if (vertex_degree_map->at(r) == j) {
                                evicted_r_set->insert(r);
                            }
                        }
                    }
                }
            }
        }
        return max_j;
    }

    uint32_t branch_bipartite_core_decomposition::find_right_core(const shared_ptr<abstract_bipartite_graph> &B,
                                                                  const shared_ptr<unordered_map<uint32_t, shared_ptr<mutex>>> &left_mutex_map,
                                                                  const shared_ptr<unordered_map<uint32_t, shared_ptr<mutex>>> &right_mutex_map,
                                                                  const shared_ptr<unordered_map<uint32_t, uint32_t>> &left_map,
                                                                  const shared_ptr<unordered_map<uint32_t, uint32_t>> &right_map,
                                                                  const shared_ptr<unordered_map<uint32_t, shared_ptr<bipartite_core_branch_store_index>>> &left_index_map,
                                                                  const shared_ptr<unordered_map<uint32_t, shared_ptr<bipartite_core_branch_store_index>>> &right_index_map,
                                                                  uint32_t k) {
        auto evicted_r_set = make_shared<unordered_set<uint32_t>>();

        uint32_t max_j = k;
        while (!right_map->empty()) {
            uint32_t j = right_map->begin()->second;

            for (const auto &[r, r_degree]: *right_map) {
                if (r_degree < j) {
                    j = r_degree;

                    evicted_r_set->clear();
                    evicted_r_set->insert(r);
                } else if (r_degree == j) {
                    evicted_r_set->insert(r);
                }
            }
            max_j = j;

            while (!evicted_r_set->empty()) {
                auto r = *evicted_r_set->begin();
                evicted_r_set->erase(r);

                right_map->erase(r);

                right_mutex_map->at(r)->lock();
                for (uint32_t index = k; index <= j; ++index) {
                    right_index_map->at(r)->insert(index, k);
                }
                right_mutex_map->at(r)->unlock();

                for (const auto &[l, e]: *B->get_right_vertex(r)->get_edge_map()) {
                    if (left_map->count(l) && left_map->at(l) >= k) {
                        --left_map->at(l);
                        if (left_map->at(l) < k) {
                            left_mutex_map->at(l)->lock();
                            left_index_map->at(l)->insert(k, j);
                            left_mutex_map->at(l)->unlock();

                            for (const auto &[r2, e2]: *B->get_left_vertex(l)->get_edge_map()) {
                                if (right_map->count(r2) && right_map->at(r2) > j) {
                                    --right_map->at(r2);
                                    if (right_map->at(r2) == j) {
                                        evicted_r_set->emplace(r2);
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
        return max_j;
    }
}