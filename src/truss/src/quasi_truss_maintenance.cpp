
#include "truss/quasi_truss_maintenance.h"

namespace scnu {

    void quasi_truss_maintenance::candidate_graph_finding(const shared_ptr<abstract_graph> &G,
                                                          const shared_ptr<unordered_map<shared_ptr<abstract_edge>, shared_ptr<mutex>>> &edge_mutex_map,
                                                          const shared_ptr<unordered_map<shared_ptr<abstract_edge>, uint32_t>> &edge_rank_map,
                                                          const shared_ptr<unordered_map<shared_ptr<abstract_edge>, uint32_t>> &edge_truss_map,
                                                          const shared_ptr<unordered_set<shared_ptr<abstract_edge>>> &edge_set,
                                                          const shared_ptr<unordered_set<shared_ptr<abstract_edge>>> &candidate_edge_set,
                                                          const shared_ptr<unordered_map<shared_ptr<abstract_edge>, uint32_t>> &candidate_edge_support_map,
                                                          const shared_ptr<unordered_map<shared_ptr<abstract_edge>, uint32_t>> &TS,
                                                          uint32_t k,
                                                          const shared_ptr<thread_pool> &pool) {
        auto clear_set = make_shared<unordered_set<shared_ptr<abstract_edge>>>();

        auto thread_number = pool->get_thread_number();
        auto global_mutex = make_shared<mutex>();
        /**
         * @brief find the k-insert graph
         */
        auto current_edge_set = make_shared<unordered_set<shared_ptr<abstract_edge>>>();
        {
            auto removed_edge_set = make_shared<unordered_set<shared_ptr<abstract_edge>>>();
            auto location_vector = pool->split_task(edge_set);
            for (uint32_t i = 0; i < thread_number; ++i) {
                pool->submit_task([=] {
                    auto sub_current_edge_set = make_shared<unordered_set<shared_ptr<abstract_edge>>>();
                    auto sub_removed_edge_set = make_shared<unordered_set<shared_ptr<abstract_edge>>>();

                    auto &sub_begin = *location_vector->at(i);
                    auto &sub_end = *location_vector->at(i + 1);

                    for (auto iter = sub_begin; iter != sub_end; ++iter) {
                        auto &e = *iter;

                        if (edge_truss_map->at(e) < k - 1) {
                            sub_removed_edge_set->insert(e);
                        } else {
                            if (edge_truss_map->at(e) == k - 1) {
                                clear_set->insert(e);
                                TS->at(e) = edge_support_computation(G, e, edge_truss_map, k - 1);
                                if (TS->at(e) >= k - 2) {
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
            for (const auto &e: *removed_edge_set) {
                edge_set->erase(e);
            }
        }

        /**
         * @brief find a set of candidate edges
         */
        auto visited_edge_set = make_shared<unordered_set<shared_ptr<abstract_edge>>>();
        auto evicted_edge_set = make_shared<unordered_set<shared_ptr<abstract_edge>>>();
        while (!current_edge_set->empty()) {
            auto next_edge_set = make_shared<unordered_set<shared_ptr<abstract_edge>>>();
            auto layer_candidate_edge_set = make_shared<unordered_set<shared_ptr<abstract_edge>>>();
            auto invalid_edge_set = make_shared<unordered_set<shared_ptr<abstract_edge>>>();
            auto location_vector = pool->split_task(current_edge_set);
            for (uint32_t i = 0; i < thread_number; ++i) {
                pool->submit_task([=] {
                    auto sub_candidate_edge_set = make_shared<unordered_set<shared_ptr<abstract_edge>>>();
                    auto sub_next_edge_set = make_shared<unordered_set<shared_ptr<abstract_edge>>>();
                    auto sub_invalid_edge_set = make_shared<unordered_set<shared_ptr<abstract_edge>>>();

                    auto &sub_begin = *location_vector->at(i);
                    auto &sub_end = *location_vector->at(i + 1);

                    for (auto iter = sub_begin; iter != sub_end; ++iter) {
                        auto &e1 = *iter;

                        auto u = e1->get_source_vertex_id();
                        auto v = e1->get_destination_vertex_id();

                        if (G->get_vertex(u)->get_degree() > G->get_vertex(v)->get_degree()) {
                            swap(u, v);
                        }

                        uint32_t e1_count = 0;
                        auto e1_set = make_shared<unordered_set<shared_ptr<abstract_edge>>>();

                        for (const auto &[w, e2]: *G->get_vertex(u)->get_edge_map()) {
                            if (w == v || evicted_edge_set->count(e2) || edge_truss_map->at(e2) < k - 1) {
                                continue;
                            }

                            if (edge_truss_map->at(e2) == k - 1) {
                                if (TS->at(e2) == 0) {
                                    TS->at(e2) = edge_support_computation(G, e2, edge_truss_map, k - 1);
                                }
                                if (TS->at(e2) < k - 2) {
                                    continue;
                                }
                            }

                            auto e3 = G->get_edge(v, w);
                            if (!e3 || evicted_edge_set->count(e3) || edge_truss_map->at(e3) < k - 1) {
                                continue;
                            }

                            if (edge_truss_map->at(e3) == k - 1) {
                                if (TS->at(e3) == 0) {
                                    TS->at(e3) = edge_support_computation(G, e3, edge_truss_map, k - 1);
                                }
                                if (TS->at(e3) < k - 2) {
                                    continue;
                                }
                            }

                            ++e1_count;

                            if (edge_truss_map->at(e2) == k - 1 && !candidate_edge_set->count(e2) &&
                                !current_edge_set->count(e2) && !visited_edge_set->count(e2)) {
                                e1_set->insert(e2);
                            }


                            if (edge_truss_map->at(e3) == k - 1 && !candidate_edge_set->count(e3) &&
                                !current_edge_set->count(e3) && !visited_edge_set->count(e3)) {
                                e1_set->insert(e3);
                            }
                        }
                        candidate_edge_support_map->at(e1) = e1_count;
                        if (e1_count >= k - 2) {
                            sub_candidate_edge_set->insert(e1);
                            sub_next_edge_set->merge(*e1_set);
                        } else {
                            sub_invalid_edge_set->insert(e1);
                        }
                    }

                    global_mutex->lock();
                    layer_candidate_edge_set->merge(*sub_candidate_edge_set);
                    next_edge_set->merge(*sub_next_edge_set);
                    invalid_edge_set->merge(*sub_invalid_edge_set);
                    global_mutex->unlock();
                });
            }
            pool->barrier();
            candidate_edge_set->merge(*layer_candidate_edge_set);
            remove_unsatisfied_edges(G, edge_mutex_map, edge_rank_map, edge_truss_map, invalid_edge_set,
                                     candidate_edge_set, candidate_edge_support_map, evicted_edge_set, k, pool);
            visited_edge_set->merge(*current_edge_set);
            swap(*current_edge_set, *next_edge_set);
        }

        for (const auto &e: *clear_set) {
            TS->at(e) = 0;
        }
    }


    void quasi_truss_maintenance::candidate_graph_finding2(const shared_ptr<abstract_graph> &G,
                                                           const shared_ptr<unordered_map<shared_ptr<abstract_edge>, shared_ptr<mutex>>> &edge_mutex_map,
                                                           const shared_ptr<unordered_map<shared_ptr<abstract_edge>, uint32_t>> &edge_rank_map,
                                                           const shared_ptr<unordered_map<shared_ptr<abstract_edge>, uint32_t>> &edge_truss_map,
                                                           const shared_ptr<unordered_set<shared_ptr<abstract_edge>>> &edge_set,
                                                           const shared_ptr<unordered_set<shared_ptr<abstract_edge>>> &candidate_edge_set,
                                                           const shared_ptr<unordered_map<shared_ptr<abstract_edge>, uint32_t>> &candidate_edge_support_map,
                                                           const shared_ptr<unordered_map<shared_ptr<abstract_edge>, uint32_t>> &TS,
                                                           uint32_t k,
                                                           const shared_ptr<thread_pool> &pool) {
        auto thread_number = pool->get_thread_number();
        auto global_mutex = make_shared<mutex>();
        /**
         * @brief find the k-insert graph
         */
        auto current_edge_set = make_shared<unordered_set<shared_ptr<abstract_edge>>>();
        {
            auto removed_edge_set = make_shared<unordered_set<shared_ptr<abstract_edge>>>();
            auto location_vector = pool->split_task(edge_set);
            for (uint32_t i = 0; i < thread_number; ++i) {
                pool->submit_task([=] {
                    auto sub_current_edge_set = make_shared<unordered_set<shared_ptr<abstract_edge>>>();
                    auto sub_removed_edge_set = make_shared<unordered_set<shared_ptr<abstract_edge>>>();

                    auto &sub_begin = *location_vector->at(i);
                    auto &sub_end = *location_vector->at(i + 1);

                    for (auto iter = sub_begin; iter != sub_end; ++iter) {
                        auto &e = *iter;

                        if (edge_truss_map->at(e) < k - 1) {
                            sub_removed_edge_set->insert(e);
                        } else {
                            if (edge_truss_map->at(e) == k - 1 && TS->at(e) >= k - 2) {
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
            for (const auto &e: *removed_edge_set) {
                edge_set->erase(e);
            }
        }

        /**
         * @brief find a set of candidate edges
         */
        auto visited_edge_set = make_shared<unordered_set<shared_ptr<abstract_edge>>>();
        auto evicted_edge_set = make_shared<unordered_set<shared_ptr<abstract_edge>>>();
        while (!current_edge_set->empty()) {
            auto next_edge_set = make_shared<unordered_set<shared_ptr<abstract_edge>>>();
            auto layer_candidate_edge_set = make_shared<unordered_set<shared_ptr<abstract_edge>>>();
            auto invalid_edge_set = make_shared<unordered_set<shared_ptr<abstract_edge>>>();
            auto location_vector = pool->split_task(current_edge_set);
            for (uint32_t i = 0; i < thread_number; ++i) {
                pool->submit_task([=] {
                    auto sub_candidate_edge_set = make_shared<unordered_set<shared_ptr<abstract_edge>>>();
                    auto sub_next_edge_set = make_shared<unordered_set<shared_ptr<abstract_edge>>>();
                    auto sub_invalid_edge_set = make_shared<unordered_set<shared_ptr<abstract_edge>>>();

                    auto &sub_begin = *location_vector->at(i);
                    auto &sub_end = *location_vector->at(i + 1);

                    for (auto iter = sub_begin; iter != sub_end; ++iter) {
                        auto &e1 = *iter;

                        auto u = e1->get_source_vertex_id();
                        auto v = e1->get_destination_vertex_id();

                        uint32_t e1_count = 0;
                        auto e1_set = make_shared<unordered_set<shared_ptr<abstract_edge>>>();

                        for (const auto &[w, e2]: *G->get_vertex(u)->get_edge_map()) {
                            if (w == v || evicted_edge_set->count(e2) || edge_truss_map->at(e2) < k - 1 ||
                                (edge_truss_map->at(e2) == k - 1 && TS->at(e2) < k - 2)) {
                                continue;
                            }

                            auto e3 = G->get_edge(v, w);
                            if (!e3 || evicted_edge_set->count(e3) || edge_truss_map->at(e3) < k - 1 ||
                                (edge_truss_map->at(e3) == k - 1 && TS->at(e3) < k - 2)) {
                                continue;
                            }

                            ++e1_count;

                            if (edge_truss_map->at(e2) == k - 1 && !candidate_edge_set->count(e2) &&
                                !current_edge_set->count(e2) && !visited_edge_set->count(e2)) {
                                e1_set->insert(e2);
                            }


                            if (edge_truss_map->at(e3) == k - 1 && !candidate_edge_set->count(e3) &&
                                !current_edge_set->count(e3) && !visited_edge_set->count(e3)) {
                                e1_set->insert(e3);
                            }
                        }
                        if (e1_count >= k - 2) {
                            sub_candidate_edge_set->insert(e1);
                            candidate_edge_support_map->at(e1) = e1_count;
                            sub_next_edge_set->merge(*e1_set);
                        } else {
                            sub_invalid_edge_set->insert(e1);
                        }
                    }

                    global_mutex->lock();
                    layer_candidate_edge_set->merge(*sub_candidate_edge_set);
                    next_edge_set->merge(*sub_next_edge_set);
                    invalid_edge_set->merge(*sub_invalid_edge_set);
                    global_mutex->unlock();
                });
            }
            pool->barrier();
            candidate_edge_set->merge(*layer_candidate_edge_set);
            remove_unsatisfied_edges(G, edge_mutex_map, edge_rank_map, edge_truss_map, invalid_edge_set,
                                     candidate_edge_set, candidate_edge_support_map, evicted_edge_set, k, pool);
            visited_edge_set->merge(*current_edge_set);
            swap(*current_edge_set, *next_edge_set);
        }
    }

    uint32_t quasi_truss_maintenance::edge_support_computation(const shared_ptr<abstract_graph> &G,
                                                               const shared_ptr<abstract_edge> &e,
                                                               const shared_ptr<unordered_map<shared_ptr<abstract_edge>, uint32_t>> &edge_truss_map,
                                                               uint32_t k) {

        const auto &e1 = e;
        uint32_t e1_support = 0;

        auto u = e1->get_source_vertex_id();
        auto v = e1->get_destination_vertex_id();

        for (const auto &[w, e2]: *G->get_vertex(u)->get_edge_map()) {
            if (w == v || edge_truss_map->at(e2) < k) {
                continue;
            }

            auto e3 = G->get_edge(v, w);
            if (!e3 || edge_truss_map->at(e3) < k) {
                continue;
            }

            ++e1_support;
        }
        return e1_support;
    }


    /**
     * @details compute the support of a given edge in k-wing
     * @param G
     * @param e
     * @param edge_wing_map
     * @param TS
     * @param k
     * @return
     */
    uint32_t quasi_truss_maintenance::edge_support_computation(const shared_ptr<abstract_graph> &G,
                                                               const shared_ptr<abstract_edge> &e,
                                                               const shared_ptr<unordered_map<shared_ptr<abstract_edge>, uint32_t>> &edge_wing_map,
                                                               const shared_ptr<unordered_map<shared_ptr<abstract_edge>, uint32_t>> &TS,
                                                               uint32_t k) {
        if (TS->at(e) == 0) {
            const auto &e1 = e;
            uint32_t e1_support = 0;

            auto u = e1->get_source_vertex_id();
            auto v = e1->get_destination_vertex_id();

            for (const auto &[w, e2]: *G->get_vertex(u)->get_edge_map()) {
                if (w == v || edge_wing_map->at(e2) < k) {
                    continue;
                }
                auto e3 = G->get_edge(v, w);
                if (!e3 || edge_wing_map->at(e3) < k) {
                    continue;
                }

                ++e1_support;
            }
            TS->at(e1) = e1_support;
        }
        return TS->at(e);
    }

    void
    quasi_truss_maintenance::init(const shared_ptr<unordered_map<shared_ptr<abstract_edge>, uint32_t>> &edge_truss_map,
                                  const shared_ptr<unordered_map<shared_ptr<abstract_edge>, shared_ptr<mutex>>> &edge_mutex_map,
                                  const shared_ptr<unordered_map<shared_ptr<abstract_edge>, uint32_t>> &edge_rank_map,
                                  const shared_ptr<unordered_map<shared_ptr<abstract_edge>, uint32_t>> &edge_support_map,
                                  const shared_ptr<thread_pool> &pool) {
        pool->submit_task([=] {
            for (const auto &[e, truss_number]: *edge_truss_map) {
                edge_mutex_map->insert({e, make_shared<mutex>()});
            }
        });

        pool->submit_task([=] {
            for (const auto &[e, truss_number]: *edge_truss_map) {
                edge_rank_map->insert({e, UINT32_MAX});
            }
        });

        pool->submit_task([=] {
            for (const auto &[e, truss_number]: *edge_truss_map) {
                edge_rank_map->insert({e, 0});
            }
        });

        pool->barrier();
    }

    void
    quasi_truss_maintenance::init(const shared_ptr<unordered_map<shared_ptr<abstract_edge>, uint32_t>> &edge_truss_map,
                                  const shared_ptr<unordered_map<shared_ptr<abstract_edge>, shared_ptr<mutex>>> &edge_mutex_map,
                                  const shared_ptr<unordered_map<shared_ptr<abstract_edge>, uint32_t>> &edge_rank_map,
                                  const shared_ptr<unordered_map<shared_ptr<abstract_edge>, uint32_t>> &edge_support_map,
                                  const shared_ptr<unordered_map<shared_ptr<abstract_edge>, uint32_t>> &TS,
                                  const shared_ptr<thread_pool> &pool) {
        pool->submit_task([=] {
            for (const auto &[e, truss_number]: *edge_truss_map) {
                edge_mutex_map->insert({e, make_shared<mutex>()});
            }
        });

        pool->submit_task([=] {
            for (const auto &[e, truss_number]: *edge_truss_map) {
                edge_rank_map->insert({e, UINT32_MAX});
            }
        });

        pool->submit_task([=] {
            for (const auto &[e, truss_number]: *edge_truss_map) {
                edge_rank_map->insert({e, 0});
            }
        });

        pool->submit_task([=] {
            for (const auto &[e, truss_number]: *edge_truss_map) {
                TS->insert({e, 0});
            }
        });
        pool->barrier();
    }

    /**
     * @details a basic insert method
     * @param G
     * @param edge_set
     * @param edge_truss_map
     * @param previous_k_max
     * @param thread_number
     */
    void quasi_truss_maintenance::insert(const shared_ptr<abstract_graph> &G,
                                         const shared_ptr<unordered_map<shared_ptr<abstract_edge>, shared_ptr<mutex>>> &edge_mutex_map,
                                         const shared_ptr<unordered_set<shared_ptr<abstract_edge>>> &edge_set,
                                         const shared_ptr<unordered_map<shared_ptr<abstract_edge>, uint32_t>> &edge_rank_map,
                                         const shared_ptr<unordered_map<shared_ptr<abstract_edge>, uint32_t>> &candidate_edge_support_map,
                                         const shared_ptr<unordered_map<shared_ptr<abstract_edge>, uint32_t>> &edge_truss_map,
                                         const shared_ptr<unordered_map<shared_ptr<abstract_edge>, uint32_t>> &TS,
                                         const shared_ptr<uint32_t> &previous_k_max,
                                         const shared_ptr<thread_pool> &pool) {
        auto thread_number = pool->get_thread_number();
        auto global_mutex = make_shared<mutex>();

        G->insert_edge_collection(edge_set);
        auto rank_id = make_shared<uint32_t>(0);
        auto affected_edge_set = make_shared<unordered_set<shared_ptr<abstract_edge>>>();
        {
            set_edge_rank(edge_set, edge_rank_map, rank_id);
            auto location_vector = pool->split_task(edge_set);
            for (uint32_t i = 0; i < thread_number; ++i) {
                pool->submit_task([=] {
                    auto sub_affected_set = make_shared<unordered_set<shared_ptr<abstract_edge>>>();

                    auto &sub_begin = *location_vector->at(i);
                    auto &sub_end = *location_vector->at(i + 1);

                    for (auto iter = sub_begin; iter != sub_end; ++iter) {
                        const auto &e1 = *iter;
                        auto u = e1->get_source_vertex_id();
                        auto v = e1->get_destination_vertex_id();

                        for (const auto &[w, e2]: *G->get_vertex(u)->get_edge_map()) {
                            if (w == v || edge_rank_map->at(e2) < edge_rank_map->at(e1)) {
                                continue;
                            }
                            auto e3 = G->get_edge(v, w);
                            if (!e3 || edge_rank_map->at(e3) < edge_rank_map->at(e1)) {
                                continue;
                            }

                            sub_affected_set->insert(e1);
                            sub_affected_set->insert(e2);
                            sub_affected_set->insert(e3);
                        }
                    }

                    global_mutex->lock();
                    affected_edge_set->merge(*sub_affected_set);
                    global_mutex->unlock();
                });
            }
            pool->barrier();
            for (const auto &e: *edge_set) {
                edge_rank_map->at(e) = UINT32_MAX;
                if (affected_edge_set->count(e)) {
                    edge_truss_map->insert({e, 3});
                } else {
                    edge_truss_map->insert({e, 2});
                }
            }
        }


        /**
         * @brief update the remainder k-wings
         */
        uint32_t k = 2;
        while (!affected_edge_set->empty()) {
            if (k <= *previous_k_max) {
                auto candidate_edge_set = make_shared<unordered_set<shared_ptr<abstract_edge>>>();
                candidate_graph_finding(G, edge_mutex_map, edge_rank_map, edge_truss_map, affected_edge_set,
                                        candidate_edge_set, candidate_edge_support_map, TS, k, pool);
                for (const auto &e: *candidate_edge_set) {
                    edge_truss_map->at(e) = k;
                    TS->at(e) = candidate_edge_support_map->at(e);
                    candidate_edge_support_map->at(e) = 0;
                }
                /**
                  * @brief continue the loop
                  */
                ++k;
            } else {
                auto candidate_edge_set = make_shared<unordered_set<shared_ptr<abstract_edge>>>();
                candidate_graph_finding(G, edge_mutex_map, edge_rank_map, edge_truss_map, affected_edge_set,
                                        candidate_edge_set, candidate_edge_support_map, TS, k, pool);
                *previous_k_max = partial_truss_decomposition(G, edge_mutex_map, edge_rank_map, edge_truss_map,
                                                              candidate_edge_set, candidate_edge_support_map, k, pool);
                break;
            }
        }
    }

    /**
    * @details a basic insert method
    * @param G
    * @param edge_set
    * @param edge_truss_map
    * @param previous_k_max
    * @param thread_number
    */
    void quasi_truss_maintenance::insert2(const shared_ptr<abstract_graph> &G,
                                          const shared_ptr<unordered_map<shared_ptr<abstract_edge>, shared_ptr<mutex>>> &edge_mutex_map,
                                          const shared_ptr<unordered_set<shared_ptr<abstract_edge>>> &edge_set,
                                          const shared_ptr<unordered_map<shared_ptr<abstract_edge>, uint32_t>> &edge_rank_map,
                                          const shared_ptr<unordered_map<shared_ptr<abstract_edge>, uint32_t>> &candidate_edge_support_map,
                                          const shared_ptr<unordered_map<shared_ptr<abstract_edge>, uint32_t>> &edge_truss_map,
                                          const shared_ptr<unordered_map<shared_ptr<abstract_edge>, uint32_t>> &TS,
                                          const shared_ptr<uint32_t> &previous_k_max,
                                          const shared_ptr<thread_pool> &pool) {
        auto thread_number = pool->get_thread_number();
        auto global_mutex = make_shared<mutex>();

        G->insert_edge_collection(edge_set);
        auto affected_edge_set = make_shared<unordered_set<shared_ptr<abstract_edge>>>();
        {
            set_edge_rank(edge_set, edge_rank_map, 0);
            auto location_vector = pool->split_task(edge_set);
            for (uint32_t i = 0; i < thread_number; ++i) {
                pool->submit_task([=] {
                    auto sub_affected_set = make_shared<unordered_set<shared_ptr<abstract_edge>>>();

                    auto &sub_begin = *location_vector->at(i);
                    auto &sub_end = *location_vector->at(i + 1);

                    for (auto iter = sub_begin; iter != sub_end; ++iter) {
                        const auto &e1 = *iter;
                        auto u = e1->get_source_vertex_id();
                        auto v = e1->get_destination_vertex_id();

                        for (const auto &[w, e2]: *G->get_vertex(u)->get_edge_map()) {
                            if (w == v || edge_rank_map->at(e2) < edge_rank_map->at(e1)) {
                                continue;
                            }

                            auto e3 = G->get_edge(v, w);
                            if (!e3 || edge_rank_map->at(e3) < edge_rank_map->at(e1)) {
                                continue;
                            }

                            edge_mutex_map->at(e1)->lock();
                            ++TS->at(e1);
                            edge_mutex_map->at(e1)->unlock();

                            if (edge_truss_map->at(e2) <= 1) {
                                edge_mutex_map->at(e2)->lock();
                                ++TS->at(e2);
                                edge_mutex_map->at(e2)->unlock();
                            }

                            if (edge_truss_map->at(e3) <= 1) {
                                edge_mutex_map->at(e3)->lock();
                                ++TS->at(e3);
                                edge_mutex_map->at(e3)->unlock();
                            }


                            sub_affected_set->insert(e1);
                            sub_affected_set->insert(e2);
                            sub_affected_set->insert(e3);
                        }
                    }

                    global_mutex->lock();
                    affected_edge_set->merge(*sub_affected_set);
                    global_mutex->unlock();
                });
            }
            pool->barrier();
            for (const auto &e: *edge_set) {
                edge_rank_map->at(e) = UINT32_MAX;
                if (affected_edge_set->count(e)) {
                    edge_truss_map->insert({e, 3});
                } else {
                    edge_truss_map->insert({e, 2});
                }
            }
        }

        /**
         * @brief update the remainder k-wings
         */
        uint32_t k = 2;
        while (!affected_edge_set->empty()) {
            if (k <= *previous_k_max) {
                auto candidate_edge_set = make_shared<unordered_set<shared_ptr<abstract_edge>>>();
                candidate_graph_finding2(G, edge_mutex_map, edge_rank_map, edge_truss_map, affected_edge_set,
                                         candidate_edge_set, candidate_edge_support_map, TS, k, pool);
                update_edge_truss_support(G, edge_mutex_map, edge_rank_map, edge_truss_map, TS, candidate_edge_set, k,
                                          pool);
                for (const auto &e: *candidate_edge_set) {
                    edge_truss_map->at(e) = k;
                    TS->at(e) = candidate_edge_support_map->at(e);
                    candidate_edge_support_map->at(e) = 0;
                }
                /**
                  * @brief continue the loop
                  */
                ++k;
            } else {
                auto candidate_edge_set = make_shared<unordered_set<shared_ptr<abstract_edge>>>();
                candidate_graph_finding2(G, edge_mutex_map, edge_rank_map, edge_truss_map, affected_edge_set,
                                         candidate_edge_set, candidate_edge_support_map, TS, k, pool);
                *previous_k_max = partial_truss_decomposition2(G, edge_mutex_map, edge_rank_map, edge_truss_map,
                                                               candidate_edge_set, candidate_edge_support_map, TS, k,
                                                               pool);
                break;
            }
        }
    }


    uint32_t quasi_truss_maintenance::partial_truss_decomposition(const shared_ptr<abstract_graph> &G,
                                                                  const shared_ptr<unordered_map<shared_ptr<abstract_edge>, shared_ptr<mutex>>> &edge_mutex_map,
                                                                  const shared_ptr<unordered_map<shared_ptr<abstract_edge>, uint32_t>> &edge_rank_map,
                                                                  const shared_ptr<unordered_map<shared_ptr<abstract_edge>, uint32_t>> &edge_truss_map,
                                                                  const shared_ptr<unordered_set<shared_ptr<abstract_edge>>> &candidate_edge_set,
                                                                  const shared_ptr<unordered_map<shared_ptr<abstract_edge>, uint32_t>> &candidate_edge_support_map,
                                                                  uint32_t k,
                                                                  const shared_ptr<thread_pool> &pool) {
        uint32_t k_max = 0;
        while (!candidate_edge_set->empty()) {
            k_max = k;

            auto current_edge_set = make_shared<unordered_set<shared_ptr<abstract_edge>>>();
            for (const auto &e: *candidate_edge_set) {
                edge_truss_map->at(e) = k;
                if (candidate_edge_support_map->at(e) < k - 1) {
                    current_edge_set->insert(e);
                }
            }
            auto evicted_edge_set = make_shared<unordered_set<shared_ptr<abstract_edge>>>();
            remove_unsatisfied_edges(G, edge_mutex_map, edge_rank_map, edge_truss_map, current_edge_set,
                                     candidate_edge_set, candidate_edge_support_map, evicted_edge_set, k + 1, pool);
            ++k;
        }
        return k_max;
    }

    uint32_t quasi_truss_maintenance::partial_truss_decomposition2(const shared_ptr<abstract_graph> &G,
                                                                   const shared_ptr<unordered_map<shared_ptr<abstract_edge>, shared_ptr<mutex>>> &edge_mutex_map,
                                                                   const shared_ptr<unordered_map<shared_ptr<abstract_edge>, uint32_t>> &edge_rank_map,
                                                                   const shared_ptr<unordered_map<shared_ptr<abstract_edge>, uint32_t>> &edge_truss_map,
                                                                   const shared_ptr<unordered_set<shared_ptr<abstract_edge>>> &candidate_edge_set,
                                                                   const shared_ptr<unordered_map<shared_ptr<abstract_edge>, uint32_t>> &candidate_edge_support_map,
                                                                   const shared_ptr<unordered_map<shared_ptr<abstract_edge>, uint32_t>> &TS,
                                                                   uint32_t k,
                                                                   const shared_ptr<thread_pool> &pool) {
        uint32_t k_max = 0;
        while (!candidate_edge_set->empty()) {
            k_max = k;

            auto current_edge_set = make_shared<unordered_set<shared_ptr<abstract_edge>>>();
            for (const auto &e: *candidate_edge_set) {
                edge_truss_map->at(e) = k;
                TS->at(e) = candidate_edge_support_map->at(e);

                if (candidate_edge_support_map->at(e) < k - 1) {
                    current_edge_set->insert(e);
                }
            }
            auto evicted_edge_set = make_shared<unordered_set<shared_ptr<abstract_edge>>>();
            remove_unsatisfied_edges(G, edge_mutex_map, edge_rank_map, edge_truss_map, current_edge_set,
                                     candidate_edge_set, candidate_edge_support_map, evicted_edge_set, k + 1, pool);
            ++k;
        }
        return k_max;
    }

    /**
     * @details a top-down maintenance method
     * @param G
     * @param edge_set
     * @param edge_truss_map
     * @param previous_k_max
     * @param thread_number
     */
    void quasi_truss_maintenance::remove(const shared_ptr<abstract_graph> &G,
                                         const shared_ptr<unordered_map<shared_ptr<abstract_edge>, shared_ptr<mutex>>> &edge_mutex_map,
                                         const shared_ptr<unordered_set<shared_ptr<abstract_edge>>> &edge_set,
                                         const shared_ptr<unordered_map<shared_ptr<abstract_edge>, uint32_t>> &edge_rank_map,
                                         const shared_ptr<unordered_map<shared_ptr<abstract_edge>, uint32_t>> &edge_truss_map,
                                         const shared_ptr<unordered_map<shared_ptr<abstract_edge>, uint32_t>> &TS,
                                         const shared_ptr<uint32_t> &previous_k_max,
                                         const shared_ptr<thread_pool> &pool) {
        auto thread_number = pool->get_thread_number();
        auto global_mutex = make_shared<mutex>();

        auto affected_edge_set = make_shared<unordered_set<shared_ptr<abstract_edge>>>();
        auto max_k = make_shared<uint32_t>(1);
        {
            set_edge_rank(edge_set, edge_rank_map, 0);
            auto location_vector = pool->split_task(edge_set);
            for (uint32_t i = 0; i < thread_number; ++i) {
                pool->submit_task([=] {
                    auto sub_affected_set = make_shared<unordered_set<shared_ptr<abstract_edge>>>();
                    uint32_t sub_max_k = 1;

                    auto &sub_begin = *location_vector->at(i);
                    auto &sub_end = *location_vector->at(i + 1);

                    for (auto iter = sub_begin; iter != sub_end; ++iter) {
                        auto &e1 = *iter;
                        auto u = e1->get_source_vertex_id();
                        auto v = e1->get_destination_vertex_id();

                        for (const auto &[w, e2]: *G->get_vertex(u)->get_edge_map()) {
                            if (w == v || edge_rank_map->at(e2) < edge_rank_map->at(e1)) {
                                continue;
                            }

                            auto e3 = G->get_edge(v, w);
                            if (!e3 || edge_rank_map->at(e3) < edge_rank_map->at(e1)) {
                                continue;
                            }

                            if (!edge_set->count(e2)) {
                                sub_affected_set->insert(e2);
                                if (edge_truss_map->at(e2) > sub_max_k) {
                                    sub_max_k = edge_truss_map->at(e2);
                                }
                            }

                            if (!edge_set->count(e3)) {
                                sub_affected_set->insert(e3);
                                if (edge_truss_map->at(e3) > sub_max_k) {
                                    sub_max_k = edge_truss_map->at(e3);
                                }
                            }
                        }
                    }

                    global_mutex->lock();
                    affected_edge_set->merge(*sub_affected_set);
                    if (sub_max_k > *max_k) {
                        *max_k = sub_max_k;
                    }
                    global_mutex->unlock();
                });
            }
            pool->barrier();

            for (const auto &e: *edge_set) {
                edge_rank_map->erase(e);
                G->remove_edge(e);
                edge_truss_map->erase(e);
            }
        }

        /**
         * @brief indicate k_max_flag is updated or not
         */
        auto update_flag = false;
        if (*max_k == *previous_k_max) {
            update_flag = true;
        }

        auto k = static_cast<int64_t>(*max_k);
        while (k >= 1) {
            auto current_removed_edge_set = make_shared<unordered_set<shared_ptr<abstract_edge>>>();
            update_single_truss(G, edge_mutex_map, edge_rank_map, edge_truss_map, affected_edge_set, TS,
                                current_removed_edge_set, k, pool);
            for (const auto &e: *current_removed_edge_set) {
                edge_truss_map->at(e) = k - 1;
                affected_edge_set->insert(e);
            }
            --k;
        }

        if (update_flag) {
            auto current_k_max = make_shared<uint32_t>(0);
            auto location_vector = pool->split_task(edge_truss_map);
            for (uint32_t i = 0; i < thread_number; ++i) {
                pool->submit_task([=] {
                    auto &sub_begin = *location_vector->at(i);
                    auto &sub_end = *location_vector->at(i + 1);

                    uint32_t sub_k_max = 0;

                    for (auto iter = sub_begin; iter != sub_end; ++iter) {
                        auto &[e, e_wing_number] = *iter;
                        if (e_wing_number > sub_k_max) {
                            sub_k_max = e_wing_number;
                        }
                    }

                    global_mutex->lock();
                    if (sub_k_max > *current_k_max) {
                        *current_k_max = sub_k_max;
                    }
                    global_mutex->unlock();
                });
            }
            pool->barrier();
            *previous_k_max = *current_k_max;
        }
    }

    void quasi_truss_maintenance::remove2(const shared_ptr<abstract_graph> &G,
                                          const shared_ptr<unordered_map<shared_ptr<abstract_edge>, shared_ptr<mutex>>> &edge_mutex_map,
                                          const shared_ptr<unordered_set<shared_ptr<abstract_edge>>> &edge_set,
                                          const shared_ptr<unordered_map<shared_ptr<abstract_edge>, uint32_t>> &edge_rank_map,
                                          const shared_ptr<unordered_map<shared_ptr<abstract_edge>, uint32_t>> &edge_truss_map,
                                          const shared_ptr<unordered_map<shared_ptr<abstract_edge>, uint32_t>> &TS,
                                          const shared_ptr<uint32_t> &previous_k_max,
                                          const shared_ptr<thread_pool> &pool) {
        auto thread_number = pool->get_thread_number();
        auto global_mutex = make_shared<mutex>();

        auto affected_edge_set = make_shared<unordered_set<shared_ptr<abstract_edge>>>();
        auto max_k = make_shared<uint32_t>(1);
        {
            set_edge_rank(edge_set, edge_rank_map, 0);
            auto location_vector = pool->split_task(edge_set);
            for (uint32_t i = 0; i < thread_number; ++i) {
                pool->submit_task([=] {
                    auto sub_affected_set = make_shared<unordered_set<shared_ptr<abstract_edge>>>();
                    uint32_t sub_max_k = 1;

                    auto &sub_begin = *location_vector->at(i);
                    auto &sub_end = *location_vector->at(i + 1);

                    for (auto iter = sub_begin; iter != sub_end; ++iter) {
                        auto &e1 = *iter;
                        auto u = e1->get_source_vertex_id();
                        auto v = e1->get_destination_vertex_id();

                        for (const auto &[w, e2]: *G->get_vertex(u)->get_edge_map()) {
                            if (w == v || edge_rank_map->at(e2) < edge_rank_map->at(e1)) {
                                continue;
                            }

                            auto e3 = G->get_edge(v, w);
                            if (!e3 || edge_rank_map->at(e3) < edge_rank_map->at(e1)) {
                                continue;
                            }

                            auto min_truss_number = std::min(
                                    {edge_truss_map->at(e1), edge_truss_map->at(e2), edge_truss_map->at(e3)});

                            if (!edge_set->count(e2) && edge_truss_map->at(e2) == min_truss_number) {
                                edge_mutex_map->at(e2)->lock();
                                --TS->at(e2);
                                edge_mutex_map->at(e2)->unlock();


                                if (TS->at(e2) < edge_truss_map->at(e2)) {
                                    sub_affected_set->insert(e2);
                                }

                                if (edge_truss_map->at(e2) > sub_max_k) {
                                    sub_max_k = edge_truss_map->at(e2);
                                }
                            }


                            if (!edge_set->count(e3) && edge_truss_map->at(e3) == min_truss_number) {
                                edge_mutex_map->at(e3)->lock();
                                --TS->at(e3);
                                edge_mutex_map->at(e3)->unlock();

                                if (TS->at(e3) < edge_truss_map->at(e3)) {
                                    sub_affected_set->insert(e3);
                                }

                                if (edge_truss_map->at(e3) > sub_max_k) {
                                    sub_max_k = edge_truss_map->at(e3);
                                }
                            }
                        }
                    }

                    global_mutex->lock();
                    affected_edge_set->merge(*sub_affected_set);
                    if (sub_max_k > *max_k) {
                        *max_k = sub_max_k;
                    }
                    global_mutex->unlock();
                });
            }
            pool->barrier();

            for (const auto &e: *edge_set) {
                edge_rank_map->erase(e);
                G->remove_edge(e);
                edge_truss_map->erase(e);
            }
        }

        /**
         * @brief indicate k_max_flag is updated or not
         */
        auto update_flag = false;
        if (*max_k == *previous_k_max) {
            update_flag = true;
        }

        auto k = static_cast<int64_t>(*max_k);
        while (k >= 1) {
            auto current_removed_edge_set = make_shared<unordered_set<shared_ptr<abstract_edge>>>();
            update_single_truss2(G, edge_mutex_map, edge_rank_map, edge_truss_map, affected_edge_set, TS,
                                 current_removed_edge_set, k, pool);
            for (const auto &e: *current_removed_edge_set) {
                edge_truss_map->at(e) = k - 1;
                affected_edge_set->insert(e);
            }
            --k;
        }

        if (update_flag) {
            auto current_k_max = make_shared<uint32_t>(0);
            auto location_vector = pool->split_task(edge_truss_map);
            for (uint32_t i = 0; i < thread_number; ++i) {
                pool->submit_task([=] {
                    auto &sub_begin = *location_vector->at(i);
                    auto &sub_end = *location_vector->at(i + 1);

                    uint32_t sub_k_max = 0;

                    for (auto iter = sub_begin; iter != sub_end; ++iter) {
                        auto &[e, e_wing_number] = *iter;
                        if (e_wing_number > sub_k_max) {
                            sub_k_max = e_wing_number;
                        }
                    }

                    global_mutex->lock();
                    if (sub_k_max > *current_k_max) {
                        *current_k_max = sub_k_max;
                    }
                    global_mutex->unlock();
                });
            }
            pool->barrier();
            *previous_k_max = *current_k_max;
        }
    }

    void quasi_truss_maintenance::remove_unsatisfied_edges(const shared_ptr<abstract_graph> &G,
                                                           const shared_ptr<unordered_map<shared_ptr<abstract_edge>, shared_ptr<mutex>>> &edge_mutex_map,
                                                           const shared_ptr<unordered_map<shared_ptr<abstract_edge>, uint32_t>> &edge_rank_map,
                                                           const shared_ptr<unordered_map<shared_ptr<abstract_edge>, uint32_t>> &edge_truss_map,
                                                           const shared_ptr<unordered_set<shared_ptr<abstract_edge>>> &current_edge_set,
                                                           const shared_ptr<unordered_set<shared_ptr<abstract_edge>>> &candidate_edge_set,
                                                           const shared_ptr<unordered_map<shared_ptr<abstract_edge>, uint32_t>> &candidate_edge_support_map,
                                                           const shared_ptr<unordered_set<shared_ptr<abstract_edge>>> &evicted_edge_set,
                                                           uint32_t k,
                                                           const shared_ptr<thread_pool> &pool) {
        auto thread_number = pool->get_thread_number();
        auto global_mutex = make_shared<mutex>();

        auto current_evicted_set = make_shared<unordered_set<shared_ptr<abstract_edge>>>();

        auto rand_id = make_shared<uint32_t>(0);
        while (!current_edge_set->empty()) {
            auto next_edge_set = make_shared<unordered_set<shared_ptr<abstract_edge>>>();
            set_edge_rank(current_edge_set, edge_rank_map, rand_id);
            auto location_vector = pool->split_task(current_edge_set);
            for (uint32_t i = 0; i < thread_number; ++i) {
                pool->submit_task([=] {
                    auto sub_next_edge_set = make_shared<unordered_set<shared_ptr<abstract_edge>>>();

                    auto &sub_begin = *location_vector->at(i);
                    auto &sub_end = *location_vector->at(i + 1);

                    for (auto iter = sub_begin; iter != sub_end; ++iter) {
                        auto &e1 = *iter;

                        auto u = e1->get_source_vertex_id();
                        auto v = e1->get_destination_vertex_id();
                        for (const auto &[w, e2]: *G->get_vertex(u)->get_edge_map()) {
                            if (w == v || edge_truss_map->at(e2) < k - 1 ||
                                edge_rank_map->at(e2) < edge_rank_map->at(e1) ||
                                (edge_truss_map->at(e2) == k - 1 && !candidate_edge_set->count(e2))) {
                                continue;
                            }

                            auto e3 = G->get_edge(v, w);
                            if (!e3 || edge_truss_map->at(e3) < k - 1 ||
                                edge_rank_map->at(e3) < edge_rank_map->at(e1) ||
                                (edge_truss_map->at(e3) == k - 1 && !candidate_edge_set->count(e3))) {
                                continue;
                            }

                            if (edge_truss_map->at(e2) == k - 1 && candidate_edge_support_map->at(e2) >= k - 2) {
                                edge_mutex_map->at(e2)->lock();
                                --candidate_edge_support_map->at(e2);
                                edge_mutex_map->at(e2)->unlock();
                                if (candidate_edge_support_map->at(e2) < k && !current_edge_set->count(e2)) {
                                    sub_next_edge_set->insert(e2);
                                }
                            }

                            if (edge_truss_map->at(e3) == k - 1 && candidate_edge_support_map->at(e3) >= k - 2) {
                                edge_mutex_map->at(e3)->lock();
                                --candidate_edge_support_map->at(e3);
                                edge_mutex_map->at(e3)->unlock();
                                if (candidate_edge_support_map->at(e3) < k && !current_edge_set->count(e3)) {
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
            pool->barrier();
            for (const auto &e: *current_edge_set) {
                candidate_edge_set->erase(e);
                candidate_edge_support_map->at(e) = 0;
            }
            current_evicted_set->merge(*current_edge_set);
            swap(*current_edge_set, *next_edge_set);
        }

        for (const auto &e: *current_evicted_set) {
            edge_rank_map->at(e) = UINT32_MAX;
        }
        evicted_edge_set->merge(*current_evicted_set);
    }

    void
    quasi_truss_maintenance::set_edge_rank(const shared_ptr<unordered_set<shared_ptr<scnu::abstract_edge>>> &edge_set,
                                           const shared_ptr<unordered_map<shared_ptr<scnu::abstract_edge>, uint32_t>> &edge_rank_map,
                                           const shared_ptr<uint32_t> &rank_id) {
        for (const auto &e: *edge_set) {
            *rank_id = *rank_id + 1;
            edge_rank_map->at(e) = *rank_id;
        }
    }

    void quasi_truss_maintenance::update_edge_truss_support(const shared_ptr<abstract_graph> &G,
                                                            const shared_ptr<unordered_map<shared_ptr<abstract_edge>, shared_ptr<mutex>>> &edge_mutex_map,
                                                            const shared_ptr<unordered_map<shared_ptr<abstract_edge>, uint32_t>> &edge_rank_map,
                                                            const shared_ptr<unordered_map<shared_ptr<abstract_edge>, uint32_t>> &edge_truss_map,
                                                            const shared_ptr<unordered_map<shared_ptr<abstract_edge>, uint32_t>> &TS,
                                                            const shared_ptr<unordered_set<shared_ptr<abstract_edge>>> &candidate_edge_set,
                                                            uint32_t k,
                                                            const shared_ptr<thread_pool> &pool) {
        set_edge_rank(candidate_edge_set, edge_rank_map, 0);
        auto thread_number = pool->get_thread_number();
        auto global_mutex = make_shared<mutex>();
        auto location_vector = pool->split_task(candidate_edge_set);
        for (uint32_t i = 0; i < thread_number; ++i) {
            pool->submit_task([=] {
                auto &sub_begin = *location_vector->at(i);
                auto &sub_end = *location_vector->at(i + 1);

                for (auto iter = sub_begin; iter != sub_end; ++iter) {
                    auto &e1 = *iter;
                    auto u = e1->get_source_vertex_id();
                    auto v = e1->get_destination_vertex_id();

                    for (const auto &[w, e2]: *G->get_vertex(u)->get_edge_map()) {
                        if (w == v || edge_truss_map->at(e2) < k || edge_rank_map->at(e2) < edge_rank_map->at(e1)) {
                            continue;
                        }

                        auto e3 = G->get_edge(v, w);
                        if (!e3 || edge_truss_map->at(e3) < k || edge_rank_map->at(e3) < edge_rank_map->at(e1)) {
                            continue;
                        }

                        if (edge_truss_map->at(e2) == k && !candidate_edge_set->count(e2)) {
                            edge_mutex_map->at(e2)->lock();
                            ++TS->at(e2);
                            edge_mutex_map->at(e2)->unlock();
                        }

                        if (edge_truss_map->at(e3) == k && !candidate_edge_set->count(e3)) {
                            edge_mutex_map->at(e3)->lock();
                            ++TS->at(e3);
                            edge_mutex_map->at(e3)->unlock();
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

    void quasi_truss_maintenance::update_single_truss(const shared_ptr<abstract_graph> &G,
                                                      const shared_ptr<unordered_map<shared_ptr<abstract_edge>, shared_ptr<mutex>>> &edge_mutex_map,
                                                      const shared_ptr<unordered_map<shared_ptr<abstract_edge>, uint32_t>> &edge_rank_map,
                                                      const shared_ptr<unordered_map<shared_ptr<abstract_edge>, uint32_t>> &edge_truss_map,
                                                      const shared_ptr<unordered_set<shared_ptr<abstract_edge>>> &edge_set,
                                                      const shared_ptr<unordered_map<shared_ptr<abstract_edge>, uint32_t>> &TS,
                                                      const shared_ptr<unordered_set<shared_ptr<abstract_edge>>> &evicted_edge_set,
                                                      uint32_t k,
                                                      const shared_ptr<thread_pool> &pool) {
        auto thread_number = pool->get_thread_number();
        auto global_mutex = make_shared<mutex>();

        auto current_edge_set = make_shared<unordered_set<shared_ptr<abstract_edge>>>();
        {
            auto removed_edge_set = make_shared<unordered_set<shared_ptr<abstract_edge>>>();
            auto location_vector = pool->split_task(edge_set);
            for (uint32_t i = 0; i < thread_number; ++i) {
                pool->submit_task([=] {
                    auto sub_current_edge_set = make_shared<unordered_set<shared_ptr<abstract_edge>>>();
                    auto sub_removed_edge_set = make_shared<unordered_set<shared_ptr<abstract_edge>>>();

                    auto &sub_begin = *location_vector->at(i);
                    auto &sub_end = *location_vector->at(i + 1);

                    for (auto iter = sub_begin; iter != sub_end; ++iter) {
                        auto &e = *iter;
                        if (edge_truss_map->at(e) > k) {
                            sub_removed_edge_set->insert(e);
                        } else {
                            if (edge_truss_map->at(e) == k) {
                                TS->at(e) = edge_support_computation(G, e, edge_truss_map, k);

                                if (TS->at(e) < k - 2) {
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
            for (const auto &e: *removed_edge_set) {
                edge_set->erase(e);
            }
        }

        auto clear_set = make_shared<unordered_set<shared_ptr<abstract_edge>>>();
        /**
         * @brief a set of edges revising their rank
         */
        auto rand_id = make_shared<uint32_t>(0);
        while (!current_edge_set->empty()) {
            set_edge_rank(current_edge_set, edge_rank_map, rand_id);
            auto next_edge_set = make_shared<unordered_set<shared_ptr<abstract_edge>>>();
            auto location_vector = pool->split_task(current_edge_set);
            for (uint32_t i = 0; i < thread_number; ++i) {
                pool->submit_task([=] {
                    auto sub_next_edge_set = make_shared<unordered_set<shared_ptr<abstract_edge>>>();
                    auto sub_clear_set = make_shared<unordered_set<shared_ptr<abstract_edge>>>();

                    auto &sub_begin = *location_vector->at(i);
                    auto &sub_end = *location_vector->at(i + 1);

                    for (auto iter = sub_begin; iter != sub_end; ++iter) {
                        auto &e1 = *iter;

                        auto u = e1->get_source_vertex_id();
                        auto v = e1->get_destination_vertex_id();
                        for (const auto &[r2, e2]: *G->get_vertex(u)->get_edge_map()) {
                            if (r2 == v || edge_truss_map->at(e2) < k ||
                                edge_rank_map->at(e2) < edge_rank_map->at(e1)) {
                                continue;
                            }
                            auto e3 = G->get_edge(v, r2);
                            if (!e3 || edge_truss_map->at(e3) < k || edge_rank_map->at(e3) < edge_rank_map->at(e1)) {
                                continue;
                            }

                            if (!current_edge_set->count(e2) && edge_truss_map->at(e2) == k) {
                                edge_mutex_map->at(e2)->lock();
                                if (TS->at(e2) == 0) {
                                    sub_clear_set->insert(e2);
                                    TS->at(e2) = edge_support_computation(G, e2, edge_truss_map, k);
                                }
                                --TS->at(e2);
                                edge_mutex_map->at(e2)->unlock();

                                if (TS->at(e2) < k - 2) {
                                    sub_next_edge_set->insert(e2);
                                }
                            }


                            if (!current_edge_set->count(e3) && edge_truss_map->at(e3) == k) {
                                edge_mutex_map->at(e3)->lock();
                                if (TS->at(e3) == 0) {
                                    sub_clear_set->insert(e3);
                                    TS->at(e3) = edge_support_computation(G, e3, edge_truss_map, k);
                                }
                                --TS->at(e3);
                                edge_mutex_map->at(e3)->unlock();

                                if (TS->at(e3) < k - 2) {
                                    current_edge_set->insert(e3);
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

        for (const auto &e: *clear_set) {
            TS->at(e) = 0;
        }
        for (const auto &e: *evicted_edge_set) {
            edge_rank_map->at(e) = UINT32_MAX;
        }
    }

    void quasi_truss_maintenance::update_single_truss2(const shared_ptr<abstract_graph> &G,
                                                       const shared_ptr<unordered_map<shared_ptr<abstract_edge>, shared_ptr<mutex>>> &edge_mutex_map,
                                                       const shared_ptr<unordered_map<shared_ptr<abstract_edge>, uint32_t>> &edge_rank_map,
                                                       const shared_ptr<unordered_map<shared_ptr<abstract_edge>, uint32_t>> &edge_truss_map,
                                                       const shared_ptr<unordered_set<shared_ptr<abstract_edge>>> &edge_set,
                                                       const shared_ptr<unordered_map<shared_ptr<abstract_edge>, uint32_t>> &TS,
                                                       const shared_ptr<unordered_set<shared_ptr<abstract_edge>>> &evicted_edge_set,
                                                       uint32_t k,
                                                       const shared_ptr<thread_pool> &pool) {
        auto thread_number = pool->get_thread_number();
        auto global_mutex = make_shared<mutex>();

        auto current_edge_set = make_shared<unordered_set<shared_ptr<abstract_edge>>>();
        {
            auto removed_edge_set = make_shared<unordered_set<shared_ptr<abstract_edge>>>();
            auto location_vector = pool->split_task(edge_set);
            for (uint32_t i = 0; i < thread_number; ++i) {
                pool->submit_task([=] {
                    auto sub_current_edge_set = make_shared<unordered_set<shared_ptr<abstract_edge>>>();
                    auto sub_removed_edge_set = make_shared<unordered_set<shared_ptr<abstract_edge>>>();

                    auto &sub_begin = *location_vector->at(i);
                    auto &sub_end = *location_vector->at(i + 1);

                    for (auto iter = sub_begin; iter != sub_end; ++iter) {
                        auto &e = *iter;
                        if (edge_truss_map->at(e) > k) {
                            sub_removed_edge_set->insert(e);
                        } else {
                            if (edge_truss_map->at(e) == k) {
                                TS->at(e) = edge_support_computation(G, e, edge_truss_map, k);

                                if (TS->at(e) < k - 2) {
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
            for (const auto &e: *removed_edge_set) {
                edge_set->erase(e);
            }
        }
        /**
         * @brief a set of edges revising their rank
         */
        auto rand_id = make_shared<uint32_t>(0);
        while (!current_edge_set->empty()) {
            set_edge_rank(current_edge_set, edge_rank_map, rand_id);
            auto next_edge_set = make_shared<unordered_set<shared_ptr<abstract_edge>>>();
            auto location_vector = pool->split_task(current_edge_set);
            for (uint32_t i = 0; i < thread_number; ++i) {
                pool->submit_task([=] {
                    auto sub_next_edge_set = make_shared<unordered_set<shared_ptr<abstract_edge>>>();

                    auto &sub_begin = *location_vector->at(i);
                    auto &sub_end = *location_vector->at(i + 1);

                    for (auto iter = sub_begin; iter != sub_end; ++iter) {
                        auto &e1 = *iter;

                        auto u = e1->get_source_vertex_id();
                        auto v = e1->get_destination_vertex_id();
                        for (const auto &[w, e2]: *G->get_vertex(u)->get_edge_map()) {
                            if (w == v || edge_truss_map->at(e2) < k ||
                                edge_rank_map->at(e2) < edge_rank_map->at(e1)) {
                                continue;
                            }

                            auto e3 = G->get_edge(v, w);
                            if (!e3 || edge_truss_map->at(e3) < k || edge_rank_map->at(e3) < edge_rank_map->at(e1)) {
                                continue;
                            }

                            if (!current_edge_set->count(e2) && edge_truss_map->at(e2) == k && TS->at(e2) >= k - 2) {
                                edge_mutex_map->at(e2)->lock();
                                --TS->at(e2);
                                edge_mutex_map->at(e2)->unlock();

                                if (TS->at(e2) < k) {
                                    sub_next_edge_set->insert(e2);
                                }
                            }


                            if (!current_edge_set->count(e3) && edge_truss_map->at(e3) == k && TS->at(e3) >= k - 2) {
                                edge_mutex_map->at(e3)->lock();
                                --TS->at(e3);
                                edge_mutex_map->at(e3)->unlock();

                                if (TS->at(e3) < k) {
                                    current_edge_set->insert(e3);
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

        for (const auto &e: *evicted_edge_set) {
            edge_rank_map->at(e) = UINT32_MAX;
        }
    }
}





