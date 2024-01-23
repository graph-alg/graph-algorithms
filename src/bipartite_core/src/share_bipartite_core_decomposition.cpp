
#include "bipartite_core/share_bipartite_core_decomposition.h"

namespace scnu {

    void share_bipartite_core_decomposition::init(const shared_ptr<abstract_bipartite_graph> &B,
                                                  const shared_ptr<unordered_map<uint32_t, shared_ptr<mutex>>> &left_mutex_map,
                                                  const shared_ptr<unordered_map<uint32_t, shared_ptr<mutex>>> &right_mutex_map,
                                                  const shared_ptr<unordered_map<uint32_t, shared_ptr<bipartite_core_left_store_index>>> &left_index_map,
                                                  const shared_ptr<unordered_map<uint32_t, shared_ptr<bipartite_core_right_store_index>>> &right_index_map,
                                                  const shared_ptr<scnu::thread_pool> &pool) {
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

    void share_bipartite_core_decomposition::compute_alpha_core(const shared_ptr<abstract_bipartite_graph> &B,
                                                                const shared_ptr<unordered_map<uint32_t, shared_ptr<bipartite_core_left_store_index>>> &left_index_map,
                                                                const shared_ptr<unordered_map<uint32_t, shared_ptr<bipartite_core_right_store_index>>> &right_index_map,
                                                                uint32_t alpha) {
        auto left_degree_map = make_shared<unordered_map<uint32_t, uint32_t>>();
        auto right_degree_map = make_shared<unordered_map<uint32_t, uint32_t>>();

        for (const auto &[r, r_vertex]: *B->get_right_vertex_map()) {
            right_degree_map->insert({r, r_vertex->get_degree()});
        }

        if (alpha == 1) {
            for (const auto &[l, l_vertex]: *B->get_left_vertex_map()) {
                left_degree_map->insert({l, l_vertex->get_degree()});
            }
        } else {
            for (const auto &[l, l_vertex]: *B->get_left_vertex_map()) {
                if (l_vertex->get_degree() < alpha) {
                    for (const auto &[r, e]: *l_vertex->get_edge_map()) {
                        --right_degree_map->at(r);
                        if (right_degree_map->at(r) == 0) {
                            right_degree_map->erase(r);
                        }
                    }
                } else {
                    left_degree_map->insert({l, l_vertex->get_degree()});
                }
            }
        }

        auto evicted_r_set = make_shared<unordered_set<uint32_t> >();

        while (!left_degree_map->empty() && !right_degree_map->empty()) {
            uint32_t beta = right_degree_map->begin()->second;

            for (const auto&[r,r_degree]:*right_degree_map) {
                if (r_degree < beta) {
                    beta = r_degree;

                    evicted_r_set->clear();
                    evicted_r_set->insert(r);
                } else if(r_degree == beta){
                    evicted_r_set->insert(r);
                }
            }

            while (!evicted_r_set->empty()) {
                auto r = *evicted_r_set->begin();
                evicted_r_set->erase(r);

                right_degree_map->erase(r);

                for (auto j = 1; j <= beta; ++j) {
                    right_index_map->at(r)->insert(j, alpha);
                }

                auto r_vertex = B->get_right_vertex(r);
                for (const auto &[l, e]: *r_vertex->get_edge_map()) {
                    if (left_degree_map->count(l) && left_degree_map->at(l) >= alpha) {
                        --left_degree_map->at(l);
                        if (left_degree_map->at(l) < alpha) {

                            left_degree_map->erase(l);
                            left_index_map->at(l)->insert(alpha, beta);

                            for (const auto &[r2, e2]: *B->get_left_vertex(l)->get_edge_map()) {
                                if (right_degree_map->count(r2) && right_degree_map->at(r2) > beta) {
                                    --right_degree_map->at(r2);
                                    if (right_degree_map->at(r2) == beta) {
                                        evicted_r_set->insert(r2);
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
    }


    void share_bipartite_core_decomposition::compute_alpha_core(const shared_ptr<abstract_bipartite_graph> &B,
                                                                const shared_ptr<unordered_map<uint32_t, shared_ptr<mutex>>> &left_mutex_map,
                                                                const shared_ptr<unordered_map<uint32_t, shared_ptr<mutex>>> &right_mutex_map,
                                                                const shared_ptr<unordered_map<uint32_t, shared_ptr<bipartite_core_left_store_index>>> &left_index_map,
                                                                const shared_ptr<unordered_map<uint32_t, shared_ptr<bipartite_core_right_store_index>>> &right_index_map,
                                                                uint32_t alpha) {
        auto vertex_degree_map = make_shared<unordered_map<uint32_t, uint32_t>>();
        auto l_set = make_shared<unordered_set<uint32_t>>();
        auto r_set = make_shared<unordered_set<uint32_t>>();

        for (const auto &[r, r_vertex]: *B->get_right_vertex_map()) {
            r_set->insert(r);
            vertex_degree_map->insert({r, r_vertex->get_degree()});
        }

        if (alpha == 1) {
            for (const auto &[l, l_vertex]: *B->get_left_vertex_map()) {
                l_set->insert(l);
                vertex_degree_map->insert({l, l_vertex->get_degree()});
            }
        } else {
            for (const auto &[l, l_vertex]: *B->get_left_vertex_map()) {
                if (l_vertex->get_degree() < alpha) {
                    for (const auto &[r, e]: *l_vertex->get_edge_map()) {
                        --vertex_degree_map->at(r);
                        if (vertex_degree_map->at(r) == 0) {
                            r_set->erase(r);
                            vertex_degree_map->erase(r);
                        }
                    }
                } else {
                    l_set->insert(l);
                    vertex_degree_map->insert({l, l_vertex->get_degree()});
                }
            }
        }
        auto evicted_r_set = make_shared<unordered_set<uint32_t> >();

        while (!r_set->empty()) {
            uint32_t beta = UINT32_MAX;

            for (const auto &r: *r_set) {
                auto r_degree = vertex_degree_map->at(r);
                if (r_degree < beta) {
                    beta = r_degree;

                    evicted_r_set->clear();
                    evicted_r_set->insert(r);
                } else if (r_degree == beta) {
                    evicted_r_set->insert(r);
                }
            }

            while (!evicted_r_set->empty()) {
                auto r = *evicted_r_set->begin();
                evicted_r_set->erase(r);

                r_set->erase(r);
                vertex_degree_map->erase(r);

                right_mutex_map->at(r)->lock();
                for (auto j = 1; j <= beta; ++j) {
                    right_index_map->at(r)->insert(j, alpha);
                }
                right_mutex_map->at(r)->unlock();

                for (const auto &[l, e]: *B->get_right_vertex(r)->get_edge_map()) {
                    if (l_set->count(l)) {
                        --vertex_degree_map->at(l);
                        if (vertex_degree_map->at(l) < alpha) {
                            l_set->erase(l);
                            vertex_degree_map->erase(l);

                            left_mutex_map->at(l)->lock();
                            left_index_map->at(l)->insert(alpha, beta);
                            left_mutex_map->at(l)->unlock();

                            if (alpha > 1) {
                                for (const auto &[r2, e2]: *B->get_left_vertex(l)->get_edge_map()) {
                                    if (r_set->count(r2)) {
                                        --vertex_degree_map->at(r2);
                                        if (vertex_degree_map->at(r2) <= beta) {
                                            evicted_r_set->insert(r2);
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
    }


    void share_bipartite_core_decomposition::compute_alpha_core2(const shared_ptr<abstract_bipartite_graph> &B,
                                                                 const shared_ptr<unordered_map<uint32_t, shared_ptr<mutex>>> &left_mutex_map,
                                                                 const shared_ptr<unordered_map<uint32_t, shared_ptr<mutex>>> &right_mutex_map,
                                                                 const shared_ptr<unordered_map<uint32_t, shared_ptr<bipartite_core_left_store_index>>> &left_index_map,
                                                                 const shared_ptr<unordered_map<uint32_t, shared_ptr<bipartite_core_right_store_index>>> &right_index_map,
                                                                 uint32_t alpha) {
        auto left_degree_map = make_shared<unordered_map<uint32_t, uint32_t>>();
        auto right_degree_map = make_shared<unordered_map<uint32_t, uint32_t>>();

        for (const auto &[r, r_vertex]: *B->get_right_vertex_map()) {
            right_degree_map->insert({r, r_vertex->get_degree()});
        }

        if (alpha == 1) {
            for (const auto &[l, l_vertex]: *B->get_left_vertex_map()) {
                left_degree_map->insert({l, l_vertex->get_degree()});
            }
        } else {
            for (const auto &[l, l_vertex]: *B->get_left_vertex_map()) {
                if (l_vertex->get_degree() < alpha) {
                    for (const auto &[r, e]: *l_vertex->get_edge_map()) {
                        --right_degree_map->at(r);
                        if (right_degree_map->at(r) == 0) {
                            right_degree_map->erase(r);
                        }
                    }
                } else {
                    left_degree_map->insert({l, l_vertex->get_degree()});
                }
            }
        }
        auto evicted_r_set = make_shared<unordered_set<uint32_t> >();

        while (!right_degree_map->empty()) {
            uint32_t beta = right_degree_map->begin()->second;

            for (const auto &[r, r_degree]: *right_degree_map) {
                if (r_degree < beta) {
                    beta = r_degree;

                    evicted_r_set->clear();
                    evicted_r_set->insert(r);
                } else if (r_degree == beta) {
                    evicted_r_set->insert(r);
                }
            }

            while (!evicted_r_set->empty()) {
                auto r = *evicted_r_set->begin();
                evicted_r_set->erase(r);
                right_degree_map->erase(r);

                right_mutex_map->at(r)->lock();
                for (auto j = 1; j <= beta; ++j) {
                    right_index_map->at(r)->insert(j, alpha);
                }
                right_mutex_map->at(r)->unlock();

                for (const auto &[l, e]: *B->get_right_vertex(r)->get_edge_map()) {
                    if (left_degree_map->count(l) && left_degree_map->at(l) >= alpha) {
                        --left_degree_map->at(l);
                        if (left_degree_map->at(l) < alpha) {

                            left_degree_map->erase(l);

                            left_mutex_map->at(l)->lock();
                            left_index_map->at(l)->insert(alpha, beta);
                            left_mutex_map->at(l)->unlock();

                            if (alpha > 1) {
                                for (const auto &[r2, e2]: *B->get_left_vertex(l)->get_edge_map()) {
                                    if (right_degree_map->count(r2) && right_degree_map->at(r2) > beta) {
                                        --right_degree_map->at(r2);
                                        if (right_degree_map->at(r2) == beta) {
                                            evicted_r_set->insert(r2);
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
    }

    void share_bipartite_core_decomposition::compute_beta_core(const shared_ptr<abstract_bipartite_graph> &B,
                                                               const shared_ptr<unordered_map<uint32_t, shared_ptr<bipartite_core_left_store_index>>> &left_index_map,
                                                               const shared_ptr<unordered_map<uint32_t, shared_ptr<bipartite_core_right_store_index>>> &right_index_map,
                                                               uint32_t beta) {
        auto left_degree_map = make_shared<unordered_map<uint32_t, uint32_t>>();
        auto right_degree_map = make_shared<unordered_map<uint32_t, uint32_t>>();

        for (const auto &[l, l_vertex]: *B->get_left_vertex_map()) {
            left_degree_map->insert({l, l_vertex->get_degree()});
        }

        if (beta == 1) {
            for (const auto &[r, r_vertex]: *B->get_right_vertex_map()) {
                right_degree_map->insert({r, r_vertex->get_degree()});
            }
        } else {
            for (const auto &[r, r_vertex]: *B->get_right_vertex_map()) {
                if (r_vertex->get_degree() < beta) {
                    for (const auto &[l, e]: *r_vertex->get_edge_map()) {
                        --left_degree_map->at(l);
                        if (left_degree_map->at(l) == 0) {
                            left_degree_map->erase(l);
                        }
                    }
                } else {
                    right_degree_map->insert({r, r_vertex->get_degree()});
                }
            }
        }

        auto evicted_l_set = make_shared<unordered_set<uint32_t> >();

        while (!left_degree_map->empty()) {
            uint32_t alpha = left_degree_map->begin()->second;
            for (const auto &[l, l_degree]: *left_degree_map) {
                if (l_degree < alpha) {
                    alpha = l_degree;

                    evicted_l_set->clear();
                    evicted_l_set->insert(l);
                } else if (l_degree == alpha) {
                    evicted_l_set->insert(l);
                }
            }

            while (!evicted_l_set->empty()) {
                auto l = *evicted_l_set->begin();
                evicted_l_set->erase(l);
                left_degree_map->erase(l);

                for (auto i = 1; i <= alpha; ++i) {
                    left_index_map->at(l)->insert(i, beta);
                }

                auto l_vertex = B->get_left_vertex(l);
                for (const auto &[r, e]: *l_vertex->get_edge_map()) {
                    if (right_degree_map->count(r) && right_degree_map->at(r) >= beta) {
                        --right_degree_map->at(r);
                        if (right_degree_map->at(r) < beta) {

                            right_degree_map->erase(r);
                            right_index_map->at(r)->insert(beta, alpha);

                            for (const auto &[l2, e2]: *B->get_right_vertex(r)->get_edge_map()) {
                                if (left_degree_map->count(l2) && left_degree_map->at(l2) > alpha) {
                                    --left_degree_map->at(l2);
                                    if (left_degree_map->at(l2) == alpha) {
                                        evicted_l_set->insert(l2);
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
    }

    void share_bipartite_core_decomposition::compute_beta_core(const shared_ptr<abstract_bipartite_graph> &B,
                                                               const shared_ptr<unordered_map<uint32_t, shared_ptr<mutex>>> &left_mutex_map,
                                                               const shared_ptr<unordered_map<uint32_t, shared_ptr<mutex>>> &right_mutex_map,
                                                               const shared_ptr<unordered_map<uint32_t, shared_ptr<bipartite_core_left_store_index>>> &left_index_map,
                                                               const shared_ptr<unordered_map<uint32_t, shared_ptr<bipartite_core_right_store_index>>> &right_index_map,
                                                               uint32_t beta) {
        auto vertex_degree_map = make_shared<unordered_map<uint32_t, uint32_t>>();

        auto l_set = make_shared<unordered_set<uint32_t>>();
        auto r_set = make_shared<unordered_set<uint32_t>>();

        for (const auto &[l, l_vertex]: *B->get_left_vertex_map()) {
            l_set->insert(l);
            vertex_degree_map->insert({l, l_vertex->get_degree()});
        }

        if (beta == 1) {
            for (const auto &[r, r_vertex]: *B->get_right_vertex_map()) {
                r_set->insert(r);
                vertex_degree_map->insert({r, r_vertex->get_degree()});
            }
        } else {
            for (const auto &[r, r_vertex]: *B->get_right_vertex_map()) {
                if (r_vertex->get_degree() < beta) {
                    for (const auto &[l, e]: *r_vertex->get_edge_map()) {
                        --vertex_degree_map->at(l);
                        if (vertex_degree_map->at(l) == 0) {
                            l_set->erase(l);
                            vertex_degree_map->erase(l);
                        }
                    }
                } else {
                    r_set->insert(r);
                    vertex_degree_map->insert({r, r_vertex->get_degree()});
                }
            }
        }

        auto evicted_l_set = make_shared<unordered_set<uint32_t> >();

        while (!l_set->empty()) {
            uint32_t alpha = UINT32_MAX;
            for (const auto &l: *l_set) {
                auto l_degree = vertex_degree_map->at(l);
                if (l_degree < alpha) {
                    alpha = l_degree;

                    evicted_l_set->clear();
                    evicted_l_set->insert(l);
                } else if (l_degree == alpha) {
                    evicted_l_set->insert(l);
                }
            }

            while (!evicted_l_set->empty()) {
                auto l = *evicted_l_set->begin();
                evicted_l_set->erase(l);

                l_set->erase(l);
                vertex_degree_map->erase(l);

                left_mutex_map->at(l)->lock();
                for (auto i = 1; i <= alpha; ++i) {
                    left_index_map->at(l)->insert(i, beta);
                }
                left_mutex_map->at(l)->unlock();

                for (const auto &[r, e]: *B->get_left_vertex(l)->get_edge_map()) {
                    if (vertex_degree_map->count(r)) {
                        --vertex_degree_map->at(r);
                        if (vertex_degree_map->at(r) < beta) {
                            r_set->erase(r);
                            vertex_degree_map->erase(r);

                            right_mutex_map->at(r)->lock();
                            right_index_map->at(r)->insert(beta, alpha);
                            right_mutex_map->at(r)->unlock();

                            if (beta > 1) {
                                for (const auto &[l2, e2]: *B->get_right_vertex(r)->get_edge_map()) {
                                    if (l_set->count(l2)) {
                                        --vertex_degree_map->at(l2);
                                        if (vertex_degree_map->at(l2) <= alpha) {
                                            evicted_l_set->insert(l2);
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
    }

    void share_bipartite_core_decomposition::compute_beta_core2(const shared_ptr<abstract_bipartite_graph> &B,
                                                                const shared_ptr<unordered_map<uint32_t, shared_ptr<mutex>>> &left_mutex_map,
                                                                const shared_ptr<unordered_map<uint32_t, shared_ptr<mutex>>> &right_mutex_map,
                                                                const shared_ptr<unordered_map<uint32_t, shared_ptr<bipartite_core_left_store_index>>> &left_index_map,
                                                                const shared_ptr<unordered_map<uint32_t, shared_ptr<bipartite_core_right_store_index>>> &right_index_map,
                                                                uint32_t beta) {

        auto left_degree_map = make_shared<unordered_map<uint32_t, uint32_t>>();
        auto right_degree_map = make_shared<unordered_map<uint32_t, uint32_t>>();

        for (const auto &[l, l_vertex]: *B->get_left_vertex_map()) {
            left_degree_map->insert({l, l_vertex->get_degree()});
        }

        if (beta == 1) {
            for (const auto &[r, r_vertex]: *B->get_right_vertex_map()) {
                right_degree_map->insert({r, r_vertex->get_degree()});
            }
        } else {
            for (const auto &[r, r_vertex]: *B->get_right_vertex_map()) {
                if (r_vertex->get_degree() < beta) {
                    for (const auto &[l, e]: *r_vertex->get_edge_map()) {
                        --left_degree_map->at(l);
                        if (left_degree_map->at(l) == 0) {
                            left_degree_map->erase(l);
                        }
                    }
                } else {
                    right_degree_map->insert({r, r_vertex->get_degree()});
                }
            }
        }

        auto evicted_l_set = make_shared<unordered_set<uint32_t> >();

        while (!left_degree_map->empty()) {
            uint32_t alpha = left_degree_map->begin()->second;
            for (const auto &[l, l_degree]: *left_degree_map) {
                if (l_degree < alpha) {
                    alpha = l_degree;

                    evicted_l_set->clear();
                    evicted_l_set->insert(l);
                } else if (l_degree == alpha) {
                    evicted_l_set->insert(l);
                }
            }

            while (!evicted_l_set->empty()) {
                auto l = *evicted_l_set->begin();
                evicted_l_set->erase(l);
                left_degree_map->erase(l);

                left_mutex_map->at(l)->lock();
                for (auto i = 1; i <= alpha; ++i) {
                    left_index_map->at(l)->insert(i, beta);
                }
                left_mutex_map->at(l)->unlock();

                for (const auto &[r, e]: *B->get_left_vertex(l)->get_edge_map()) {
                    if (right_degree_map->count(r) && right_degree_map->at(r) >= beta) {
                        --right_degree_map->at(r);
                        if (right_degree_map->at(r) < beta) {

                            right_degree_map->erase(r);
                            right_mutex_map->at(r)->lock();
                            right_index_map->at(r)->insert(beta, alpha);
                            right_mutex_map->at(r)->unlock();

                            if (beta > 1) {
                                for (const auto &[l2, e2]: *B->get_right_vertex(r)->get_edge_map()) {
                                    if (left_degree_map->count(l2) && left_degree_map->at(l2) > alpha) {
                                        --left_degree_map->at(l2);
                                        if (left_degree_map->at(l2) == alpha) {
                                            evicted_l_set->insert(l2);
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
    }

    uint32_t share_bipartite_core_decomposition::decompose(const shared_ptr<abstract_bipartite_graph> &B,
                                                           const shared_ptr<unordered_map<uint32_t, shared_ptr<bipartite_core_left_store_index>>> &left_index_map,
                                                           const shared_ptr<unordered_map<uint32_t, shared_ptr<bipartite_core_right_store_index>>> &right_index_map) {
        for (const auto &[l, l_vertex]: *B->get_left_vertex_map()) {
            left_index_map->insert({l, make_shared<bipartite_core_left_store_index>()});
        }
        for (const auto &[r, r_vertex]: *B->get_right_vertex_map()) {
            right_index_map->insert({r, make_shared<bipartite_core_right_store_index>()});
        }

        uint32_t max_delta = find_core2(B);

        for (auto delta = 1; delta <= max_delta; ++delta) {
            compute_alpha_core(B,  left_index_map, right_index_map, delta);
            compute_beta_core(B, left_index_map, right_index_map, delta);
        }
        return max_delta;
    }



    uint32_t share_bipartite_core_decomposition::decompose(const shared_ptr<abstract_bipartite_graph> &B,
                                                           const shared_ptr<unordered_map<uint32_t, shared_ptr<mutex>>> &left_mutex_map,
                                                           const shared_ptr<unordered_map<uint32_t, shared_ptr<mutex>>> &right_mutex_map,
                                                           const shared_ptr<unordered_map<uint32_t, shared_ptr<bipartite_core_left_store_index>>> &left_index_map,
                                                           const shared_ptr<unordered_map<uint32_t, shared_ptr<bipartite_core_right_store_index>>> &right_index_map,
                                                           const shared_ptr<thread_pool> &pool) {


        auto max_delta = find_core(B);

        for (auto delta = 1; delta <= max_delta; ++delta) {
            pool->submit_task([=] {
                compute_alpha_core(B, left_mutex_map, right_mutex_map, left_index_map, right_index_map, delta);
            });
        }

        for (auto delta = 1; delta <= max_delta; ++delta) {
            pool->submit_task([=] {
                compute_beta_core(B, left_mutex_map, right_mutex_map, left_index_map, right_index_map, delta);
            });
        }
        pool->barrier();

        return max_delta;
    }

    uint32_t share_bipartite_core_decomposition::decompose2(const shared_ptr<abstract_bipartite_graph> &B,
                                                            const shared_ptr<unordered_map<uint32_t, shared_ptr<mutex>>> &left_mutex_map,
                                                            const shared_ptr<unordered_map<uint32_t, shared_ptr<mutex>>> &right_mutex_map,
                                                            const shared_ptr<unordered_map<uint32_t, shared_ptr<bipartite_core_left_store_index>>> &left_index_map,
                                                            const shared_ptr<unordered_map<uint32_t, shared_ptr<bipartite_core_right_store_index>>> &right_index_map,
                                                            const shared_ptr<thread_pool> &pool) {
        auto max_delta = find_core2(B);

        for (auto delta = 1; delta <= max_delta; ++delta) {
            pool->submit_task([=] {
                compute_alpha_core2(B, left_mutex_map, right_mutex_map, left_index_map, right_index_map, delta);
            });
        }

        for (auto delta = 1; delta <= max_delta; ++delta) {
            pool->submit_task([=] {
                compute_beta_core2(B, left_mutex_map, right_mutex_map, left_index_map, right_index_map, delta);
            });
        }
        pool->barrier();

        return max_delta;
    }


    uint32_t share_bipartite_core_decomposition::find_core(const shared_ptr<abstract_bipartite_graph> &B) {
        auto vertex_degree_map = make_shared<unordered_map<uint32_t, uint32_t>>();
        auto l_set = make_shared<unordered_set<uint32_t>>();
        auto r_set = make_shared<unordered_set<uint32_t>>();

        for (const auto &[l, l_vertex]: *B->get_left_vertex_map()) {
            l_set->insert(l);
            vertex_degree_map->insert({l, l_vertex->get_degree()});
        }
        for (const auto &[r, r_vertex]: *B->get_right_vertex_map()) {
            r_set->insert(r);
            vertex_degree_map->insert({r, r_vertex->get_degree()});
        }

        auto evicted_l_set = make_shared<unordered_set<uint32_t>>();
        auto evicted_r_set = make_shared<unordered_set<uint32_t>>();

        uint32_t max_delta = 1;
        for (uint32_t delta = 1; !vertex_degree_map->empty() && !vertex_degree_map->empty(); ++delta) {
            max_delta = delta;
            for (const auto &l: *l_set) {
                if (vertex_degree_map->at(l) <= delta) {
                    evicted_l_set->insert(l);
                }
            }

            for (const auto &r: *r_set) {
                if (vertex_degree_map->at(r) <= delta) {
                    evicted_r_set->insert(r);
                }
            }

            while (!evicted_l_set->empty() || !evicted_r_set->empty()) {
                while (!evicted_l_set->empty()) {
                    auto l = *evicted_l_set->begin();
                    evicted_l_set->erase(l);
                    l_set->erase(l);
                    vertex_degree_map->erase(l);

                    for (const auto &[r, e]: *B->get_left_vertex(l)->get_edge_map()) {
                        if (r_set->count(r)) {
                            --vertex_degree_map->at(r);
                            if (vertex_degree_map->at(r) <= delta) {
                                evicted_r_set->insert(r);
                            }
                        }
                    }
                }

                while (!evicted_r_set->empty()) {
                    auto r = *evicted_r_set->begin();
                    evicted_r_set->erase(r);
                    r_set->erase(r);
                    vertex_degree_map->erase(r);

                    for (const auto &[l, e]: *B->get_right_vertex(r)->get_edge_map()) {
                        if (l_set->count(l)) {
                            --vertex_degree_map->at(l);
                            if (vertex_degree_map->at(l) <= delta) {
                                evicted_l_set->insert(l);
                            }
                        }
                    }
                }
            }
        }
        return max_delta;
    }

    uint32_t share_bipartite_core_decomposition::find_core2(const shared_ptr<abstract_bipartite_graph> &B) {
        auto left_degree_map = make_shared<unordered_map<uint32_t, uint32_t>>();
        auto right_degree_map = make_shared<unordered_map<uint32_t, uint32_t>>();

        for (const auto &[l, l_vertex]: *B->get_left_vertex_map()) {
            left_degree_map->insert({l, l_vertex->get_degree()});
        }
        for (const auto &[r, r_vertex]: *B->get_right_vertex_map()) {
            right_degree_map->insert({r, r_vertex->get_degree()});
        }

        auto evicted_l_set = make_shared<unordered_set<uint32_t>>();
        auto evicted_r_set = make_shared<unordered_set<uint32_t>>();

        uint32_t max_delta = 1;
        for (uint32_t delta = 1; !left_degree_map->empty() && !right_degree_map->empty(); ++delta) {
            max_delta = delta;
            for (const auto &[l, l_degree]: *left_degree_map) {
                if (l_degree <= delta) {
                    evicted_l_set->insert(l);
                }
            }

            for (const auto &[r, r_degree]: *right_degree_map) {
                if (r_degree <= delta) {
                    evicted_r_set->insert(r);
                }
            }

            while (!evicted_l_set->empty() || !evicted_r_set->empty()) {
                while (!evicted_l_set->empty()) {
                    auto l = *evicted_l_set->begin();
                    evicted_l_set->erase(l);
                    left_degree_map->erase(l);

                    for (const auto &[r, e]: *B->get_left_vertex(l)->get_edge_map()) {
                        if (right_degree_map->count(r) && right_degree_map->at(r) > delta) {
                            --right_degree_map->at(r);
                            if (right_degree_map->at(r) == delta) {
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
                        if (left_degree_map->count(l) && left_degree_map->at(l) > delta) {
                            --left_degree_map->at(l);
                            if (left_degree_map->at(l) == delta) {
                                evicted_l_set->insert(l);
                            }
                        }
                    }
                }
            }
        }
        return max_delta;
    }
}
