
#include "truss/mixed_structure_truss_maintenance.h"

namespace scnu {

    uint32_t
    mixed_structure_truss_maintenance::compute_pre_truss_number(const shared_ptr<abstract_graph> &G,
                                                    const shared_ptr<unordered_map<shared_ptr<abstract_edge>, uint32_t>> &edge_truss_map,
                                                    const shared_ptr<abstract_edge>& e) {
        auto result_map = make_shared<map<uint32_t,uint32_t> >();
        auto &e1 = e;
        auto u = e1->get_source_vertex_id();
        auto v = e1->get_destination_vertex_id();

        if(G->get_vertex(u)->get_degree() > G->get_vertex(v)->get_degree()){
            swap(u, v);
        }
        for (const auto &[w,e2]:*G->get_vertex(u)->get_edge_map()) {
            if(w == v){
                continue;
            }
            auto e3 = G->get_edge(v,w);
            if(!e3){
                continue;
            }

            uint32_t min_k = std::min(edge_truss_map->at(e2),edge_truss_map->at(e3));
            if (!result_map->count(min_k)) {
                result_map->insert({min_k, 0});
            }
            ++result_map->at(min_k);
        }

        uint32_t k = 2;
        if(!result_map->empty()){
            uint32_t sum = 0;
            k = result_map->rbegin()->first;
            while (k >= 2) {
                if(result_map->count(k)){
                    sum += result_map->at(k);
                }
                if (sum >= k - 2) {
                    break;
                }
                --k;
            }
        }
        return k;
    }

    uint32_t mixed_structure_truss_maintenance::compute_truss_support(const shared_ptr<abstract_graph> &G,
                                                                      const shared_ptr<unordered_map<shared_ptr<abstract_edge>, uint32_t>> &edge_truss_map,
                                                                      const shared_ptr<::abstract_edge> &e1,
                                                                      uint32_t k) {
        uint32_t support = 0;

        auto u = e1->get_source_vertex_id();
        auto v = e1->get_destination_vertex_id();

        if (G->get_vertex(u)->get_degree() > G->get_vertex(v)->get_degree()) {
            swap(u, v);
        }

        for (const auto &[w, e2]: *G->get_vertex(u)->get_edge_map()) {
            if (w == v || edge_truss_map->at(e2) < k) {
                continue;
            }

            auto e3 = G->get_edge(v, w);
            if (!e3 || edge_truss_map->at(e3) < k) {
                continue;
            }

            ++support;
        }

        return support;
    }


    void mixed_structure_truss_maintenance::decremental_traversal(const shared_ptr<abstract_graph> &G,
                                                                  const shared_ptr<unordered_map<shared_ptr<abstract_edge>, uint32_t>> &edge_truss_map,
                                                                  const shared_ptr<map<uint32_t, shared_ptr<unordered_set<shared_ptr<abstract_edge>>>>> &M) {
        auto removed_map = make_shared<map<uint64_t, shared_ptr<unordered_set<shared_ptr<abstract_edge>>>>>();
        for (const auto &[k, k_set]: *M) {
            removed_map->insert({k, shared_ptr<unordered_set<shared_ptr<abstract_edge>>>()});
        }

        for (auto iter = M->begin(); iter != M->end(); ++iter) {
            auto k = iter->first;
            auto k_set = iter->second;

            removed_map->at(k) = make_shared<unordered_set<shared_ptr<abstract_edge>>>();

            auto S = make_shared<unordered_map<shared_ptr<abstract_edge>, uint32_t>>();
            auto Q = make_shared<unordered_set<shared_ptr<abstract_edge>>>();

            for (const auto &e: *k_set) {
                S->insert({e, compute_truss_support(G, edge_truss_map, e, k)});

                if (S->at(e) < k - 2) {
                    Q->insert(e);
                }
            }


            auto &removed_set = removed_map->at(k);

            while (!Q->empty()) {
                auto e1 = *Q->begin();
                Q->erase(e1);

                removed_set->insert(e1);

                auto u = e1->get_source_vertex_id();
                auto v = e1->get_destination_vertex_id();
                if (G->get_vertex(u)->get_degree() > G->get_vertex(v)->get_degree()) {
                    swap(u, v);
                }
                for (const auto &[w, e2]: *G->get_vertex(u)->get_edge_map()) {
                    if (w == v || edge_truss_map->at(e2) < k || removed_set->count(e2)) {
                        continue;
                    }

                    auto e3 = G->get_edge(v, w);
                    if (!e3 || edge_truss_map->at(e3) < k || removed_set->count(e3)) {
                        continue;
                    }

                    if (edge_truss_map->at(e2) == k) {
                        if (!S->count(e2)) {
                            S->insert({e2, compute_truss_support(G, edge_truss_map, e2, k)});
                        }

                        --S->at(e2);
                        if (S->at(e2) < k - 2) {
                            Q->insert(e2);
                        }
                    }

                    if (edge_truss_map->at(e3) == k) {
                        if (!S->count(e3)) {
                            S->insert({e3, compute_truss_support(G, edge_truss_map, e3, k)});
                        }

                        --S->at(e3);
                        if (S->at(e3) < k - 2) {
                            Q->insert(e3);
                        }
                    }
                }
            }
        }

        for(const auto &[k, k_set]:*removed_map){
            for(const auto &e:*k_set){
                edge_truss_map->at(e) = k - 1;
            }
        }
    }

    void mixed_structure_truss_maintenance::decremental_traversal(const shared_ptr<abstract_graph> &G,
                                                                  const shared_ptr<unordered_map<shared_ptr<abstract_edge>, uint32_t>> &edge_truss_map,
                                                                  const shared_ptr<map<uint32_t, shared_ptr<unordered_set<shared_ptr<abstract_edge>>>>> &M,
                                                                  const shared_ptr<thread_pool>& pool) {
        auto removed_map = make_shared<map<uint64_t, shared_ptr<unordered_set<shared_ptr<abstract_edge>>>>>();
        for(const auto &[k, k_set]:*M){
            removed_map->insert({k, shared_ptr<unordered_set<shared_ptr<abstract_edge>>>()});
        }
        for (const auto &[k, k_set]:*M) {
           pool->submit_task([=] {
               removed_map->at(k) = make_shared<unordered_set<shared_ptr<abstract_edge>>>();

               auto S = make_shared<unordered_map<shared_ptr<abstract_edge>, uint32_t>>();
               auto Q = make_shared<unordered_set<shared_ptr<abstract_edge>>>();

               for (const auto &e: *k_set) {
                   S->insert({e, compute_truss_support(G, edge_truss_map, e, k)});

                   if (S->at(e) < k - 2) {
                       Q->insert(e);
                   }
               }

               auto &removed_set = removed_map->at(k);
               while (!Q->empty()) {
                   auto e1 = *Q->begin();
                   Q->erase(e1);

                   removed_set->insert(e1);

                   auto u = e1->get_source_vertex_id();
                   auto v = e1->get_destination_vertex_id();
                   if (G->get_vertex(u)->get_degree() > G->get_vertex(v)->get_degree()) {
                       swap(u, v);
                   }
                   for (const auto &[w, e2]: *G->get_vertex(u)->get_edge_map()) {
                       if (w == v || edge_truss_map->at(e2) < k || removed_set->count(e2)) {
                           continue;
                       }

                       auto e3 = G->get_edge(v, w);
                       if (!e3 || edge_truss_map->at(e3) < k || removed_set->count(e3)) {
                           continue;
                       }

                       if (edge_truss_map->at(e2) == k) {
                           if (!S->count(e2)) {
                               S->insert({e2, compute_truss_support(G, edge_truss_map, e2, k)});
                           }

                           --S->at(e2);
                           if (S->at(e2) < k - 2) {
                               Q->insert(e2);
                           }
                       }

                       if (edge_truss_map->at(e3) == k) {
                           if (!S->count(e3)) {
                               S->insert({e3, compute_truss_support(G, edge_truss_map, e3, k)});
                           }

                           --S->at(e3);
                           if (S->at(e3) < k - 2) {
                               Q->insert(e3);
                           }
                       }
                   }
               }
            });
        }
        pool->barrier();

        for (const auto &[k, k_set]: *removed_map) {
            pool->submit_task([=] {
                for (const auto &e: *k_set) {
                    edge_truss_map->at(e) = k - 1;
                }
            });
        }
        pool->barrier();
    }


    void mixed_structure_truss_maintenance::incremental_traversal(const shared_ptr<abstract_graph> &G,
                                                                  const shared_ptr<unordered_map<shared_ptr<abstract_edge>, uint32_t>> &edge_truss_map,
                                                                  const shared_ptr<map<uint32_t, shared_ptr<unordered_set<shared_ptr<abstract_edge>>>>> &M) {
        auto inserted_map = make_shared<map<uint64_t, shared_ptr<unordered_set<shared_ptr<abstract_edge>>>>>();
        for(const auto &[k, k_set]:*M){
            inserted_map->insert({k, shared_ptr<unordered_set<shared_ptr<abstract_edge>>>()});
        }
        for (auto iter = M->rbegin(); iter!=M->rend(); ++iter) {
            auto k = iter->first;
            auto k_set = iter->second;

            auto V = make_shared<unordered_set<shared_ptr<abstract_edge>>>();
            auto Q = make_shared<unordered_set<shared_ptr<abstract_edge>>>();
            for (const auto &e: *k_set) {
                Q->insert(e);
                V->insert(e);
            }

            inserted_map->at(k) = make_shared<unordered_set<shared_ptr<abstract_edge>>>();
            auto &candidate_set =  inserted_map->at(k);


            auto S = make_shared<unordered_map<shared_ptr<abstract_edge>, uint32_t>>();
            auto evicted_set = make_shared<unordered_set<shared_ptr<abstract_edge>>>();

            while (!Q->empty()) {
                auto e1 = *Q->begin();
                Q->erase(e1);

                candidate_set->insert(e1);
                S->insert({e1, 0});

                auto e1_set = make_shared<unordered_set<shared_ptr<abstract_edge>>>();

                auto u = e1->get_source_vertex_id();
                auto v = e1->get_destination_vertex_id();

                if (G->get_vertex(u)->get_degree() > G->get_vertex(v)->get_degree()) {
                    swap(u, v);
                }
                for (const auto &[w, e2]: *G->get_vertex(u)->get_edge_map()) {
                    if (w == v || edge_truss_map->at(e2) < k) {
                        continue;
                    }

                    auto e3 = G->get_edge(v, w);
                    if (!e3 || edge_truss_map->at(e3) < k) {
                        continue;
                    }

                    ++S->at(e1);

                    if (edge_truss_map->at(e2) == k && !V->count(e2)) {
                        e1_set->insert(e2);
                    }

                    if (edge_truss_map->at(e3) == k && !V->count(e3)) {
                        e1_set->insert(e3);
                    }
                }

                if (S->at(e1) > k - 2) {
                    for (const auto &e: *e1_set) {
                        Q->insert(e);
                        V->insert(e);
                    }
                } else {
                    evicted_set->insert(e1);
                }
            }

            while (!evicted_set->empty()) {
                auto e1 = *evicted_set->begin();
                evicted_set->erase(e1);

                candidate_set->erase(e1);
                S->erase(e1);

                auto u = e1->get_source_vertex_id();
                auto v = e1->get_destination_vertex_id();
                if (G->get_vertex(u)->get_degree() > G->get_vertex(v)->get_degree()) {
                    swap(u, v);
                }
                for (const auto &[w, e2]: *G->get_vertex(u)->get_edge_map()) {
                    if (w == v || edge_truss_map->at(e2) < k || edge_truss_map->at(e2) == k && !candidate_set->count(e2)) {
                        continue;
                    }

                    auto e3 = G->get_edge(v, w);
                    if (!e3 || edge_truss_map->at(e3) < k || edge_truss_map->at(e3) == k && !candidate_set->count(e3)) {
                        continue;
                    }

                    if (candidate_set->count(e2)) {
                        --S->at(e2);
                        if (S->at(e2) <= k - 2) {
                            evicted_set->insert(e2);
                        }
                    }
                    if (candidate_set->count(e3)) {
                        --S->at(e3);
                        if (S->at(e3) <= k - 2) {
                            evicted_set->insert(e3);
                        }
                    }
                }
            }
        }

        for(const auto &[k, k_set]:*inserted_map){
            for (const auto &e: *k_set) {
                edge_truss_map->at(e) = k + 1;
            }
        }
    }

    void mixed_structure_truss_maintenance::incremental_traversal(const shared_ptr<abstract_graph> &G,
                                                                  const shared_ptr<unordered_map<shared_ptr<abstract_edge>, uint32_t>> &edge_truss_map,
                                                                  const shared_ptr<map<uint32_t, shared_ptr<unordered_set<shared_ptr<abstract_edge>>>>> &M,
                                                                  const shared_ptr<thread_pool>& pool) {
        auto inserted_map = make_shared<map<uint64_t, shared_ptr<unordered_set<shared_ptr<abstract_edge>>>>>();
        for(const auto &[k, k_set]:*M){
            inserted_map->insert({k, shared_ptr<unordered_set<shared_ptr<abstract_edge>>>()});
        }
        for (const auto &[k, k_set]:*M) {
            pool->submit_task([=] {
                auto V = make_shared<unordered_set<shared_ptr<abstract_edge>>>();
                auto Q = make_shared<unordered_set<shared_ptr<abstract_edge>>>();
                for (const auto &e: *k_set) {
                    Q->insert(e);
                    V->insert(e);
                }

                inserted_map->at(k) = make_shared<unordered_set<shared_ptr<abstract_edge>>>();
                auto &candidate_set =  inserted_map->at(k);


                auto S = make_shared<unordered_map<shared_ptr<abstract_edge>, uint32_t>>();
                auto evicted_set = make_shared<unordered_set<shared_ptr<abstract_edge>>>();

                while (!Q->empty()) {
                    auto e1 = *Q->begin();
                    Q->erase(e1);

                    candidate_set->insert(e1);
                    S->insert({e1, 0});

                    auto e1_set = make_shared<unordered_set<shared_ptr<abstract_edge>>>();

                    auto u = e1->get_source_vertex_id();
                    auto v = e1->get_destination_vertex_id();

                    if (G->get_vertex(u)->get_degree() > G->get_vertex(v)->get_degree()) {
                        swap(u, v);
                    }
                    for (const auto &[w, e2]: *G->get_vertex(u)->get_edge_map()) {
                        if (w == v || edge_truss_map->at(e2) < k) {
                            continue;
                        }

                        auto e3 = G->get_edge(v, w);
                        if (!e3 || edge_truss_map->at(e3) < k) {
                            continue;
                        }

                        ++S->at(e1);

                        if (edge_truss_map->at(e2) == k && !V->count(e2)) {
                            e1_set->insert(e2);
                        }

                        if (edge_truss_map->at(e3) == k && !V->count(e3)) {
                            e1_set->insert(e3);
                        }
                    }

                    if (S->at(e1) > k - 2) {
                        for (const auto &e: *e1_set) {
                            Q->insert(e);
                            V->insert(e);
                        }
                    } else {
                        evicted_set->insert(e1);
                    }
                }

                while (!evicted_set->empty()) {
                    auto e1 = *evicted_set->begin();
                    evicted_set->erase(e1);

                    candidate_set->erase(e1);
                    S->erase(e1);

                    auto u = e1->get_source_vertex_id();
                    auto v = e1->get_destination_vertex_id();
                    if (G->get_vertex(u)->get_degree() > G->get_vertex(v)->get_degree()) {
                        swap(u, v);
                    }
                    for (const auto &[w, e2]: *G->get_vertex(u)->get_edge_map()) {
                        if (w == v || edge_truss_map->at(e2) < k || edge_truss_map->at(e2) == k && !candidate_set->count(e2)) {
                            continue;
                        }

                        auto e3 = G->get_edge(v, w);
                        if (!e3 || edge_truss_map->at(e3) < k || edge_truss_map->at(e3) == k && !candidate_set->count(e3)) {
                            continue;
                        }

                        if (candidate_set->count(e2)) {
                            --S->at(e2);
                            if (S->at(e2) <= k - 2) {
                                evicted_set->insert(e2);
                            }
                        }
                        if (candidate_set->count(e3)) {
                            --S->at(e3);
                            if (S->at(e3) <= k - 2) {
                                evicted_set->insert(e3);
                            }
                        }
                    }
                }
            });
        }
        pool->barrier();

        for (const auto &[k, k_set]: *inserted_map) {
            pool->submit_task([=] {
                for (const auto &e: *k_set) {
                    edge_truss_map->at(e) = k + 1;
                }
            });
        }
        pool->barrier();
    }

    void mixed_structure_truss_maintenance::insert(const shared_ptr<abstract_graph> &G,
                                                   const shared_ptr<unordered_set<shared_ptr<abstract_edge>>> &EI,
                                                   const shared_ptr<unordered_map<shared_ptr<abstract_edge>, uint32_t>> &edge_truss_map) {
        while (!EI->empty()) {
            /**
             * @brief suppose the inserted vertex set is empty
             */
            G->insert_edge_collection(EI);

            auto E_triangle_MS = make_shared<unordered_set<shared_ptr<abstract_edge>>>();
            auto E_MS = make_shared<unordered_set<shared_ptr<abstract_edge>>>();

            for (const auto &e1: *EI) {
                auto e1_set = make_shared<unordered_set<shared_ptr<abstract_edge>>>();
                e1_set->insert(e1);

                auto u = e1->get_source_vertex_id();
                auto v = e1->get_destination_vertex_id();

                if (G->get_vertex(u)->get_degree() > G->get_vertex(v)->get_degree()) {
                    swap(u, v);
                }
                for (const auto &[w, e2]: *G->get_vertex(u)->get_edge_map()) {
                    if (w == v) {
                        continue;
                    }

                    auto e3 = G->get_edge(v, w);
                    if (!e3) {
                        continue;
                    }

                    e1_set->insert(e2);
                    e1_set->insert(e3);
                }

                bool flag = true;
                for (const auto &e: *e1_set) {
                    if (E_triangle_MS->count(e)) {
                        flag = false;
                        break;
                    }
                }

                if (flag) {
                    E_MS->insert(e1);
                    E_triangle_MS->merge(*e1_set);
                }
            }

            G->remove_edge_collection(EI);

            for (const auto &e: *E_MS) {
                EI->erase(e);
                G->insert_edge(e);
            }

            auto M = make_shared<map<uint32_t, shared_ptr<unordered_set<shared_ptr<abstract_edge>>>>>();
            for (const auto &e1: *E_MS) {
                edge_truss_map->insert({e1, compute_pre_truss_number(G, edge_truss_map, e1)});

                auto u = e1->get_source_vertex_id();
                auto v = e1->get_destination_vertex_id();
                if (G->get_vertex(u)->get_degree() > G->get_vertex(v)->get_degree()) {
                    swap(u, v);
                }
                for (const auto &[w, e2]: *G->get_vertex(u)->get_edge_map()) {
                    if (w == v) {
                        continue;
                    }

                    auto e3 = G->get_edge(v, w);
                    if (!e3) {
                        continue;
                    }

                    auto k = min(edge_truss_map->at(e2), edge_truss_map->at(e3));
                    if (k <= edge_truss_map->at(e1)) {
                        if (edge_truss_map->at(e2) == k) {
                            if (!M->count(k)) {
                                M->insert({k, make_shared<unordered_set<shared_ptr<abstract_edge>>>()});
                            }
                            M->at(k)->insert(e2);
                        }
                        if (edge_truss_map->at(e3) == k) {
                            if (!M->count(k)) {
                                M->insert({k, make_shared<unordered_set<shared_ptr<abstract_edge>>>()});
                            }
                            M->at(k)->insert(e3);
                        }
                    }
                }
            }
            incremental_traversal(G, edge_truss_map, M);
        }
    }

    uint32_t mixed_structure_truss_maintenance::insert(const shared_ptr<abstract_graph> &G,
                                                       const shared_ptr<unordered_set<shared_ptr<abstract_edge>>> &EI,
                                                       const shared_ptr<unordered_map<shared_ptr<abstract_edge>, uint32_t>> &edge_truss_map,
                                                       const shared_ptr<thread_pool> &pool) {
        uint32_t loop_count = 0;
        while (!EI->empty()) {
            /**
             * @brief suppose the inserted vertex set is empty
             */
            G->insert_edge_collection(EI);
            auto E_MS = make_shared<unordered_set<shared_ptr<abstract_edge>>>();
            auto E_triangle_MS = make_shared<unordered_set<shared_ptr<abstract_edge>>>();

            for (const auto &e1: *EI) {
                auto e1_set = make_shared<unordered_set<shared_ptr<abstract_edge>>>();
                e1_set->insert(e1);

                auto u = e1->get_source_vertex_id();
                auto v = e1->get_destination_vertex_id();

                if (G->get_vertex(u)->get_degree() > G->get_vertex(v)->get_degree()) {
                    swap(u, v);
                }
                for (const auto &[w, e2]: *G->get_vertex(u)->get_edge_map()) {
                    if (w == v) {
                        continue;
                    }

                    auto e3 = G->get_edge(v, w);
                    if (!e3) {
                        continue;
                    }

                    e1_set->insert(e2);
                    e1_set->insert(e3);
                }

                bool flag = true;
                for (const auto &e: *e1_set) {
                    if (E_triangle_MS->count(e)) {
                        flag = false;
                        break;
                    }
                }

                if (flag) {
                    E_MS->insert(e1);
                    E_triangle_MS->merge(*e1_set);
                }
            }

            G->remove_edge_collection(EI);

            for (const auto &e: *E_MS) {
                EI->erase(e);
                G->insert_edge(e);
            }

            auto M = make_shared<map<uint32_t, shared_ptr<unordered_set<shared_ptr<abstract_edge>>>>>();
            for (const auto &e1: *E_MS) {
                edge_truss_map->insert({e1, compute_pre_truss_number(G, edge_truss_map, e1)});

                auto u = e1->get_source_vertex_id();
                auto v = e1->get_destination_vertex_id();
                if (G->get_vertex(u)->get_degree() > G->get_vertex(v)->get_degree()) {
                    swap(u, v);
                }
                for (const auto &[w, e2]: *G->get_vertex(u)->get_edge_map()) {
                    if (w == v) {
                        continue;
                    }

                    auto e3 = G->get_edge(v, w);
                    if (!e3) {
                        continue;
                    }

                    auto k = min(edge_truss_map->at(e2), edge_truss_map->at(e3));
                    if (k <= edge_truss_map->at(e1)) {
                        if (edge_truss_map->at(e2) == k) {
                            if (!M->count(k)) {
                                M->insert({k, make_shared<unordered_set<shared_ptr<abstract_edge>>>()});
                            }
                            M->at(k)->insert(e2);
                        }
                        if (edge_truss_map->at(e3) == k) {
                            if (!M->count(k)) {
                                M->insert({k, make_shared<unordered_set<shared_ptr<abstract_edge>>>()});
                            }
                            M->at(k)->insert(e3);
                        }
                    }
                }
            }

            ++loop_count;
            incremental_traversal(G, edge_truss_map, M, pool);
        }
        return loop_count;
    }

    void mixed_structure_truss_maintenance::remove(const shared_ptr<abstract_graph> &G,
                                                   const shared_ptr<unordered_set<shared_ptr<abstract_edge>>> &ED,
                                                   const shared_ptr<unordered_map<shared_ptr<abstract_edge>, uint32_t>> &edge_truss_map) {
        while (!ED->empty()) {

            /**
             * @brief suppose the inserted vertex set is empty
             */
            auto E_MS = make_shared<unordered_set<shared_ptr<abstract_edge>>>();
            auto E_triangle_MS = make_shared<unordered_set<shared_ptr<abstract_edge>>>();
            for (const auto &e1: *ED) {
                auto e1_set = make_shared<unordered_set<shared_ptr<abstract_edge>>>();

                e1_set->insert(e1);

                auto u = e1->get_source_vertex_id();
                auto v = e1->get_destination_vertex_id();

                if (G->get_vertex(u)->get_degree() > G->get_vertex(v)->get_degree()) {
                    swap(u, v);
                }

                for (const auto &[w, e2]: *G->get_vertex(u)->get_edge_map()) {
                    if (w == v) {
                        continue;
                    }

                    auto e3 = G->get_edge(v, w);
                    if (!e3) {
                        continue;
                    }

                    e1_set->insert(e2);
                    e1_set->insert(e3);
                }

                bool flag = true;
                for (const auto &e: *e1_set) {
                    if (E_triangle_MS->count(e)) {
                        flag = false;
                        break;
                    }
                }

                if (flag) {
                    E_MS->insert(e1);
                    E_triangle_MS->merge(*e1_set);
                }
            }

            auto M = make_shared<map<uint32_t, shared_ptr<unordered_set<shared_ptr<abstract_edge>>>>>();

            for (const auto &e1: *E_MS) {
                auto u = e1->get_source_vertex_id();
                auto v = e1->get_destination_vertex_id();

                if (G->get_vertex(u)->get_degree() > G->get_vertex(v)->get_degree()) {
                    swap(u, v);
                }

                for (const auto &[w, e2]: *G->get_vertex(u)->get_edge_map()) {
                    if (w == v) {
                        continue;
                    }

                    auto e3 = G->get_edge(v, w);
                    if (!e3) {
                        continue;
                    }

                    auto k = min(edge_truss_map->at(e2), edge_truss_map->at(e3));
                    if (k <= edge_truss_map->at(e1)) {
                        if (edge_truss_map->at(e2) == k) {
                            if (!M->count(k)) {
                                M->insert({k, make_shared<unordered_set<shared_ptr<abstract_edge>>>()});
                            }
                            M->at(k)->insert(e2);
                        }
                        if (edge_truss_map->at(e3) == k) {
                            if (!M->count(k)) {
                                M->insert({k, make_shared<unordered_set<shared_ptr<abstract_edge>>>()});
                            }
                            M->at(k)->insert(e3);
                        }
                    }
                }
            }

            for (const auto &e: *E_MS) {
                ED->erase(e);
                G->remove_edge(e);
            }

            decremental_traversal(G, edge_truss_map, M);

            for (const auto &e: *E_MS) {
                edge_truss_map->erase(e);
            }
        }
    }

    uint32_t mixed_structure_truss_maintenance::remove(const shared_ptr<abstract_graph> &G,
                                                       const shared_ptr<unordered_set<shared_ptr<abstract_edge>>> &ED,
                                                       const shared_ptr<unordered_map<shared_ptr<abstract_edge>, uint32_t>> &edge_truss_map,
                                                       const shared_ptr<thread_pool> &pool) {
        uint32_t loop_count = 0;
        while (!ED->empty()) {

            /**
             * @brief suppose the inserted vertex set is empty
             */
            auto E_MS = make_shared<unordered_set<shared_ptr<abstract_edge>>>();
            auto E_triangle_MS = make_shared<unordered_set<shared_ptr<abstract_edge>>>();
            for (const auto &e1: *ED) {
                auto e1_set = make_shared<unordered_set<shared_ptr<abstract_edge>>>();

                e1_set->insert(e1);

                auto u = e1->get_source_vertex_id();
                auto v = e1->get_destination_vertex_id();

                if (G->get_vertex(u)->get_degree() > G->get_vertex(v)->get_degree()) {
                    swap(u, v);
                }

                for (const auto &[w, e2]: *G->get_vertex(u)->get_edge_map()) {
                    if (w == v) {
                        continue;
                    }

                    auto e3 = G->get_edge(v, w);
                    if (!e3) {
                        continue;
                    }

                    e1_set->insert(e2);
                    e1_set->insert(e3);
                }

                bool flag = true;
                for (const auto &e: *e1_set) {
                    if (E_triangle_MS->count(e)) {
                        flag = false;
                        break;
                    }
                }

                if (flag) {
                    E_MS->insert(e1);
                    E_triangle_MS->merge(*e1_set);
                }
            }

            auto M = make_shared<map<uint32_t, shared_ptr<unordered_set<shared_ptr<abstract_edge>>>>>();

            for (const auto &e1: *E_MS) {
                auto u = e1->get_source_vertex_id();
                auto v = e1->get_destination_vertex_id();

                if (G->get_vertex(u)->get_degree() > G->get_vertex(v)->get_degree()) {
                    swap(u, v);
                }

                for (const auto &[w, e2]: *G->get_vertex(u)->get_edge_map()) {
                    if (w == v) {
                        continue;
                    }

                    auto e3 = G->get_edge(v, w);
                    if (!e3) {
                        continue;
                    }

                    auto k = min(edge_truss_map->at(e2), edge_truss_map->at(e3));
                    if (k <= edge_truss_map->at(e1)) {
                        if (edge_truss_map->at(e2) == k) {
                            if (!M->count(k)) {
                                M->insert({k, make_shared<unordered_set<shared_ptr<abstract_edge>>>()});
                            }
                            M->at(k)->insert(e2);
                        }
                        if (edge_truss_map->at(e3) == k) {
                            if (!M->count(k)) {
                                M->insert({k, make_shared<unordered_set<shared_ptr<abstract_edge>>>()});
                            }
                            M->at(k)->insert(e3);
                        }
                    }
                }
            }

            for (const auto &e: *E_MS) {
                ED->erase(e);
                G->remove_edge(e);
            }

            ++loop_count;

            decremental_traversal(G, edge_truss_map, M, pool);

            for (const auto &e: *E_MS) {
                edge_truss_map->erase(e);
            }
        }

        return loop_count;
    }
} 