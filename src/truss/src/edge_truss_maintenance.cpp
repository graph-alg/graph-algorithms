
#include "truss/edge_truss_maintenance.h"

namespace scnu{

    shared_ptr<unordered_set<uint32_t>> edge_truss_maintenance::get_common_neighbor_set(const shared_ptr<abstract_graph>&G,
                                                                                        uint32_t u,
                                                                                        uint32_t v) {
        auto result_set = make_shared<unordered_set<uint32_t>>();
        if(G->get_vertex(u)->get_degree() > G->get_vertex(v)->get_degree()){
            swap(u,v);
        }
        auto u_vertex = G->get_vertex(u);
        auto v_vertex = G->get_vertex(v);
        for (const auto &[w,e]: *u_vertex->get_edge_map()) {
            if (v_vertex->get_edge(w)) {
                result_set->insert(w);
            }
        }
        return result_set;
    }


    void edge_truss_maintenance::insert(const shared_ptr<abstract_graph>& G,
                                        const shared_ptr<abstract_edge> &e0,
                                        const shared_ptr<unordered_map<shared_ptr<abstract_edge>, uint32_t>>& edge_truss_map) {

        G->insert_edge(e0);
        auto update_edge_truss_map = make_shared<unordered_map<shared_ptr<abstract_edge>, uint32_t>>();
        auto p = get_k1_k2(G, e0, edge_truss_map);
        auto k1 = p.first;
        auto k2 = p.second;
        auto k_max = k2 - 1;
        edge_truss_map->insert({e0, k1});
        auto l_map = make_shared<unordered_map<uint32_t, shared_ptr<unordered_set<shared_ptr<abstract_edge>>>>>();
        for (auto k = 2; k <= k_max; ++k) {
            l_map->insert({k, make_shared<unordered_set<shared_ptr<abstract_edge>>>() });
        }
        auto u = e0->get_source_vertex_id();
        auto v = e0->get_destination_vertex_id();
        auto W = get_common_neighbor_set(G, u, v);
        for (auto w:*W) {
            auto e1 = G->get_edge(u, w);
            auto e2 = G->get_edge(v, w);

            auto k = std::min(edge_truss_map->at(e1), edge_truss_map->at(e2));
            if (k <= k_max) {
                if (edge_truss_map->at(e1) == k) {
                    l_map->at(k)->insert(e1);
                }
                if (edge_truss_map->at(e2) == k) {
                    l_map->at(k)->insert(e2);
                }
            }
        }
        for (auto k = k_max; k >= 2; k--) {
            auto evicted_set = l_map->at(k);
            auto s_map = make_shared<unordered_map<shared_ptr<abstract_edge>, uint32_t>>();
            while (!evicted_set->empty()) {
                auto e1 = *evicted_set->begin();
                evicted_set->erase(e1);

                s_map->insert({e1, 0});
                auto x = e1->get_source_vertex_id();
                auto y = e1->get_destination_vertex_id();
                auto Z = get_common_neighbor_set(G, x, y);
                for (auto z:*Z) {
                    auto e2 = G->get_edge(x, z);
                    auto e3 = G->get_edge(y, z);
                    if (edge_truss_map->at(e2) < k || edge_truss_map->at(e3) < k) {
                        continue;
                    }
                    ++s_map->at(e1);
                    if (edge_truss_map->at(e2) == k && !l_map->at(k)->count(e2)) {
                        evicted_set->insert(e2);
                        l_map->at(k)->insert(e2);
                    }

                    if (edge_truss_map->at(e3) == k && !l_map->at(k)->count(e3)) {
                        evicted_set->insert(e3);
                        l_map->at(k)->insert(e3);
                    }
                }
            }

            auto s_set = make_shared<unordered_set<shared_ptr<abstract_edge>>>();
            for (const auto& [e, e_support]:*s_map) {
                if (e_support <= k - 2) {
                    s_set->insert(e);
                }
            }
            while (!s_set->empty()) {
                auto e = *s_set->begin();
                s_set->erase(e);

                auto x = e->get_source_vertex_id();
                auto y = e->get_destination_vertex_id();
                l_map->at(k)->erase(e);
                auto Z = get_common_neighbor_set(G, x, y);
                for (auto z:*Z) {
                    auto xz = G->get_edge(x, z);
                    auto yz = G->get_edge(y, z);
                    if (edge_truss_map->at(xz) < k || edge_truss_map->at(yz) < k) {
                        continue;
                    }
                    if (edge_truss_map->at(xz) == k && !l_map->at(k)->count(xz)) {
                        continue;
                    }
                    if (edge_truss_map->at(yz) == k && !l_map->at(k)->count(yz)) {
                        continue;
                    }


                    if (l_map->at(k)->count(xz)) {
                        --s_map->at(xz);
                        if (s_map->at(xz) <= k - 2) {
                            s_set->insert(xz);
                        }
                    }
                    if (l_map->at(k)->count(yz)) {
                        --s_map->at(yz);
                        if (s_map->at(yz) <= k - 2) {
                            s_set->insert(yz);
                        }
                    }
                }
            }
            for (const auto& e:*l_map->at(k)) {
                update_edge_truss_map->insert({e, k + 1});
            }
        }
        for (const auto &[e,truss_number]:*update_edge_truss_map) {
            edge_truss_map->at(e) = truss_number;
        }
    }

    void edge_truss_maintenance::remove(const shared_ptr<abstract_graph>& G,
                                        const shared_ptr<abstract_edge>& e0,
                                        const shared_ptr<unordered_map<shared_ptr<abstract_edge>, uint32_t>>& edge_truss_map) {
        auto k_max = edge_truss_map->at(e0);
        auto l_map = make_shared<unordered_map<uint32_t, shared_ptr<unordered_set<shared_ptr<abstract_edge>>>>>();
        for (auto k = 2; k <= k_max; ++k) {
            l_map->insert({k, make_shared<unordered_set<shared_ptr<abstract_edge>>>()});
        }
        auto u = e0->get_source_vertex_id();
        auto v = e0->get_destination_vertex_id();
        auto W = get_common_neighbor_set(G, u, v);
        for (auto w:*W) {
            auto e1 = G->get_edge(u, w);
            auto e2 = G->get_edge(v, w);

            auto k = std::min(edge_truss_map->at(e1), edge_truss_map->at(e2));
            if (k <= k_max) {
                if (edge_truss_map->at(e1) == k) {
                    l_map->at(k)->insert(e1);
                }
                if (edge_truss_map->at(e2) == k) {
                    l_map->at(k)->insert(e2);
                }
            }
        }
        G->remove_edge(e0);
        edge_truss_map->erase(e0);

        for (auto k = k_max; k >= 2; k--) {
            auto evicted_set = l_map->at(k);
            auto s_map = make_shared<unordered_map<shared_ptr<abstract_edge>, uint32_t>>();

            while (!evicted_set->empty()) {
                auto e1 = *evicted_set->begin();
                evicted_set->erase(e1);
                s_map->insert({e1, 0});
                auto u1 = e1->get_source_vertex_id();
                auto v1 = e1->get_destination_vertex_id();
                W = get_common_neighbor_set(G, u1, v1);
                for (auto w1:*W) {
                    auto e2 = G->get_edge(u1, w1);
                    auto e3 = G->get_edge(v1, w1);
                    if (edge_truss_map->at(e2) < k || edge_truss_map->at(e3) < k) {
                        continue;
                    }
                    ++s_map->at(e1);
                    if (edge_truss_map->at(e2) == k && !l_map->at(k)->count(e2)) {
                        evicted_set->insert(e2);
                        l_map->at(k)->insert(e2);
                    }

                    if (edge_truss_map->at(e3) == k && !l_map->at(k)->count(e3)) {
                        evicted_set->insert(e3);
                        l_map->at(k)->insert(e3);
                    }
                }
            }
            auto s_queue = make_shared<unordered_set<shared_ptr<abstract_edge>>>();
            for (const auto&[e, e_support]:*s_map) {
                if (e_support >= k - 2) {
                    l_map->at(k)->erase(e);
                } else {
                    s_queue->insert(e);
                }
            }
            while (!s_queue->empty()) {
                auto e = *s_queue->begin();
                s_queue->erase(e);
                edge_truss_map->at(e) = k - 1;

                auto x = e->get_source_vertex_id();
                auto y = e->get_destination_vertex_id();
                auto Z = get_common_neighbor_set(G, x, y);
                for (auto z:*Z) {
                    auto xz = G->get_edge(x, z);
                    auto yz = G->get_edge(y, z);
                    if (edge_truss_map->at(xz) < k || edge_truss_map->at(yz) < k) {
                        continue;
                    }

                    if (edge_truss_map->at(xz) == k) {
                        --s_map->at(xz);
                        if (s_map->at(xz) < k - 2 && !l_map->at(k)->count(xz)) {
                            s_queue->insert(xz);
                            l_map->at(k)->insert(xz);
                        }
                    }
                    if (edge_truss_map->at(yz) == k) {
                        --s_map->at(yz);
                        if (s_map->at(yz) < k - 2 && !l_map->at(k)->count(yz)) {
                            s_queue->insert(yz);
                            l_map->at(k)->insert(yz);
                        }
                    }
                }
            }
        }
    }


    pair<uint32_t, uint32_t> edge_truss_maintenance::get_k1_k2(const shared_ptr<abstract_graph>& G,
            const shared_ptr<abstract_edge>&e1,
            const shared_ptr<unordered_map<shared_ptr<abstract_edge>, uint32_t>>& edge_truss_map){
        auto u = e1->get_source_vertex_id();
        auto v = e1->get_destination_vertex_id();
        auto k1 = 2;
        auto k2 = 2;
        auto W = get_common_neighbor_set(G, u, v);
        auto truss_w_map = make_shared<map<uint32_t, uint32_t>>();
        for (auto w:*W) {
            auto e2 = G->get_edge(u, w);
            auto e3 = G->get_edge(v, w);
            uint32_t min_truss_number = std::min(edge_truss_map->at(e2), edge_truss_map->at(e3));
            for (auto k = 2; k <= min_truss_number; ++k) {
                if (!truss_w_map->count(k)) {
                    truss_w_map->insert({k, 0});
                }
                ++truss_w_map->at(k);
            }
        }
        int k = 3;
        while (truss_w_map->count(k - 1)) {
            if (truss_w_map->count(k) && truss_w_map->at(k) >= k - 2 && k > k1) {
                k1 = k;
            }
            if (truss_w_map->count(k - 1) && truss_w_map->at(k - 1) >= k - 2 && k > k2) {
                k2 = k;
            }
            ++k;
        }
        return std::make_pair(k1, k2);
    }
}