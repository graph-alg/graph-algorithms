
#include "core/traversal_core_maintenance.h"

namespace scnu {

    uint32_t traversal_core_maintenance::compute_core_rcd(const shared_ptr<abstract_graph> &G,
                                                          uint32_t k,
                                                          uint32_t u,
                                                          const shared_ptr<unordered_map<uint32_t, uint32_t>> &core,
                                                          const shared_ptr<unordered_map<uint32_t, shared_ptr<unordered_map<uint32_t, uint32_t>>>> &rcd,
                                                          uint32_t h) {
        uint32_t value = 0;
        auto u_vertex = G->get_vertex(u);
        for (const auto&[v, e]:*u_vertex->get_edge_map()) {
            if (core->at(v) < k) {
                continue;
            }
            if (h == 1) {
                if (!rcd->count(v)) {
                    rcd->insert({v, make_shared<unordered_map<uint32_t, uint32_t>>()});
                }
                if (!rcd->at(v)->count(h)) {
                    auto v_core_degree = basic_core_decomposition::get_core_degree(G, core, v, k);
                    rcd->at(v)->insert({h, v_core_degree});
                }
                ++value;
            } else {
                if (!rcd->count(v)) {
                    rcd->insert({v, make_shared<unordered_map<uint32_t, uint32_t>>()});
                }
                if (!rcd->at(v)->count(h - 1)) {
                    auto rcdValue = compute_core_rcd(G, k, v, core, rcd, h - 1);
                    rcd->at(v)->insert({h - 1, rcdValue});
                }
                if ((core->at(u) == core->at(v) && rcd->at(v)->at(h - 1) > core->at(u))
                    || core->at(u) < core->at(v)) {
                    ++value;
                }
            }
        }
        return value;
    }

    /**
     * @details traversal insert maintenance algorithm
     * @param G: the given graph
     * @param core: the map of vertex and its core number
     * @param rcd: residential degree
     * @param n: the number of hops (by default, n =2)
     * @param e: inserted edge
     */
    void traversal_core_maintenance::insert(const shared_ptr<abstract_graph> &G,
                                            const std::shared_ptr<unordered_map<uint32_t, uint32_t>> &core,
                                            const std::shared_ptr<unordered_map<uint32_t, shared_ptr<unordered_map<uint32_t, uint32_t>>>> &rcd,
                                            uint32_t n,
                                            const shared_ptr<abstract_edge> &e) {
        /**
         * @brief insert new vertices into core map
         */
        auto u1 = e->get_source_vertex_id();
        auto u2 = e->get_destination_vertex_id();
        if (!core->count(u1)) {
            core->insert({u1, 1});
        }
        if (!core->count(u2)) {
            core->insert({u2, 1});
        }

        /**
         * @brief insert new edge
         */
        auto r = u1;
        if (core->at(u2) < core->at(u1)) {
            r = u2;
        }
        G->insert_edge(e);

        multi_hop_prepare_rcd_insertion(G, core, rcd, n, e);

        /**
         * @brief prepare
         */
        stack<uint32_t> S;
        auto visited = make_shared<unordered_set<uint32_t>>();
        auto evicted = make_shared<unordered_set<uint32_t>>();
        auto changed = make_shared<list<uint32_t>>();
        /**
         * @remarks cd value may less than 0
         */
        auto cd = make_shared<unordered_map<uint32_t, long long>>();
        auto k = core->at(r);

        /**
         * @brief prepare pure core degree (pcd)
         */
        cd->insert({r, rcd->at(r)->at(n)});
        S.push(r);
        visited->insert(r);

        /**
         * @brief DFS traversal
         */
        while (!S.empty()) {
            auto v = S.top();
            S.pop();
            if (cd->at(v) > k) {
                auto v_vertex = G->get_vertex(v);
                for (const auto&[w, e]:*v_vertex->get_edge_map()) {
                    /**
                     * @brief maximal core degree (mcd)
                     */
                    if (core->at(w) == k && rcd->at(w)->at(n - 1) > k && !visited->count(w)) {
                        S.push(w);
                        visited->insert(w);
                        if (!cd->count(w)) {
                            cd->insert({w, 0});
                        }
                        /**
                         * @brief pure core degree (pcd)
                         */
                        cd->at(w) = cd->at(w) + rcd->at(w)->at(n);
                    }
                }
            } else {
                if (!evicted->count(v)) {
                    propagate_eviction(G, core, cd, evicted, k, v);
                }
            }
        }
        for (auto v:*visited) {
            if (!evicted->count(v)) {
                core->at(v) = core->at(v) + 1;
                changed->push_back(v);
            }
        }
        multi_hop_recompute_rcd(G, core, rcd, n, changed, true);
    }

    /**
     * @details evict unsatisfied vertices
     * @param G: the graph
     * @param core: a map of vertex and its core number
     * @param n: the number of hops
     * @param cd: a map of vetex and its maximal core degree
     * @param evicted: a set of unsatisfied vertices
     * @param k: the given k
     * @param v: a vertex
     */
    void traversal_core_maintenance::propagate_eviction(const shared_ptr<abstract_graph> &G,
                                                        const shared_ptr<unordered_map<uint32_t, uint32_t>> &core,
                                                        const shared_ptr<unordered_map<uint32_t, long long>> &cd,
                                                        const shared_ptr<unordered_set<uint32_t>> &evicted,
                                                        uint32_t k,
                                                        uint32_t v) {
        evicted->insert(v);
        auto v_vertex = G->get_vertex(v);
        for (const auto&[w, e]:*v_vertex->get_edge_map()) {
            if (core->at(w) == k) {
                /**
                 * @remarks core degree may be less than 0 due to eviction
                 */
                if (!cd->count(w)) {
                    cd->insert({w, 0});
                }
                cd->at(w) = cd->at(w) - 1;
                if (cd->at(w) == k && !evicted->count(w)) {
                    propagate_eviction(G, core, cd, evicted, k, w);
                }
            }
        }
    }

    /**
     * @details traversal core maintenance under remove case
     * @param G: the given graph
     * @param core: a map of each vertex and its core number
     * @param rcd: a map of each vertex and its residential core degree
     * @param n: the number of hops (by default: n = 1)
     * @param e: remove edge
     */
    void traversal_core_maintenance::remove(const shared_ptr<abstract_graph> &G,
                                            const std::shared_ptr<unordered_map<uint32_t, uint32_t>> &core,
                                            const std::shared_ptr<unordered_map<uint32_t, shared_ptr<unordered_map<uint32_t, uint32_t>>>> &rcd,
                                            uint32_t n,
                                            const shared_ptr<abstract_edge> &e) {
        auto u1 = e->get_source_vertex_id();
        auto u2 = e->get_destination_vertex_id();
        auto r = u1;
        if (core->at(u2) < core->at(u1)) {
            r = u2;
        }
        auto isolated_vertex_set = G->remove_edge(e);

        multi_hop_prepare_rcd_removal(G, core, rcd, n, e);

        auto visited = make_shared<unordered_set<uint32_t>>();
        auto dismissed = make_shared<unordered_set<uint32_t>>();
        /**
         * @brief record the vertices in (k-1)-core
         */
        auto changed = make_shared<list<uint32_t>>();
        auto cd = make_shared<unordered_map<uint32_t, long long>>();
        auto k = core->at(r);
        if (core->at(u1) != core->at(u2)) {
            visited->insert(r);
            /**
             * @brief maximal core degree
             */
            cd->insert({r, rcd->at(r)->at(n-1)});
            if (cd->at(r) < k) {
                propagate_dismissal(G, core, rcd, cd, dismissed, visited, k, r, n-1, changed);
            }
        } else {
            visited->insert(u1);
            /**
             * @brief maximal core degree
             */
            cd->insert({u1, rcd->at(u1)->at(n-1)});
            if (cd->at(u1) < k) {
                propagate_dismissal(G, core, rcd, cd, dismissed, visited, k, u1, n-1, changed);
            }
            visited->insert(u2);
            cd->insert({u2, rcd->at(u2)->at(n-1)});
            if (!dismissed->count(u2) && cd->at(u2) < k) {
                propagate_dismissal(G, core, rcd, cd, dismissed, visited, k, u2, n-1, changed);
            }
        }

        multi_hop_recompute_rcd(G, core, rcd, n, changed, false);
        for(const auto& u:*isolated_vertex_set)
        {
            core->erase(u);
        }
    }

    /**
     * @details recursively delete unsatisfied vertices
     * @param G: the given graph
     * @param core: a map of each vertex and its core number
     * @param rcd: a map of each vertex and its residential core degree
     * @param cd: a map of each vertex and its core degree
     * @param dismissed: a set of  dismissed vertices
     * @param visited: a set of visited vertices
     * @param k: the core number of v
     * @param v: dismissed vertex
     * @param n: the number of hops
     * @param changed: a list of vertices whose core number changed
     */
    void traversal_core_maintenance::propagate_dismissal(const shared_ptr<abstract_graph> &G,
                                                         const shared_ptr<unordered_map<uint32_t, uint32_t>> &core,
                                                         const shared_ptr<unordered_map<uint32_t, shared_ptr<unordered_map<uint32_t, uint32_t>>>> &rcd,
                                                         const shared_ptr<unordered_map<uint32_t, long long>> &cd,
                                                         const shared_ptr<unordered_set<uint32_t>> &dismissed,
                                                         const shared_ptr<unordered_set<uint32_t>> &visited,
                                                         uint32_t k,
                                                         uint32_t v,
                                                         uint32_t n,
                                                         const shared_ptr<list<uint32_t>> &changed) {
        dismissed->insert(v);
        core->at(v) = core->at(v) - 1;
        /**
         * @brief record the vertex that located in (k-1)-core
         */
        changed->push_back(v);
        auto v_vertex = G->get_vertex(v);
        /**
         * @remarks v may only connect one edge, isolated vertices does not spread its influence
         */
        if (!v_vertex) {
            return;
        }
        for (const auto &[w,e]:*v_vertex->get_edge_map()) {
            if (core->at(w) == k) {
                if (!visited->count(w)) {
                    /**
                     * @brief maximal core degree
                     */
                    if (!cd->count(w)) {
                        cd->insert({w, 0});
                    }
                    cd->at(w) = cd->at(w) + rcd->at(w)->at(n);
                    visited->insert(w);
                }
                cd->at(w) = cd->at(w) - 1;
                if (cd->at(w) < k && !dismissed->count(w)) {
                    propagate_dismissal(G, core, rcd, cd, dismissed, visited, k, w, n, changed);
                }
            }
        }
    }

    /**
     * @details prepare rcd values under insert case
     * @param G: the given graph
     * @param core: a map of vertex and its core number
     * @param rcd: the residential core degree of vertices in G
     * @param n: the number of hops
     * @param e: an inserted edge
     */
    void traversal_core_maintenance::multi_hop_prepare_rcd_insertion(const shared_ptr<abstract_graph> &G,
                                                                     const shared_ptr<unordered_map<uint32_t, uint32_t>> &core,
                                                                     const shared_ptr<unordered_map<uint32_t, shared_ptr<unordered_map<uint32_t, uint32_t>>>> &rcd,
                                                                     uint32_t n,
                                                                     const shared_ptr<abstract_edge> &e) {
        auto u1 = e->get_source_vertex_id();
        auto u2 = e->get_destination_vertex_id();

        auto r = u1;
        if (core->at(u2) < core->at(u1)) {
            r = u2;
        }
        auto k = core->at(r);
        auto roots = make_shared<unordered_set<uint32_t>>();
        auto frontiers = make_shared<unordered_map<uint32_t, shared_ptr<unordered_set<uint32_t>>>>();
        for (int h = 1; h <= n + 1; ++h) {
            frontiers->insert({h, make_shared<unordered_set<uint32_t>>()});
        }
        if (core->at(u1) != core->at(u2)) {
            /**
             * @brief insert r into rcd
             */
            if (!rcd->count(r)) {
                rcd->insert({r, make_shared<unordered_map<uint32_t, uint32_t>>()});
            }
            for (int h = 1; h <= n; ++h) {
                /**
                 * insert h into rcd
                 */
                if (!rcd->at(r)->count(h)) {
                    rcd->at(r)->insert({h, 0});
                }
                rcd->at(r)->at(h) = rcd->at(r)->at(h) + 1;
                if (h < n && rcd->at(r)->at(h) == k + 1) {
                    frontiers->at(h + 1)->insert(r);
                }
                if (h > 1) {
                    for (auto v:*frontiers->at(h)) {
                        auto v_vertex = G->get_vertex(v);
                        for (const auto&[w, e]:*v_vertex->get_edge_map()) {
                            if (core->at(w) == k) {
                                if (!rcd->count(w)) {
                                    rcd->insert({w, make_shared<unordered_map<uint32_t, uint32_t>>()});
                                }
                                if (!rcd->at(w)->count(h)) {
                                    rcd->at(w)->insert({h, 0});
                                }
                                rcd->at(w)->at(h) = rcd->at(w)->at(h) + 1;
                                if (h < n && rcd->at(w)->at(h) == k + 1) {
                                    frontiers->at(h + 1)->insert(w);
                                }
                            }
                        }
                    }
                }
            }
        } else {
            for (auto h = 1; h <= n; ++h) {
                if (h == 1) {
                    /**
                     * @brief update u1
                     */
                    if (!rcd->count(u1)) {
                        rcd->insert({u1, make_shared<unordered_map<uint32_t, uint32_t>>()});
                    }
                    if (!rcd->at(u1)->count(h)) {
                        rcd->at(u1)->insert({h, 0});
                    }
                    rcd->at(u1)->at(h) = rcd->at(u1)->at(h) + 1;
                    if (rcd->at(u1)->at(h) == k + 1) {
                        frontiers->at(h + 1)->insert(u1);
                    }
                    /**
                     * update u2
                     */
                    if (!rcd->count(u2)) {
                        rcd->insert({u2, make_shared<unordered_map<uint32_t, uint32_t>>()});
                    }
                    if (!rcd->at(u2)->count(h)) {
                        rcd->at(u2)->insert({h, 0});
                    }
                    rcd->at(u2)->at(h) = rcd->at(u2)->at(h) + 1;
                    if (rcd->at(u2)->at(h) == k + 1) {
                        frontiers->at(h + 1)->insert(u2);
                    }
                } else {
                    if (!rcd->count(u1)) {
                        rcd->insert({u1, make_shared<unordered_map<uint32_t, uint32_t>>()});
                    }
                    if (!rcd->at(u1)->count(h)) {
                        rcd->at(u1)->insert({h, 0});
                    }
                    if (rcd->at(u2)->at(h - 1) > k) {
                        rcd->at(u1)->at(h) = rcd->at(u1)->at(h) + 1;
                        if (h < n && rcd->at(u1)->at(h) == k + 1) {
                            frontiers->at(h + 1)->insert(u1);
                        }
                    }
                    if (!rcd->count(u2)) {
                        rcd->insert({u2, make_shared<unordered_map<uint32_t, uint32_t>>()});
                    }
                    if (!rcd->at(u2)->count(h)) {
                        rcd->at(u2)->insert({h, 0});
                    }
                    if (rcd->at(u1)->at(h - 1) > k) {
                        rcd->at(u2)->at(h) = rcd->at(u2)->at(h) + 1;
                        if (h < n && rcd->at(u2)->at(h) == k + 1) {
                            frontiers->at(h + 1)->insert(u2);
                        }
                    }

                    for (auto v : *frontiers->at(h)) {
                        auto v_vertex = G->get_vertex(v);
                        for (const auto&[w, e]:*v_vertex->get_edge_map()) {
                            if (!(v == u1 && w == u2) && !(v == u2 && w == u1) && core->at(w) == k) {
                                if (!rcd->count(w)) {
                                    rcd->insert({w, make_shared<unordered_map<uint32_t, uint32_t>>()});
                                }
                                if (!rcd->at(w)->count(h)) {
                                    rcd->at(w)->insert({h, 0});
                                }
                                rcd->at(w)->at(h) = rcd->at(w)->at(h) + 1;
                                if (h < n && rcd->at(w)->at(h) == k + 1) {
                                    frontiers->at(h + 1)->insert(w);
                                }
                            }
                        }
                    }
                }
            }
        }
    }

    /**
     * @details compute multi-hop rcd value
     * @param G: the graph
     * @param core: the core number of vertices in G
     * @param rcd: residential core degrees
     * @param n: the number of hops
     * @param changed: a list of vertices with updated core number
     * @param state: indicate insert or remove case
     */
    void traversal_core_maintenance::multi_hop_recompute_rcd(const shared_ptr<abstract_graph> &G,
                                                             const shared_ptr<unordered_map<uint32_t, uint32_t>> &core,
                                                             const shared_ptr<unordered_map<uint32_t, shared_ptr<unordered_map<uint32_t, uint32_t>>>> &rcd,
                                                             uint32_t n,
                                                             const shared_ptr<list<uint32_t>> &changed,
                                                             bool state) {
        unordered_set<uint32_t> visited;
        /**
         * @brief conserve vertex in changed since it will be changed
         */
        //auto changed_backup = make_shared<list<uint32_t>>();

        for (auto v:*changed) {
            visited.insert(v);
            //changed->push_back(v);
        }

        for (auto h = 1; h <= n; ++h) {
            list<uint32_t> updated;
            for (auto v:*changed) {
                auto v_vertex = G->get_vertex(v);
                if (!v_vertex) {
                    continue;
                }
                for (const auto&[w, e]:*v_vertex->get_edge_map()) {
                    /**
                     * @details two cases:  (1) state = true represents the insert case
                    (2) state = false represents the remove case
                     */
                    if (state) {
                        if (!visited.count(w) && (core->at(w) == core->at(v) || core->at(w) == core->at(v) - 1)) {
                            updated.push_back(w);
                            visited.insert(w);
                        }
                    } else {
                        if (!visited.count(w) && (core->at(w) == core->at(v) || core->at(w) == core->at(v) + 1)) {
                            updated.push_back(w);
                            visited.insert(w);
                        }
                    }
                }
            }
            /**
             * @brief copy elements from update to changed
             */
            copy(updated.begin(), updated.end(), back_inserter(*changed));
            for (auto v:*changed) {
                if (!rcd->count(v)) {
                    rcd->insert({v, make_shared<unordered_map<uint32_t, uint32_t>>()});
                }
                if (!rcd->at(v)->count(h)) {
                    rcd->at(v)->insert({h, 0});
                }
                auto rcd_value = compute_rcd(G, v, core, rcd, h);
                rcd->at(v)->at(h) = rcd_value;
            }
        }
    }

    /**
     * @details prepare residential core degree (rcd) of each vertex in G with multiple-hops
     * @param G: the given graph
     * @param core: a map of a vertex and its core number
     * @param rcd: a map of residential core degree (rcd) for each vertex with multi-hops
     * @param n: the number of hops
     */
    void traversal_core_maintenance::init_rcd(const shared_ptr<abstract_graph> &G,
                                              const shared_ptr<unordered_map<uint32_t, uint32_t>> &core,
                                              const shared_ptr<unordered_map<uint32_t, shared_ptr<unordered_map<uint32_t, uint32_t>>>> &rcd,
                                              uint32_t n) {
        for (auto h = 1; h <= n; ++h) {
            for (const auto &[v,e]:*G->get_vertex_map()) {
                if (!rcd->count(v)) {
                    rcd->insert({v, make_shared<unordered_map<uint32_t, uint32_t>>()});
                }
                auto rcd_value = compute_rcd(G, v, core, rcd, h);
                rcd->at(v)->insert({h, rcd_value});
            }
        }
    }


    int traversal_core_maintenance::compute_rcd(const shared_ptr<abstract_graph> &G,
                                                uint32_t u,
                                                const shared_ptr<unordered_map<uint32_t, uint32_t>> &core,
                                                const shared_ptr<unordered_map<uint32_t, shared_ptr<unordered_map<uint32_t, uint32_t>>>> &rcd,
                                                uint32_t h) {
        auto u_vertex = G->get_vertex(u);
        int value = 0;
        if (!u_vertex) {
            return 0;
        }
        for (const auto &v_pair:*u_vertex->get_edge_map()) {
            auto v = v_pair.first;
            if (h == 1) {
                if (core->at(u) <= core->at(v)) {
                    ++value;
                }
            } else {
                if ((core->at(u) == core->at(v) && rcd->at(v)->at(h - 1) > core->at(u)) || core->at(u) < core->at(v)) {
                    ++value;
                }
            }
        }
        return value;
    }


    /**
     * @details compute residential core degree for the remove case
     * @param G: the given graph
     * @param core: a map of each vertex and its core number
     * @param rcd: a map of each vertex and its residential core degree
     * @param n: the number of hops
     * @param e: removed edge
     */
    void traversal_core_maintenance::multi_hop_prepare_rcd_removal(const shared_ptr<abstract_graph> &G,
                                                                   const shared_ptr<unordered_map<uint32_t, uint32_t>> &core,
                                                                   const shared_ptr<unordered_map<uint32_t, shared_ptr<unordered_map<uint32_t, uint32_t>>>> &rcd,
                                                                   uint32_t n,
                                                                   const shared_ptr<abstract_edge> &e) {
        auto u1 = e->get_source_vertex_id();
        auto u2 = e->get_destination_vertex_id();
        auto r = u1;
        if (core->at(u2) < core->at(u1)) {
            r = u2;
        }
        auto k = core->at(r);
        auto frontiers = make_shared<unordered_map<uint32_t, shared_ptr<unordered_set<uint32_t>>>>();
        for (int h = 1; h <= n + 1; ++h) {
            frontiers->insert({h, make_shared<unordered_set<uint32_t>>()});
        }
        if (core->at(u1) != core->at(u2)) {
            for (int h = 1; h <= n; ++h) {
                rcd->at(r)->at(h) = rcd->at(r)->at(h) - 1;
                if (h < n && rcd->at(r)->at(h) == k) {
                    frontiers->at(h + 1)->insert(r);
                }
                if (h > 1) {
                    for (auto v: *frontiers->at(h)) {
                        auto v_vertex = G->get_vertex(v);
                        for (const auto &[w,e]:*v_vertex->get_edge_map()) {
                            if (core->at(w) == k) {
                                rcd->at(w)->at(h) = rcd->at(w)->at(h) - 1;
                                if (h < n && rcd->at(w)->at(h) == k) {
                                    frontiers->at(h + 1)->insert(w);
                                }
                            }
                        }
                    }
                }
            }
        } else {
            auto old_rcd = make_shared<unordered_map<uint32_t, shared_ptr<unordered_map<uint32_t, uint32_t>>>>();
            old_rcd->insert({u1, make_shared<unordered_map<uint32_t, uint32_t>>()});
            for (auto p:*rcd->at(u1)) {
                old_rcd->at(u1)->insert({p.first, p.second});
            }
            old_rcd->insert({u2, make_shared<unordered_map<uint32_t, uint32_t>>()});
            for (auto p:*rcd->at(u2)) {
                old_rcd->at(u2)->insert({p.first, p.second});
            }
            for (int h = 1; h <= n; ++h) {
                if (h == 1) {
                    rcd->at(u1)->at(h) = rcd->at(u1)->at(h) - 1;
                    if (rcd->at(u1)->at(h) == k) {
                        frontiers->at(h + 1)->insert(u1);
                    }
                    rcd->at(u2)->at(h) = rcd->at(u2)->at(h) - 1;
                    if (rcd->at(u2)->at(h) == k) {
                        frontiers->at(h + 1)->insert(u2);
                    }
                } else {
                    if (old_rcd->at(u2)->at(h - 1) > k) {
                        rcd->at(u1)->at(h) = rcd->at(u1)->at(h) - 1;
                        if (h < n && rcd->at(u1)->at(h) == k) {
                            frontiers->at(h + 1)->insert(u1);
                        }
                    }
                    if (old_rcd->at(u1)->at(h - 1) > k) {
                        rcd->at(u2)->at(h) = rcd->at(u2)->at(h) - 1;
                        if (h < n && rcd->at(u2)->at(h) == k) {
                            frontiers->at(h + 1)->insert(u2);
                        }
                    }
                    for (auto v:*frontiers->at(h)) {
                        auto v_vertex = G->get_vertex(v);
                        for (const auto &[w,e]:*v_vertex->get_edge_map()) {
                            if (!(v == u1 && w == u2) && !(v == u2 && w == u1) && core->at(w) == k) {
                                rcd->at(w)->at(h) = rcd->at(w)->at(h) - 1;
                                if (h < n && rcd->at(w)->at(h) == k) {
                                    frontiers->at(h + 1)->insert(w);
                                }
                            }
                        }
                    }
                }
            }
        }
    }
}