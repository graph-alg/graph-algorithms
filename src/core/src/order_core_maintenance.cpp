
#include "core/order_core_maintenance.h"

namespace scnu {

    /**
     * @details compute residential core degree
     * @param G
     * @param u
     * @param core
     * @param rcd
     * @param h
     * @return
     */
    uint32_t order_core_maintenance::compute_rcd(const shared_ptr<abstract_graph> &G,
                                                 uint32_t u,
                                                 const shared_ptr<unordered_map<uint32_t, uint32_t>> &core,
                                                 const shared_ptr<unordered_map<uint32_t, shared_ptr<unordered_map<uint32_t, uint32_t>>>> &rcd,
                                                 uint32_t h) {
        auto u_vertex = G->get_vertex(u);
        uint32_t value = 0;
        if (!u_vertex) {
            return 0;
        }
        for (const auto &[v,e]:*u_vertex->get_edge_map()) {
            if (h == 1) {
                if (core->at(v) >= core->at(u)) {
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
     * @details initialize the remaining degree for each vertex
     * @param G
     * @param core
     * @param tree
     * @param remaining_degree_map
     */
    void order_core_maintenance::initialize(const shared_ptr<abstract_graph> &G,
                                            const shared_ptr<unordered_map<uint32_t, uint32_t>> &core,
                                            const shared_ptr<unordered_map<uint32_t, shared_ptr<extend_list<double, uint32_t>>>> &tree,
                                            const shared_ptr<unordered_map<uint32_t, uint32_t>> &remaining_degree_map,
                                            const shared_ptr<unordered_map<uint32_t, uint32_t>> &candidate_degree_map,
                                            const shared_ptr<unordered_map<uint32_t,shared_ptr<unordered_map<uint32_t,uint32_t>>>>& rcd,
                                            uint32_t n) {

        /**
         * @details initialize remaining degree and candidate degree
         */
        for (const auto&[u,core_number]:*core) {
            auto u_vertex = G->get_vertex(u);
            candidate_degree_map->insert({u,0});
            auto remaining_degree = 0;
            for (const auto &[v,e]:*u_vertex->get_edge_map()) {
                if (test_order(core, tree, u, v)) {
                    ++remaining_degree;
                }
            }
            remaining_degree_map->insert({u, remaining_degree});
        }

        /**
         * @brief initialize rcd
         */
        for (auto h = 1; h <= n; ++h) {
            for (const auto &[u,u_vertex]:*G->get_vertex_map()) {
                if (!rcd->count(u)) {
                    rcd->insert({u, make_shared<unordered_map<uint32_t, uint32_t>>()});
                }
                auto rcd_value = compute_rcd(G, u, core, rcd, h);
                rcd->at(u)->insert({h, rcd_value});
            }
        }
    }

    /**
    * @details initialize the remaining degree for each vertex
    * @param G
    * @param core
    * @param tree
    * @param remaining_degree_map
    */
    void order_core_maintenance::initialize(const shared_ptr<abstract_graph> &G,
                                            const shared_ptr<unordered_map<uint32_t, uint32_t>> &core,
                                            const shared_ptr<gadget::Treap> &tree,
                                            const shared_ptr<vector<long long>>& root,
                                            const shared_ptr<unordered_map<uint32_t, uint32_t>> &remaining_degree_map,
                                            const shared_ptr<unordered_map<uint32_t, uint32_t>> &candidate_degree_map,
                                            const shared_ptr<unordered_map<uint32_t,shared_ptr<unordered_map<uint32_t,uint32_t>>>>& rcd,
                                            uint32_t n) {

        /**
         * @details initialize remaining degree and candidate degree
         */
        for (const auto&[u,core_number]:*core) {
            auto u_vertex = G->get_vertex(u);
            candidate_degree_map->insert({u,0});
            auto remaining_degree = 0;
            for (const auto &[v,e]:*u_vertex->get_edge_map()) {
                if (test_order(core, tree, u, v)) {
                    ++remaining_degree;
                }
            }
            remaining_degree_map->insert({u, remaining_degree});
        }

        /**
         * @brief initialize rcd
         */
        for (auto h = 1; h <= n; ++h) {
            for (const auto &[u,u_vertex]:*G->get_vertex_map()) {
                if (!rcd->count(u)) {
                    rcd->insert({u, make_shared<unordered_map<uint32_t, uint32_t>>()});
                }
                auto rcd_value = compute_rcd(G, u, core, rcd, h);
                rcd->at(u)->insert({h, rcd_value});
            }
        }
    }

    void order_core_maintenance::multi_hop_prepare_rcd_insertion(const shared_ptr<abstract_graph> &G,
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
      * @details compute residential core degree for the remove case
      * @param G: the given graph
      * @param core: a map of each vertex and its core number
      * @param rcd: a map of each vertex and its residential core degree
      * @param n: the number of hops
      * @param e: removed edge
      */
    void order_core_maintenance::multi_hop_prepare_rcd_removal(const shared_ptr<abstract_graph> &G,
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

    /**
     * @details compute multi-hop rcd value
     * @param G: the graph
     * @param core: the core number of vertices in G
     * @param rcd: residential core degrees
     * @param n: the number of hops
     * @param changed: a list of vertices with updated core number
     * @param state: indicate insert or remove case
     */
    void order_core_maintenance::multi_hop_recompute_rcd(const shared_ptr<abstract_graph> &G,
                                                         const shared_ptr<unordered_map<uint32_t, uint32_t>> &core,
                                                         const shared_ptr<unordered_map<uint32_t, shared_ptr<unordered_map<uint32_t, uint32_t>>>> &rcd,
                                                         uint32_t n,
                                                         const shared_ptr<list<uint32_t>> &changed,
                                                         bool state) {
        unordered_set<uint32_t> visited;
        /**
         * @brief conserve vertex in changed since it will be changed
         */
        auto changed_backup = make_shared<list<uint32_t>>();

        for (auto v:*changed) {
            visited.insert(v);
            changed_backup->push_back(v);
        }

        for (auto h = 1; h <= n; ++h) {
            list<uint32_t> updated;
            for (auto v:*changed_backup) {
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
            copy(updated.begin(), updated.end(), back_inserter(*changed_backup));
            for (auto v:*changed_backup) {
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
     * @details maintenance core order list under the insert case
     * @param G
     * @param e
     * @param core
     * @param k_order
     * @param node_map
     * @param tree
     * @param remaining_degree_map
     * @param candidate_degree_map
     * @param rcd
     * @param n
     */
    void order_core_maintenance::insert(const shared_ptr<abstract_graph> &G,
                                        const shared_ptr<abstract_edge> &e,
                                        const shared_ptr<unordered_map<uint32_t, uint32_t>> &core,
                                        const shared_ptr<unordered_map<uint32_t, shared_ptr<double_list<uint32_t>>>> &k_order,
                                        const shared_ptr<unordered_map<uint32_t, shared_ptr<double_node<uint32_t>>>> &node_map,
                                        const shared_ptr<unordered_map<uint32_t, shared_ptr<extend_list<double, uint32_t>>>> &tree,
                                        const shared_ptr<unordered_map<uint32_t, uint32_t>> &remaining_degree_map,
                                        const shared_ptr<unordered_map<uint32_t, uint32_t>> &candidate_degree_map,
                                        const shared_ptr<unordered_map<uint32_t, shared_ptr<unordered_map<uint32_t, uint32_t> >>> &rcd,
                                        uint32_t n
    ) {
        auto u = e->get_source_vertex_id();
        auto v = e->get_destination_vertex_id();

        /**
         * @brief prepare skip map
         */
        auto B = make_shared<map<double, shared_ptr<double_node<uint32_t>>>>();

        /**
         * @brief deal with vertices are not contained in the graph
         */
        if (!k_order->count(1)) {
            k_order->insert({1, make_shared<double_list<uint32_t>>()});
            tree->insert({1, make_shared<extend_list<double,uint32_t>>()});
        }
        if (!core->count(u)) {
            core->insert({u, 1});
            node_map->insert({u, make_shared<double_node<uint32_t>>(u)});

            k_order->at(1)->left_insert(node_map->at(u));
            tree->at(1)->left_insert(u);

            remaining_degree_map->insert({u, 0});
            candidate_degree_map->insert({u, 0});

            rcd->insert({u, make_shared<unordered_map<uint32_t, uint32_t>>()});
        }
        if (!core->count(v)) {
            core->insert({v, 1});
            node_map->insert({v, make_shared<double_node<uint32_t>>(v)});

            k_order->at(1)->left_insert(node_map->at(v));
            tree->at(1)->left_insert(v);

            remaining_degree_map->insert({v, 0});
            candidate_degree_map->insert({v, 0});

            rcd->insert({v, make_shared<unordered_map<uint32_t, uint32_t>>()});
        }

        /**
         * @brief update rcd value for u and v
         */
        if (core->at(u) <= core->at(v)) {
            if (!rcd->count(u)) {
                rcd->insert({u, make_shared<unordered_map<uint32_t, uint32_t>>()});
                rcd->at(u)->insert({1, 1});
            }
        }
        if (core->at(v) <= core->at(u)) {
            if (!rcd->count(v)) {
                rcd->insert({v, make_shared<unordered_map<uint32_t, uint32_t>>()});
                rcd->at(v)->insert({1, 1});
            }
        }

        auto K = min(core->at(u), core->at(v));
        G->insert_edge(e);

        if (!test_order(core, tree, u, v)) {
            swap(u, v);
        }

        ++remaining_degree_map->at(u);
        if (remaining_degree_map->at(u) <= K) {
            return;
        }

        B->insert({tree->at(K)->find_key(u).value(), node_map->at(u)});

        auto swap = make_shared<vector<uint32_t>>();// conserve the unaffected vertices
        auto k_list = k_order->at(K);
        auto new_k_list = make_shared<double_list<uint32_t>>();
        auto Vc = make_shared<list<uint32_t>>();
        auto Vc_set = make_shared<unordered_set<uint32_t>>();
        auto p = k_list->get_head();
        while (p) {
            auto vi = p->get_value();

            B->erase(tree->at(K)->find_key(vi).value());
            if (candidate_degree_map->at(vi) + remaining_degree_map->at(vi) > K) {
                /**
                 * @brief conserve the next of p
                 */
                auto q = p->get_next();
                Vc_set->insert(vi);
                /**
                 * @remarks note the vertex order
                 */
                Vc->push_front(vi);

                auto vi_vertex = G->get_vertex(vi);
                for (const auto &[w,e]:*vi_vertex->get_edge_map()) {
                    if (core->at(w) == K && test_order(core, tree, vi, w)) {
                        ++candidate_degree_map->at(w);
                        /**
                         * @brief avoid repeatedly insert
                         */
                        if (!B->count(tree->at(K)->find_key(w).value())) {
                            B->insert({tree->at(K)->find_key(w).value(), node_map->at(w)});
                        }
                    }
                }
                /**
                 * @brief move p to next node
                 */
                p = q;
            } else if (candidate_degree_map->at(vi) == 0) {
                shared_ptr<double_node<uint32_t>> q;
                shared_ptr<double_node<uint32_t>> prior_q;
                /**
                 * @details two cases: (1) get the prior node of the first item of B when B is nonempty (2) skip to the rear of this list
                 */
                if (!B->empty()) {
                    q = B->begin()->second;
                    prior_q = q->get_prior();
                } else {
                    prior_q = k_list->get_rear();
                }
                /**
                 * @brief move nodes from old list to the new list
                 */
                new_k_list->right_insert(p, prior_q);
                p = q;
            } else {
                auto q = p->get_next();
                new_k_list->right_insert(p);
                remaining_degree_map->at(vi) = remaining_degree_map->at(vi) + candidate_degree_map->at(vi);
                candidate_degree_map->at(vi) = 0;
                remove_candidates(G, Vc, Vc_set, node_map, tree, B, new_k_list, swap, vi, K, core,
                                  remaining_degree_map,
                                  candidate_degree_map);
                auto pivot_node = tree->at(K)->find(vi);
                p = q;
            }
        }

        k_order->at(K) = new_k_list;

        for (auto w:*swap) {
            auto node = node_map->at(w);

            auto node_prior = node->get_prior();
            /**
             * @brief readjust the order of nodes who do not change their core number
             */
            tree->at(K)->remove(w);
            auto w_extend_node = make_shared<extend_node<double,uint32_t>>(0, w);
            auto prior_extend_node = tree->at(K)->find(node_prior->get_value());

            tree->at(K)->insert_after(w_extend_node, prior_extend_node);
        }

        for (auto w:*Vc) {
            if (!Vc_set->count(w)) {
                continue;
            }

            if (!k_order->count(K + 1)) {
                k_order->insert({K + 1, make_shared<double_list<uint32_t>>()});
                tree->insert({K + 1, make_shared<extend_list<double,uint32_t>>()});
            }
            k_order->at(K + 1)->left_insert(node_map->at(w));

            candidate_degree_map->at(w) = 0;
            core->at(w) = K + 1;

            tree->at(K)->remove(w);
            auto w_extend_node = make_shared<extend_node<double,uint32_t>>(0, w);
            tree->at(K + 1)->left_insert(w_extend_node);
            /**
             * @brief update maximal core degree
             */
            auto rcd_value = compute_rcd(G, w, core, rcd, n);
            rcd->at(w)->insert({n, rcd_value});
        }
    }

    /**
     * @details maintenance core order list under the insert case
     * @param G
     * @param e
     * @param core
     * @param k_order
     * @param node_map
     * @param tree
     * @param remaining_degree_map
     * @param candidate_degree_map
     * @param rcd
     * @param n
     */
    void order_core_maintenance::insert(const shared_ptr<abstract_graph> &G,
                                        const shared_ptr<abstract_edge> &e,
                                        const shared_ptr<unordered_map<uint32_t, uint32_t>> &core,
                                        const shared_ptr<unordered_map<uint32_t, shared_ptr<double_list<uint32_t>>>> &k_order,
                                        const shared_ptr<unordered_map<uint32_t, shared_ptr<double_node<uint32_t>>>> &node_map,
                                        const shared_ptr<gadget::Treap> &tree,
                                        const shared_ptr<vector<long long>> &root,
                                        const shared_ptr<unordered_map<uint32_t, uint32_t>> &remaining_degree_map,
                                        const shared_ptr<unordered_map<uint32_t, uint32_t>> &candidate_degree_map,
                                        const shared_ptr<unordered_map<uint32_t, shared_ptr<unordered_map<uint32_t, uint32_t> >>> &rcd,
                                        uint32_t n
    ) {
        auto u = e->get_source_vertex_id();
        auto v = e->get_destination_vertex_id();

        /**
         * @brief prepare skip map
         */
        auto B = make_shared<map<double, shared_ptr<double_node<uint32_t>>>>();

        /**
         * @brief deal with vertices are not contained in the graph
         */
        if (!k_order->count(1)) {
            k_order->insert({1, make_shared<double_list<uint32_t>>()});
        }
        if (!core->count(u)) {
            core->insert({u, 1});
            node_map->insert({u, make_shared<double_node<uint32_t>>(u)});

            k_order->at(1)->left_insert(node_map->at(u));
            tree->Insert(u,true,root->at(1));

            remaining_degree_map->insert({u, 0});
            candidate_degree_map->insert({u, 0});
        }
        if (!core->count(v)) {
            core->insert({v, 1});
            node_map->insert({v, make_shared<double_node<uint32_t>>(v)});

            k_order->at(1)->left_insert(node_map->at(v));
            tree->Insert(v,true,root->at(1));

            remaining_degree_map->insert({v, 0});
            candidate_degree_map->insert({v, 0});
        }


//        if (core->at(u) <= core->at(v)) {
//            if (!rcd->count(u)) {
//                rcd->insert({u, make_shared<unordered_map<uint32_t, uint32_t>>()});
//                rcd->at(u)->insert({1, 1});
//            }
//
//        }
//        if (core->at(v) <= core->at(u)) {
//            if (!rcd->count(v)) {
//                rcd->insert({v, make_shared<unordered_map<uint32_t, uint32_t>>()});
//                rcd->at(v)->insert({1, 1});
//            }
//        }
        /**
         * @brief update rcd value for u and v
         */
        auto K = min(core->at(u), core->at(v));
        G->insert_edge(e);

        multi_hop_prepare_rcd_insertion(G, core, rcd, n, e);

        if (!test_order(core, tree, u, v)) {
            swap(u, v);
        }

        ++remaining_degree_map->at(u);
        if (remaining_degree_map->at(u) <= K) {
            return;
        }

        B->insert({tree->Rank(u), node_map->at(u)});

        auto swap = make_shared<vector<uint32_t>>();// conserve the unaffected vertices
        auto k_list = k_order->at(K);
        auto new_k_list = make_shared<double_list<uint32_t>>();
        auto Vc = make_shared<list<uint32_t>>();
        auto Vc_set = make_shared<unordered_set<uint32_t>>();
        auto p = k_list->get_head();
        while (p) {
            auto vi = p->get_value();

            B->erase(tree->Rank(vi));
            if (candidate_degree_map->at(vi) + remaining_degree_map->at(vi) > K) {
                /**
                 * @brief conserve the next of p
                 */
                auto q = p->get_next();
                Vc_set->insert(vi);
                /**
                 * @remarks note the vertex order
                 */
                Vc->push_front(vi);

                auto vi_vertex = G->get_vertex(vi);
                for (const auto &[w,e]:*vi_vertex->get_edge_map()) {
                    if (core->at(w) == K && test_order(core, tree, vi, w)) {
                        ++candidate_degree_map->at(w);
                        /**
                         * @brief avoid repeatedly insert
                         */
                        if (!B->count(tree->Rank(w))) {
                            B->insert({tree->Rank(w), node_map->at(w)});
                        }
                    }
                }
                /**
                 * @brief move p to next node
                 */
                p = q;
            } else if (candidate_degree_map->at(vi) == 0) {
                shared_ptr<double_node<uint32_t>> q;
                shared_ptr<double_node<uint32_t>> prior_q;
                /**
                 * @details two cases: (1) get the prior node of the first item of B when B is nonempty (2) skip to the rear of this list
                 */
                if (!B->empty()) {
                    q = B->begin()->second;
                    prior_q = q->get_prior();
                } else {
                    prior_q = k_list->get_rear();
                }
                /**
                 * @brief move nodes from old list to the new list
                 */
                new_k_list->right_insert(p, prior_q);
                p = q;
            } else {
                auto q = p->get_next();
                new_k_list->right_insert(p);
                remaining_degree_map->at(vi) = remaining_degree_map->at(vi) + candidate_degree_map->at(vi);
                candidate_degree_map->at(vi) = 0;
                remove_candidates(G, Vc, Vc_set, node_map, tree,root, B, new_k_list, swap, vi, K, core,
                                  remaining_degree_map,
                                  candidate_degree_map);
                p = q;
            }
        }

        k_order->at(K) = new_k_list;

        for (auto w:*swap) {
            auto node = node_map->at(w);

            auto node_prior = node->get_prior();
            /**
             * @brief readjust the order of nodes who do not change their core number
             */
            tree->Delete(w,root->at(K));
            tree->InsertAfter(w,node->get_prior()->get_value(),root->at(K));
        }

        auto changed = make_shared<list<uint32_t>>();
        for (auto w:*Vc) {
            if (!Vc_set->count(w)) {
                continue;
            }

            if (!k_order->count(K + 1)) {
                k_order->insert({K + 1, make_shared<double_list<uint32_t>>()});
            }
            k_order->at(K + 1)->left_insert(node_map->at(w));

            candidate_degree_map->at(w) = 0;
            core->at(w) = K + 1;

            tree->Delete(w,root->at(K));
            tree->Insert(w, true,root->at(K+1));
            /**
             * @brief update maximal core degree
             */
            changed->push_back(w);
        }

        multi_hop_recompute_rcd(G, core, rcd, n, changed, true);
    }

    void order_core_maintenance::order_remove(const shared_ptr<abstract_graph> &G,
                                              const shared_ptr<abstract_edge> &e,
                                              const shared_ptr<unordered_map<uint32_t, uint32_t>> &core,
                                              const shared_ptr<unordered_map<uint32_t, shared_ptr<double_list<uint32_t>>>> &core_vertex_order,
                                              const shared_ptr<unordered_map<uint32_t, shared_ptr<double_node<uint32_t>>>> &node_map,
                                              const shared_ptr<unordered_map<uint32_t, shared_ptr<extend_list<double, uint32_t>>>> &core_order_query_map,
                                              const shared_ptr<unordered_map<uint32_t, uint32_t>> &remaining_degree_map,
                                              const shared_ptr<unordered_map<uint32_t,shared_ptr<unordered_map<uint32_t,uint32_t>>>>& rcd,
                                              uint32_t n) {
        auto u = e->get_source_vertex_id();
        auto v = e->get_destination_vertex_id();

        auto K = std::min(core->at(u), core->at(v));
        auto isolated_vertex_set = G->remove_edge(e);


        auto r = u;
        if (core->at(v) < core->at(r)) {
            r = v;
        }
        auto visited = make_shared<unordered_set<uint32_t>>();
        auto dismissed = make_shared<unordered_set<uint32_t>>();
        /**
         * @brief record the vertices in (k-1)core
         */
        auto Vc = make_shared<list<uint32_t>>() ;
        auto cd = make_shared<unordered_map<uint32_t, long long>>() ;
        auto k = core->at(r);

        multi_hop_prepare_rcd_removal(G, core, rcd, n, e);

        if (core->at(u) != core->at(v)) {
            visited->insert(r);
            //mcd
            cd->insert({r, rcd->at(r)->at(n)});
            if (cd->at(r) < k) {
                propagate_dismissal(G, core, rcd, cd, dismissed, visited, k, r, n, Vc);
            }
        } else {
            visited->insert(u);
            //mcd
            cd->insert({u, rcd->at(u)->at(n)});
            if (cd->at(u) < k) {
                propagate_dismissal(G, core, rcd, cd, dismissed, visited, k, u, n, Vc);
            }
            visited->insert(v);
            cd->insert({v, rcd->at(v)->at(n)});
            if (!dismissed->count(v) && cd->at(v) < k) {
                propagate_dismissal(G, core, rcd, cd, dismissed, visited, k, v, n, Vc);
            }
        }

        multi_hop_recompute_rcd(G, core, rcd, n, Vc, false);

        for (const auto& w:*Vc) {
            auto w_node = core_order_query_map->at(K)->remove(w);
            if (K > 1) {
                core_order_query_map->at(K - 1)->right_insert(w_node);
            }
        }

        for (const auto& w:*Vc) {
            remaining_degree_map->at(w) = 0;
            auto w_vertex = G->get_vertex(w);
            if (!w_vertex) {
                continue;
            }
            for (const auto& e_pair:*w_vertex->get_edge_map()) {
                auto w1 = e_pair.first;
                if (core->at(w1) == K && test_order(core, core_order_query_map, w1, w)) {
                    --remaining_degree_map->at(w);
                }
                if (core->at(w1) >= K) {
                    ++remaining_degree_map->at(w);
                }
            }
            if (K > 1) {
                core_vertex_order->at(k - 1)->right_insert(node_map->at(w));
            } else {
                core_vertex_order->at(K)->remove(w);
                node_map->erase(w);
            }
        }

        for(const auto&w:*isolated_vertex_set)
        {
            core->erase(w);
        }
    }

    void order_core_maintenance::remove(const shared_ptr<abstract_graph> &G,
                                        const shared_ptr<abstract_edge> &e,
                                        const shared_ptr<unordered_map<uint32_t, uint32_t>> &core,
                                        const shared_ptr<unordered_map<uint32_t, shared_ptr<double_list<uint32_t>>>> &core_vertex_order,
                                        const shared_ptr<unordered_map<uint32_t, shared_ptr<double_node<uint32_t>>>> &node_map,
                                        const shared_ptr<gadget::Treap> &tree,
                                        const shared_ptr<vector<long long>>& root,
                                        const shared_ptr<unordered_map<uint32_t, uint32_t>> &remaining_degree_map,
                                        const shared_ptr<unordered_map<uint32_t,shared_ptr<unordered_map<uint32_t,uint32_t>>>>& rcd,
                                        uint32_t n) {
        auto u = e->get_source_vertex_id();
        auto v = e->get_destination_vertex_id();

        auto K = std::min(core->at(u), core->at(v));
        auto isolated_vertex_set = G->remove_edge(e);


        auto r = u;
        if (core->at(v) < core->at(r)) {
            r = v;
        }
        auto visited = make_shared<unordered_set<uint32_t>>();
        auto dismissed = make_shared<unordered_set<uint32_t>>();
        /**
         * @brief record the vertices in (k-1)core
         */
        auto Vc = make_shared<list<uint32_t>>() ;
        auto cd = make_shared<unordered_map<uint32_t, long long>>() ;
        auto k = core->at(r);

        multi_hop_prepare_rcd_removal(G, core, rcd, n, e);

        if (core->at(u) != core->at(v)) {
            visited->insert(r);
            //mcd
            cd->insert({r, rcd->at(r)->at(n)});
            if (cd->at(r) < k) {
                propagate_dismissal(G, core, rcd, cd, dismissed, visited, k, r, n, Vc);
            }
        } else {
            visited->insert(u);
            //mcd
            cd->insert({u, rcd->at(u)->at(n)});
            if (cd->at(u) < k) {
                propagate_dismissal(G, core, rcd, cd, dismissed, visited, k, u, n, Vc);
            }
            visited->insert(v);
            cd->insert({v, rcd->at(v)->at(n)});
            if (!dismissed->count(v) && cd->at(v) < k) {
                propagate_dismissal(G, core, rcd, cd, dismissed, visited, k, v, n, Vc);
            }
        }

        multi_hop_recompute_rcd(G, core, rcd, n, Vc, false);

        for (const auto& w:*Vc) {
            tree->Delete(w,root->at(K));
            tree->Insert(w, false,root->at(K-1));
        }

        for (const auto& w:*Vc) {
            remaining_degree_map->at(w) = 0;
            auto w_vertex = G->get_vertex(w);
            if (!w_vertex) {
                continue;
            }
            for (const auto& e_pair:*w_vertex->get_edge_map()) {
                auto w1 = e_pair.first;
                if (core->at(w1) == K && test_order(core, tree, w1, w)) {
                    --remaining_degree_map->at(w);
                }
                if (core->at(w1) >= K) {
                    ++remaining_degree_map->at(w);
                }
            }
            if (K > 1) {
                core_vertex_order->at(k - 1)->right_insert(node_map->at(w));
            } else {
                core_vertex_order->at(K)->remove(w);
                node_map->erase(w);
            }
        }

        for(const auto&w:*isolated_vertex_set)
        {
            core->erase(w);
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
    void order_core_maintenance::propagate_dismissal(const shared_ptr<abstract_graph> &G,
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
     * @details remove unsatisfied vertices from the candidate vertex set
     * @param G
     * @param Vc
     * @param Vc_set
     * @param node_map
     * @param order_query_map
     * @param B
     * @param new_k_list
     * @param swap
     * @param w
     * @param K
     * @param core
     * @param remaining_degree_map
     * @param candidate_degree_map
     */
    void order_core_maintenance::remove_candidates(const shared_ptr<abstract_graph> &G,
                                                   const shared_ptr<list<uint32_t>> &Vc,
                                                   const shared_ptr<unordered_set<uint32_t>> &Vc_set,
                                                   const shared_ptr<unordered_map<uint32_t, shared_ptr<double_node<uint32_t>>>> &node_map,
                                                   const shared_ptr<unordered_map<uint32_t, shared_ptr<extend_list<double, uint32_t>>>> &order_query_map,
                                                   const shared_ptr<map<double, shared_ptr<double_node<uint32_t>>>> &B,
                                                   const shared_ptr<double_list<uint32_t>> &new_k_list,
                                                   const shared_ptr<vector<uint32_t>> &swap,
                                                   uint32_t w,
                                                   uint32_t K,
                                                   const shared_ptr<unordered_map<uint32_t, uint32_t>> &core,
                                                   const shared_ptr<unordered_map<uint32_t, uint32_t>> &remaining_degree_map,
                                                   const shared_ptr<unordered_map<uint32_t, uint32_t>> &candidate_degree_map
    ) {
        queue<uint32_t> Q;
        unordered_set<uint32_t> visited;
        auto w_vertex = G->get_vertex(w);
        for (const auto &[w1,e]:*w_vertex->get_edge_map()) {
            if (Vc_set->count(w1)) {
                //w1 behind w
                --remaining_degree_map->at(w1);
                if (remaining_degree_map->at(w1) + candidate_degree_map->at(w1) <= K) {
                    Q.push(w1);
                    visited.insert(w1);
                }
            }
        }
        auto pivot_node = order_query_map->at(K)->find(w);
        while (!Q.empty()) {
            auto w1 = Q.front();
            Q.pop();

            remaining_degree_map->at(w1) += candidate_degree_map->at(w1);
            candidate_degree_map->at(w1) = 0;

            //remove w1 from Vc
            Vc_set->erase(w1);
            swap->push_back(w1);

            new_k_list->right_insert(node_map->at(w1));
            B->erase(order_query_map->at(K)->find_key(w1).value());

            auto w1_vertex = G->get_vertex(w1);
            for (const auto &[w2,e]:*w1_vertex->get_edge_map()) {
                if (core->at(w2) == K) {
                    if (test_order(core, order_query_map, w, w2)) {
                        --candidate_degree_map->at(w2);
                        if (candidate_degree_map->at(w2) == 0) {
                            B->erase(order_query_map->at(K)->find_key(w2).value());
                        }
                    } else if (test_order(core, order_query_map, w1, w2) && Vc_set->count(w2)) {
                        --candidate_degree_map->at(w2);
                        if (candidate_degree_map->at(w2) + remaining_degree_map->at(w2) <= K
                            && !visited.count(w2)) {
                            Q.push(w2);
                            visited.insert(w2);
                        }
                    } else if (Vc_set->count(w2)) {
                        --remaining_degree_map->at(w2);
                        if (candidate_degree_map->at(w2) + remaining_degree_map->at(w2) <= K
                            && !visited.count(w2)) {
                            Q.push(w2);
                            visited.insert(w2);
                        }
                    }
                }
            }
        }
    }

    /**
     * @details remove unsatisfied vertices from the candidate vertex set
     * @param G
     * @param Vc
     * @param Vc_set
     * @param node_map
     * @param tree
     * @param B
     * @param new_k_list
     * @param swap
     * @param w
     * @param K
     * @param core
     * @param remaining_degree_map
     * @param candidate_degree_map
     */
    void order_core_maintenance::remove_candidates(const shared_ptr<abstract_graph> &G,
                                                   const shared_ptr<list<uint32_t>> &Vc,
                                                   const shared_ptr<unordered_set<uint32_t>> &Vc_set,
                                                   const shared_ptr<unordered_map<uint32_t, shared_ptr<double_node<uint32_t>>>> &node_map,
                                                   const shared_ptr<gadget::Treap> &tree,
                                                   const shared_ptr<vector<long long>> &root,
                                                   const shared_ptr<map<double, shared_ptr<double_node<uint32_t>>>> &B,
                                                   const shared_ptr<double_list<uint32_t>> &new_k_list,
                                                   const shared_ptr<vector<uint32_t>> &swap,
                                                   uint32_t w,
                                                   uint32_t K,
                                                   const shared_ptr<unordered_map<uint32_t, uint32_t>> &core,
                                                   const shared_ptr<unordered_map<uint32_t, uint32_t>> &remaining_degree_map,
                                                   const shared_ptr<unordered_map<uint32_t, uint32_t>> &candidate_degree_map
    ) {
        queue<uint32_t> Q;
        unordered_set<uint32_t> visited;
        auto w_vertex = G->get_vertex(w);
        for (const auto &[w1,e]:*w_vertex->get_edge_map()) {
            if (Vc_set->count(w1)) {
                //w1 behind w
                --remaining_degree_map->at(w1);
                if (remaining_degree_map->at(w1) + candidate_degree_map->at(w1) <= K) {
                    Q.push(w1);
                    visited.insert(w1);
                }
            }
        }
        while (!Q.empty()) {
            auto w1 = Q.front();
            Q.pop();

            remaining_degree_map->at(w1) += candidate_degree_map->at(w1);
            candidate_degree_map->at(w1) = 0;

            //remove w1 from Vc
            Vc_set->erase(w1);
            swap->push_back(w1);

            new_k_list->right_insert(node_map->at(w1));
            B->erase(tree->Rank(w1));

            auto w1_vertex = G->get_vertex(w1);
            for (const auto &[w2,e]:*w1_vertex->get_edge_map()) {
                if (core->at(w2) == K) {
                    if (test_order(core, tree, w, w2)) {
                        --candidate_degree_map->at(w2);
                        if (candidate_degree_map->at(w2) == 0) {
                            B->erase(tree->Rank(w2));
                        }
                    } else if (test_order(core, tree, w1, w2) && Vc_set->count(w2)) {
                        --candidate_degree_map->at(w2);
                        if (candidate_degree_map->at(w2) + remaining_degree_map->at(w2) <= K
                            && !visited.count(w2)) {
                            Q.push(w2);
                            visited.insert(w2);
                        }
                    } else if (Vc_set->count(w2)) {
                        --remaining_degree_map->at(w2);
                        if (candidate_degree_map->at(w2) + remaining_degree_map->at(w2) <= K
                            && !visited.count(w2)) {
                            Q.push(w2);
                            visited.insert(w2);
                        }
                    }
                }
            }
        }
    }



    bool order_core_maintenance::test_order(const shared_ptr<unordered_map<uint32_t, uint32_t>> &core,
                                            const shared_ptr<unordered_map<uint32_t, shared_ptr<extend_list<double, uint32_t>>>> &order_query_map,
                                            uint32_t u,
                                            uint32_t v) {
        if (core->at(u) < core->at(v)) {
            return true;
        } else if (core->at(u) > core->at(v)) {
            return false;
        } else {
            auto k = core->at(u);
            return order_query_map->at(k)->find_key(u).value() < order_query_map->at(k)->find_key(v).value();
        }
    }

    bool order_core_maintenance::test_order(const shared_ptr<unordered_map<uint32_t, uint32_t>> &core,
                                            const shared_ptr<gadget::Treap> &tree,
                                            uint32_t u,
                                            uint32_t v) {
        if (core->at(u) < core->at(v)) {
            return true;
        } else if (core->at(u) > core->at(v)) {
            return false;
        } else {
            auto k = core->at(u);
            return tree->Rank(u) < tree->Rank(v);
        }
    }
}