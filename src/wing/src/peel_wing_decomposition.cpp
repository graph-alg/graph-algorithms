
#include "wing/peel_wing_decomposition.h"

namespace scnu
{
    /**
     * @details decompose a given G
     * @param G
     * @param edge_wing_map
     * @param wing_edge_map
     * @return
     */
    uint32_t peel_wing_decomposition::decompose(const std::shared_ptr<abstract_bipartite_graph> &G,
                                                const std::shared_ptr<unordered_map<shared_ptr<abstract_bipartite_edge>, uint32_t>> &edge_wing_map,
                                                const std::shared_ptr<unordered_map<uint32_t, shared_ptr<unordered_set<shared_ptr<abstract_bipartite_edge>>>>> &wing_edge_map) {
        auto edge_set = G->get_edge_set();

        auto vertex_priority_map = vertex_priority_computation(G);

        auto edge_support_map  = edge_support_computation(G,edge_set,vertex_priority_map);
        auto support_edge_map = make_shared<unordered_map<uint32_t,shared_ptr<unordered_set<shared_ptr<abstract_bipartite_edge>>>>>();

        uint32_t max_support = 0;
        for(const auto&[e,support]:*edge_support_map){
            if(support > max_support){
                max_support = support;
            }
            if(!support_edge_map->count(support)){
                support_edge_map->insert({support, make_shared<unordered_set<shared_ptr<abstract_bipartite_edge>>>()});
            }
            support_edge_map->at(support)->insert(e);
        }
        for(uint32_t support = 0;support <= max_support;++support){
            if(!support_edge_map->count(support)){
                support_edge_map->insert({support, make_shared<unordered_set<shared_ptr<abstract_bipartite_edge>>>()});
            }
        }

        uint32_t k = 0;
        uint32_t k_max = 0;
        while (!edge_set->empty()) {
            auto current_edge_set = support_edge_map->at(k);
            if(!current_edge_set->empty())
            {
                k_max = k;
            }
            while (!current_edge_set->empty()) {
                // e1 = (l1, r1)
                auto e1 = *current_edge_set->begin();
                current_edge_set->erase(e1);

                edge_set->erase(e1);
                edge_wing_map->insert({e1, k});

                auto l1 = e1->get_left_vertex_id();
                auto r1 = e1->get_right_vertex_id();

                auto l1_vertex = G->get_left_vertex(l1);
                auto r1_vertex = G->get_right_vertex(r1);

                //e2 = (l1,r2)
                for (const auto&[r2, e2]: *l1_vertex->get_edge_map()) {
                    if (r2 == r1 || !edge_set->count(e2)) {
                        continue;
                    }

                    //e3 = (l2,r1)
                    for (const auto&[l2, e3] : *r1_vertex->get_edge_map()) {
                        if (l2 == l1 || !edge_set->count(e3)) {
                            continue;
                        }

                        //e4 = (l2, r2)
                        auto e4 = G->get_edge(l2, r2);
                        if (!e4 || !edge_set->count(e4)) {
                            continue;
                        }

                        support_edge_map->at(edge_support_map->at(e2))->erase(e2);
                        edge_support_map->at(e2) = edge_support_map->at(e2) > k ? edge_support_map->at(e2) - 1 : k;
                        if (edge_support_map->at(e2) == k) {
                            current_edge_set->insert(e2);
                        } else {
                            support_edge_map->at(edge_support_map->at(e2))->insert(e2);
                        }

                        support_edge_map->at(edge_support_map->at(e3))->erase(e3);
                        edge_support_map->at(e3) = edge_support_map->at(e3) > k ? edge_support_map->at(e3) - 1 : k;
                        if (edge_support_map->at(e3) == k) {
                            current_edge_set->insert(e3);
                        } else {
                            support_edge_map->at(edge_support_map->at(e3))->insert(e3);
                        }


                        support_edge_map->at(edge_support_map->at(e4))->erase(e4);
                        edge_support_map->at(e4) = edge_support_map->at(e4) > k ? edge_support_map->at(e4) - 1 : k;

                        if (edge_support_map->at(e4) == k) {
                            current_edge_set->insert(e4);
                        } else {
                            support_edge_map->at(edge_support_map->at(e4))->insert(e4);
                        }
                    }
                }
            }
            ++k;
        }
        return k_max;
    }


    /**
     * @details compute the support (the number of butterflies) of each edge in the given graph
     * @param G
     * @param vertex_priority_map
     * @param thread_count
     * @return
     */
    shared_ptr<unordered_map<shared_ptr<abstract_bipartite_edge>,uint32_t>>
    peel_wing_decomposition::edge_support_computation(const shared_ptr<abstract_bipartite_graph> &G,
                                                      const shared_ptr<unordered_set<shared_ptr<abstract_bipartite_edge>>> &edge_set,
                                                      const shared_ptr<unordered_map<uint32_t, uint32_t>> &vertex_priority_map) {
        auto edge_support_map = make_shared<unordered_map<shared_ptr<abstract_bipartite_edge>,uint32_t>>();
        for (const auto &e:*edge_set) {
            edge_support_map->insert({e, 0});
        }

        /**
         * @brief compute wedges from left vertices
         */
        for (const auto &p:*G->get_left_vertex_map()) {
            auto&[l1, l1_vertex] = p;
            auto wedge_count_map = make_shared<unordered_map<uint32_t, uint32_t>>();

            for (const auto&[r1, e]:*l1_vertex->get_edge_map()) {
                if (vertex_priority_map->at(r1) < vertex_priority_map->at(l1)) {
                    auto r1_vertex = G->get_right_vertex(r1);
                    for (const auto&[l2, e]:*r1_vertex->get_edge_map()) {
                        if (vertex_priority_map->at(l2) < vertex_priority_map->at(l1)) {
                            if (!wedge_count_map->count(l2)) {
                                wedge_count_map->insert({l2, 0});
                            }
                            ++wedge_count_map->at(l2);
                        }
                    }
                }
            }
            for (const auto&[r1, l1r1_edge]:*l1_vertex->get_edge_map()) {
                if (vertex_priority_map->at(r1) < vertex_priority_map->at(l1)) {
                    auto r1_vertex = G->get_right_vertex(r1);
                    for (const auto&[l2, l2r1_edge]:*r1_vertex->get_edge_map()) {
                        if (vertex_priority_map->at(l2) < vertex_priority_map->at(l1)) {
                            if (wedge_count_map->at(l2) > 1) {
                                auto delta = wedge_count_map->at(l2) - 1;

                                edge_support_map->at(l1r1_edge) += delta;
                                edge_support_map->at(l2r1_edge) += delta;
                            }
                        }
                    }
                }
            }
        }
        /**
         * @brief compute wedges from right vertices
         */
        for (const auto &p:*G->get_right_vertex_map()) {
            auto&[r1, r1_vertex] = p;
            auto wedge_count_map = make_shared<unordered_map<uint32_t, uint32_t>>();

            for (const auto&[l1, l1r1_edge]:*r1_vertex->get_edge_map()) {
                if (vertex_priority_map->at(l1) < vertex_priority_map->at(r1)) {
                    auto l1_vertex = G->get_left_vertex(l1);
                    for (const auto&[r2, l1r2_edge]:*l1_vertex->get_edge_map()) {
                        if (vertex_priority_map->at(r2) < vertex_priority_map->at(r1)) {
                            if (!wedge_count_map->count(r2)) {
                                wedge_count_map->insert({r2, 0});
                            }
                            ++wedge_count_map->at(r2);
                        }
                    }
                }
            }
            for (const auto&[l1, l1r1_edge]:*r1_vertex->get_edge_map()) {
                if (vertex_priority_map->at(l1) < vertex_priority_map->at(r1)) {
                    auto l1_vertex = G->get_left_vertex(l1);
                    for (const auto&[r2, l1r2_edge]:*l1_vertex->get_edge_map()) {
                        if (vertex_priority_map->at(r2) < vertex_priority_map->at(r1)) {
                            if (wedge_count_map->at(r2) > 1) {
                                auto delta = wedge_count_map->at(r2) - 1;

                                edge_support_map->at(l1r1_edge) += delta;

                                edge_support_map->at(l1r2_edge) += delta;
                            }
                        }
                    }
                }
            }
        }

        return edge_support_map;
    }


    shared_ptr<unordered_map<uint32_t,uint32_t>>
    peel_wing_decomposition::vertex_priority_computation(const shared_ptr<abstract_bipartite_graph> &G) {
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
        std::sort(pair_vector->begin(), pair_vector->end(), pair_compare);
        auto p = make_shared<unordered_map<uint32_t, uint32_t>>();
        for (uint32_t i = 0; i < pair_vector->size(); ++i) {
            p->insert({pair_vector->at(i).first, i + 1});
        }
        return p;
    }
}

