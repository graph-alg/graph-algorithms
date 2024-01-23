
#include "truss/bin_sort_truss_decomposition.h"

namespace scnu{

    void bin_sort_truss_decomposition::init(const shared_ptr<abstract_graph> &G,
                                            const shared_ptr<unordered_map<shared_ptr<abstract_edge>, uint32_t>> &edge_support_map,
                                            const shared_ptr<unordered_map<shared_ptr<abstract_edge>, uint32_t>> &edge_truss_map) {
        auto edge_set = G->get_edge_set();
        for(const auto &e:*edge_set){
            edge_support_map->insert({e, 0});
            edge_truss_map->insert({e, 0});
        }
    }

    void bin_sort_truss_decomposition::init(const shared_ptr<abstract_graph> &G,
                                            const shared_ptr<unordered_map<shared_ptr<abstract_edge>, uint32_t>> &edge_support_map,
                                            const shared_ptr<unordered_map<shared_ptr<abstract_edge>, uint32_t>> &edge_truss_map,
                                            const shared_ptr<unordered_map<shared_ptr<abstract_edge>, uint32_t>> &rem) {
        auto edge_set = G->get_edge_set();
        for(const auto &e:*edge_set){
            edge_support_map->insert({e, 0});
            edge_truss_map->insert({e, 0});
            rem->insert({e, 0});
        }
    }

    uint32_t bin_sort_truss_decomposition::decompose(const shared_ptr<abstract_graph>& G,
                                                     const shared_ptr<unordered_map<shared_ptr<abstract_edge>, uint32_t>>& edge_support_map,
                                                     const shared_ptr<unordered_map<shared_ptr<abstract_edge>, uint32_t>>& edge_truss_map) {
        auto edge_set = G->get_edge_set();
        edge_support_computation(G,edge_set, edge_support_map);

        auto support_edge_map = make_shared<unordered_map<uint32_t, shared_ptr<unordered_set<shared_ptr<abstract_edge>>>>>();
        for (const auto &[e, e_support]:*edge_support_map) {
            if(!support_edge_map->count(e_support)){
                support_edge_map->insert({e_support, make_shared<unordered_set<shared_ptr<abstract_edge>>>()});
            }
            support_edge_map->at(e_support)->insert(e);
        }
        uint32_t k = 2;
        uint32_t k_max = 2;
        while (!edge_set->empty()) {
            k_max = k;
            if(!support_edge_map->count(k-2)){
                ++k;
                continue;
            }
            auto evicted_set = support_edge_map->at(k-2);
            while (!evicted_set->empty()) {
                auto e1 = *evicted_set->begin();
                evicted_set->erase(e1);
                edge_set->erase(e1);

                edge_truss_map->at(e1) = k;

                auto u = e1->get_source_vertex_id();
                auto v = e1->get_destination_vertex_id();
                if(G->get_vertex(u)->get_degree() > G->get_vertex(v)->get_degree()){
                    swap(u,v);
                }
                for(const auto&[w,e2]:* G->get_vertex(u)->get_edge_map()){
                    if(w == v ||!edge_set->count(e2))
                    {
                        continue;
                    }
                    auto e3 = G->get_edge(v, w);
                    if(!e3 || !edge_set->count(e3))
                    {
                        continue;
                    }

                    if(edge_support_map->at(e2) > k-2) {
                        support_edge_map->at(edge_support_map->at(e2))->erase(e2);
                        --edge_support_map->at(e2);

                        if (edge_support_map->at(e2) == k - 2) {
                            evicted_set->insert(e2);
                        } else {
                            if (!support_edge_map->count(edge_support_map->at(e2))) {
                                support_edge_map->insert({edge_support_map->at(e2),
                                                          make_shared<unordered_set<shared_ptr<abstract_edge>>>()});
                            }
                            support_edge_map->at(edge_support_map->at(e2))->insert(e2);
                        }
                    }

                    if (edge_support_map->at(e3) > k - 2) {
                        support_edge_map->at(edge_support_map->at(e3))->erase(e3);
                        --edge_support_map->at(e3);
                        if (edge_support_map->at(e3) == k - 2) {
                            evicted_set->insert(e3);
                        } else {
                            if (!support_edge_map->count(edge_support_map->at(e3))) {
                                support_edge_map->insert({edge_support_map->at(e3),
                                                          make_shared<unordered_set<shared_ptr<abstract_edge>>>()});
                            }
                            support_edge_map->at(edge_support_map->at(e3))->insert(e3);
                        }
                    }
                }
            }
            ++k;
        }
        return k_max;
    }

    uint32_t bin_sort_truss_decomposition::decompose(const shared_ptr<abstract_graph>& G,
                                                     const shared_ptr<unordered_map<shared_ptr<abstract_edge>, uint32_t>>& edge_support_map,
                                                     const shared_ptr<unordered_map<shared_ptr<abstract_edge>, uint32_t>>& edge_truss_map,
                                                     const shared_ptr<unordered_map<uint32_t, shared_ptr<extend_list<double, shared_ptr<abstract_edge>>>>>& truss_order_map,
                                                     const shared_ptr<unordered_map<shared_ptr<abstract_edge>, uint32_t>>& rem) {
        auto edge_set = G->get_edge_set();
        edge_support_computation(G,edge_set, edge_support_map);

        auto support_edge_map = make_shared<unordered_map<uint32_t, shared_ptr<unordered_set<shared_ptr<abstract_edge>>>>>();
        for (const auto &[e, e_support]:*edge_support_map) {
            if(!support_edge_map->count(e_support)){
                support_edge_map->insert({e_support, make_shared<unordered_set<shared_ptr<abstract_edge>>>()});
            }
            support_edge_map->at(e_support)->insert(e);
        }
        uint32_t k = 2;
        uint32_t k_max = 2;
        while (!edge_set->empty()) {
            k_max = k;
            if(!support_edge_map->count(k-2)){
                ++k;
                continue;
            }
            if(!truss_order_map->count(k)){
                truss_order_map->insert({k, make_shared<extend_list<double, shared_ptr<abstract_edge>>>()});
            }
            auto evicted_set = support_edge_map->at(k-2);
            while (!evicted_set->empty()) {
                auto e1 = *evicted_set->begin();
                evicted_set->erase(e1);
                edge_set->erase(e1);

                edge_truss_map->at(e1) = k;
                truss_order_map->at(k)->push_back(e1);

                auto u = e1->get_source_vertex_id();
                auto v = e1->get_destination_vertex_id();
                if(G->get_vertex(u)->get_degree() > G->get_vertex(v)->get_degree()){
                    swap(u,v);
                }
                for(const auto&[w,e2]:* G->get_vertex(u)->get_edge_map()){
                    if(w == v ||!edge_set->count(e2))
                    {
                        continue;
                    }
                    auto e3 = G->get_edge(v, w);
                    if(!e3 || !edge_set->count(e3))
                    {
                        continue;
                    }

                    ++rem->at(e1);

                    if(edge_support_map->at(e2) > k-2){
                        support_edge_map->at(edge_support_map->at(e2))->erase(e2);
                        --edge_support_map->at(e2);

                        if (edge_support_map->at(e2) == k - 2) {
                            evicted_set->insert(e2);
                        } else {
                            if (!support_edge_map->count(edge_support_map->at(e2))) {
                                support_edge_map->insert({edge_support_map->at(e2),
                                                          make_shared<unordered_set<shared_ptr<abstract_edge>>>()});
                            }
                            support_edge_map->at(edge_support_map->at(e2))->insert(e2);
                        }
                    }


                    if(edge_support_map->at(e3) > k-2){
                        support_edge_map->at(edge_support_map->at(e3))->erase(e3);
                        --edge_support_map->at(e3);
                        if (edge_support_map->at(e3) == k - 2) {
                            evicted_set->insert(e3);
                        } else {
                            if (!support_edge_map->count(edge_support_map->at(e3))) {
                                support_edge_map->insert({edge_support_map->at(e3),
                                                          make_shared<unordered_set<shared_ptr<abstract_edge>>>()});
                            }
                            support_edge_map->at(edge_support_map->at(e3))->insert(e3);
                        }
                    }
                }
            }
            ++k;
        }
        return k_max;
    }


    void bin_sort_truss_decomposition::edge_support_computation(const shared_ptr<abstract_graph> &G,
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
                if(!e3||visited_set->count(e3))
                {
                    continue;
                }

                ++edge_support_map->at(e1);
                ++edge_support_map->at(e2);
                ++edge_support_map->at(e3);
            }
        }
    }
}