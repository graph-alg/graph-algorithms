
#include "core/bin_sort_core_decomposition.h"

namespace scnu
{
    /**
    * @details bin sort core decomposition to get core number and core vertex map
    * @param graph
    * @param vertex_core_map
    * @return
    */
    uint32_t bin_sort_core_decomposition::decompose(const shared_ptr<abstract_graph> &graph,
                                                    const shared_ptr<unordered_map<uint32_t,uint32_t>> &vertex_core_map) {
        auto vertex_degree_map = make_shared<unordered_map<uint32_t,uint32_t>>();
        auto degree_vertex_map = make_shared<map<uint32_t,shared_ptr<unordered_set<uint32_t>>>>();
        for (const auto &[u,u_vertex]:*graph->get_vertex_map()) {
            auto degree = u_vertex->get_degree();
            vertex_degree_map->insert({u,degree});
            if (!degree_vertex_map->count(degree)) {
                degree_vertex_map->insert({degree, make_shared<unordered_set<uint32_t>>()});
            }
            degree_vertex_map->at(degree)->insert(u);
        }

        uint32_t k_max = 1;
        uint32_t k = 1;

        while(!vertex_degree_map->empty()){
            k_max = k;
            if(!degree_vertex_map->count(k)){
                ++k;
                continue;
            }
            auto vertex_set = degree_vertex_map->at(k);
            while(!vertex_set->empty()){
                auto u = *vertex_set->begin();
                vertex_set->erase(u);
                vertex_core_map->insert({u,k});

                for(const auto&[v,e]:*graph->get_vertex(u)->get_edge_map()){
                    if(!vertex_degree_map->count(v) || vertex_degree_map->at(v) <= k){
                        continue;
                    }

                    degree_vertex_map->at(vertex_degree_map->at(v))->erase(v);
                    --vertex_degree_map->at(v);
                    if(vertex_degree_map->at(v) <= k){
                        vertex_set->insert(v);
                    }else
                    {
                        if(!degree_vertex_map->count(vertex_degree_map->at(v))){
                            degree_vertex_map->insert({vertex_degree_map->at(v), make_shared<unordered_set<uint32_t>>()});
                        }
                        degree_vertex_map->at(vertex_degree_map->at(v))->insert(v);
                    }
                }

                vertex_degree_map->erase(u);
            }
            ++k;
        }
        return k_max;
    }
}

