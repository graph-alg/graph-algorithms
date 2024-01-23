
#include "wing/BE_index.h"

namespace scnu{

    BE_index::BE_index(): bloom_map(make_shared<unordered_map<pair<uint32_t,uint32_t>,shared_ptr<priority_obeyed_bloom>,hash_pair,equal_pair>>()),
                          edge_support_map(make_shared<unordered_map<shared_ptr<abstract_bipartite_edge>,uint32_t>>()),
                          edge_bloom_map(make_shared<unordered_map<shared_ptr<abstract_bipartite_edge>,shared_ptr<unordered_set<shared_ptr<priority_obeyed_bloom>>>>>()){

    }

    bool BE_index::count(uint32_t u,uint32_t v){
        return bloom_map->count({u,v});
    }

    bool BE_index::count(const shared_ptr<abstract_bipartite_edge>& e){
        return edge_support_map->count(e);
    }

    uint32_t BE_index::get_support(const shared_ptr<abstract_bipartite_edge>& e){
        return edge_support_map->count(e) ? edge_support_map->at(e):0;
    }

    shared_ptr<unordered_map<shared_ptr<abstract_bipartite_edge>,uint32_t>> BE_index::get_edge_support_map()
    {
        return edge_support_map;
    }

    shared_ptr<priority_obeyed_bloom> BE_index::get_bloom(uint32_t u,uint32_t v)
    {
        return bloom_map->count({u,v}) ? bloom_map->at({u,v}):shared_ptr<priority_obeyed_bloom>();
    }

    shared_ptr<unordered_map<pair<uint32_t,uint32_t>,shared_ptr<priority_obeyed_bloom>,hash_pair,equal_pair>>
    BE_index::get_bloom_map()
    {
        return bloom_map;
    }

    shared_ptr<unordered_set<shared_ptr<priority_obeyed_bloom>>> BE_index::get_bloom_set(const shared_ptr<abstract_bipartite_edge>& e){
        return edge_bloom_map->count(e)? edge_bloom_map->at(e):shared_ptr<unordered_set<shared_ptr<priority_obeyed_bloom>>>();
    }

    void BE_index::insert_bloom(const shared_ptr<priority_obeyed_bloom>& B)
    {
        bloom_map->insert({B->get_vertex_pair(),B});
    }

    void BE_index::link_bloom(const shared_ptr<abstract_bipartite_edge>& e,
                              const shared_ptr<priority_obeyed_bloom>& B){
        edge_bloom_map->at(e)->insert(B);
    }

    void BE_index::insert_edge(const shared_ptr<abstract_bipartite_edge>& e,
                     uint32_t support){
        edge_bloom_map->insert({e,make_shared<unordered_set<shared_ptr<priority_obeyed_bloom>>>()});
        edge_support_map->insert({e,support});
    }

    void BE_index::removal_bloom(const shared_ptr<priority_obeyed_bloom>& B){
        for(const auto&[e1,e2]:*B->get_edge_map()){
            edge_bloom_map->erase(e1);
        }
        bloom_map->erase(B->get_vertex_pair());
    }

    void BE_index::remove_edge(const shared_ptr<abstract_bipartite_edge>& e){
        if(edge_bloom_map->count(e)){
            for(const auto&B:*edge_bloom_map->at(e))
            {
                B->remove_edge(e);
            }
        }
        edge_bloom_map->erase(e);
        edge_support_map->erase(e);
    }

    void BE_index::update_support(const shared_ptr<abstract_bipartite_edge>& e,
                                  uint32_t value)
    {
        edge_support_map->at(e) = value;
    }

    void BE_index::remove_edge(const shared_ptr<abstract_bipartite_edge>& e,
                               const shared_ptr<priority_obeyed_bloom>& B)
    {
        edge_bloom_map->at(e)->erase(B);
    }
}

