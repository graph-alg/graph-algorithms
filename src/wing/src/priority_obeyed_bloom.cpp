
#include "wing/priority_obeyed_bloom.h"

namespace scnu{
    priority_obeyed_bloom::priority_obeyed_bloom(uint32_t u, uint32_t w):
                                                 vertex_pair({u,w}),
                                                 butterfly_count(0),
                                                 edge_twin_map(make_shared<unordered_map<shared_ptr<abstract_bipartite_edge>,
                                                            shared_ptr<abstract_bipartite_edge>>>())
    {

    }

    uint32_t priority_obeyed_bloom::count(const shared_ptr<abstract_bipartite_edge>& e){
        return edge_twin_map->count(e);
    }

    bool priority_obeyed_bloom::empty()
    {
        return edge_twin_map->empty();
    }

    /**
     * @details get the total butterfly count of this bloom
     * @return
     */
    uint32_t priority_obeyed_bloom::get_butterfly_count() const {
        return butterfly_count;
    }

    shared_ptr<unordered_map<shared_ptr<abstract_bipartite_edge>, shared_ptr<abstract_bipartite_edge>>>
    priority_obeyed_bloom::get_edge_map()
    {
        return edge_twin_map;
    }

    /**
     * @details compute the number of another type of vertices in this bloom (e.g., (2, k)-clique)
     * @return
     */
    uint32_t priority_obeyed_bloom::get_k() const{
         return static_cast<uint32_t>(ceil(sqrt(butterfly_count*2)));
    }

    /**
     * @details get the twin of the given edge
     * @param e
     * @return
     */
    shared_ptr<abstract_bipartite_edge> priority_obeyed_bloom::get_twin(const shared_ptr<abstract_bipartite_edge> &e) {
        return edge_twin_map->count(e) ? edge_twin_map->at(e) : shared_ptr<abstract_bipartite_edge>();
    }

    pair<uint32_t,uint32_t> priority_obeyed_bloom::get_vertex_pair(){
        return vertex_pair;
    }

    void priority_obeyed_bloom::remove_edge(const shared_ptr<abstract_bipartite_edge>&e){
        edge_twin_map->erase(e);
    }

    void priority_obeyed_bloom::set_butterfly_count(uint32_t value){
        butterfly_count = value;
    }

    void priority_obeyed_bloom::link_twin(const shared_ptr<abstract_bipartite_edge> &e,
                                          const shared_ptr<abstract_bipartite_edge> &twin_edge){
        edge_twin_map->insert({e,twin_edge});
        edge_twin_map->insert({twin_edge,e});
    }
}

