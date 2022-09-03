
#include "graph/temporal_left_vertex.h"

namespace scnu {
    /**
     * @details construct a vertex with the given left vertex id
     * @param other_l
     */
    temporal_left_vertex::temporal_left_vertex(uint32_t other_l) :
            left_vertex_id(other_l),
            edge_map(make_shared<unordered_map<uint32_t, shared_ptr<unordered_set<shared_ptr<temporal_bipartite_edge>>>>>()) {

    }

    /**
     * @details construct a vertex with the given left vertex
     * @param l_vertex
     */
    temporal_left_vertex::temporal_left_vertex(const shared_ptr<temporal_left_vertex> &l_vertex) :
            temporal_left_vertex(l_vertex->left_vertex_id) {
        auto other_edge_map = l_vertex->get_edge_map();
        for(const auto &[r,r_set]:*other_edge_map){
            edge_map->insert({r, make_shared<unordered_set<shared_ptr<temporal_bipartite_edge>>>()});
            copy(r_set->begin(), r_set->end(), inserter(*edge_map->at(r),edge_map->at(r)->end()));
        }
    }

    /**
     * @details get the id of this left vertex
     * @return
     */
    uint32_t temporal_left_vertex::get_left_vertex_id() const {
        return left_vertex_id;
    }

    /**
     * @details get the degree of this left vertex
     * @return
     */
    uint32_t temporal_left_vertex::get_neighbor_size() {
        return edge_map->size();
    }

    /**
     * @details get the adjacent edge with the given right vertex id
     * @param r
     * @return
     */
    shared_ptr<unordered_set<shared_ptr<temporal_bipartite_edge>>> temporal_left_vertex::get_edge_set(uint32_t r) {
        if (edge_map->count(r)) {
            return edge_map->at(r);
        }
        return shared_ptr<unordered_set<shared_ptr<temporal_bipartite_edge>>>();
    }

    /**
     * @details get the adjacent edge map
     * @return
     */
    shared_ptr<unordered_map<uint32_t, shared_ptr<unordered_set<shared_ptr<temporal_bipartite_edge>>>>> temporal_left_vertex::get_edge_map() {
        return edge_map;
    }

    void temporal_left_vertex::insert_edge(uint32_t r, const shared_ptr<temporal_bipartite_edge> &e) {
        if (!edge_map->count(r)) {
            edge_map->insert({r, make_shared<unordered_set<shared_ptr<temporal_bipartite_edge>>>()});
        }
        edge_map->at(r)->insert(e);
    }

    /**
     * @details insert an edge into this vertex
     * @param r
     * @param edge
     */
    void temporal_left_vertex::insert_edge_set(uint32_t r, const shared_ptr<unordered_set<shared_ptr<temporal_bipartite_edge>>> &e) {
        if (!edge_map->count(r)) {
            edge_map->insert({r, make_shared<unordered_set<shared_ptr<temporal_bipartite_edge>>>()});
        }
        copy(e->begin(), e->end(), inserter(*edge_map->at(r), edge_map->at(r)->end()));
    }

    /**
     * @details remove an e with the given right vertex id
     * @param r
     */
    void temporal_left_vertex::remove_edge(uint32_t r, const shared_ptr<temporal_bipartite_edge>& e) {
        if(edge_map->count(r)){
            edge_map->at(r)->erase(e);
            if(edge_map->at(r)->empty()){
                edge_map->erase(r);
            }
        }
    }

    /**
     * @details remove an edge with the given right vertex id
     * @param r
     */
    void temporal_left_vertex::remove_edge_set(uint32_t r) {
        edge_map->erase(r);
    }

    void temporal_left_vertex::remove_edge_set(uint32_t r,
                                               const shared_ptr<unordered_set<shared_ptr<temporal_bipartite_edge>>> &edge_set) {
        if(edge_map->count(r)){
            for(const auto &e:*edge_set){
                edge_map->at(r)->erase(e);
            }
        }

    }
}




