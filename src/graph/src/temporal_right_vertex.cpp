
#include "graph/temporal_right_vertex.h"

namespace scnu {
    /**
     * @details construct a new right vertex with the given right vertex id
     * @param right_vertex_id
     */
    temporal_right_vertex::temporal_right_vertex(uint32_t other_r) : right_vertex_id(
            other_r), edge_map(
            make_shared<unordered_map<uint32_t, shared_ptr<unordered_set<shared_ptr<temporal_bipartite_edge>>>>>()) {

    }

    /**
     * @details construct a new right vertex with the given right vertex
     * @param other_r_vertex
     */
    temporal_right_vertex::temporal_right_vertex(const shared_ptr<temporal_right_vertex> &other_r_vertex) :
            temporal_right_vertex(other_r_vertex->right_vertex_id) {
        for(const auto &[l, l_set]:*other_r_vertex->get_edge_map()){
            edge_map->insert({l, make_shared<unordered_set<shared_ptr<temporal_bipartite_edge>>>()});
            copy(l_set->begin(), l_set->end(), inserter(*edge_map->at(l),edge_map->at(l)->end()));
        }
    }

    /**
     * @details get the adjacent edge with the given left vertex id
     * @param right_vertex_id
     * @return
     */
    shared_ptr<unordered_set<shared_ptr<temporal_bipartite_edge>>> temporal_right_vertex::get_edge_set(uint32_t other_r) {
        if (edge_map->count(other_r)) {
            return edge_map->at(other_r);
        }
        return {};
    }

    /**
     * @details get the adjacent edge map of this right vertex
     * @return
     */
    shared_ptr<unordered_map<uint32_t, shared_ptr<unordered_set<shared_ptr<temporal_bipartite_edge>>>>> temporal_right_vertex::get_edge_map() {
        return edge_map;
    }


    /**
     * @details get the degree of this right vertex
     * @return
     */
    uint32_t temporal_right_vertex::get_neighbor_size() {
        return edge_map->size();
    }

    /**
     * @details get the vertex id of this vertex
     * @return
     */
    uint32_t temporal_right_vertex::get_right_vertex_id() const {
        return right_vertex_id;
    }

    /**
     * @details insert an adjacent e
     * @param l
     * @param e
     */
    void temporal_right_vertex::insert_edge(uint32_t l,
                                            const shared_ptr<temporal_bipartite_edge> &e) {
        if (!edge_map->count(l)) {
            edge_map->insert({l, make_shared<unordered_set<shared_ptr<temporal_bipartite_edge>>>()});
        }
        edge_map->at(l)->insert(e);
    }

    void temporal_right_vertex::insert_edge_set(uint32_t l,
                                            const shared_ptr<unordered_set<shared_ptr<temporal_bipartite_edge>>> &edge_set) {
        if (!edge_map->count(l)) {
            edge_map->insert({l, make_shared<unordered_set<shared_ptr<temporal_bipartite_edge>>>()});
        }
        copy(edge_set->begin(),edge_set->end(), inserter(*edge_map->at(l), edge_map->at(l)->end()));
    }

    /**
     * @details remove an adjacent edge
     * @param leftVertexId
     */
    void temporal_right_vertex::remove_edge(uint32_t l,const shared_ptr<temporal_bipartite_edge> &edge) {
        if(edge_map->count(l)){
            edge_map->at(l)->erase(edge);
            if(edge_map->at(l)->empty()){
                edge_map->erase(l);
            }
        }
    }

    void temporal_right_vertex::remove_edge_set(uint32_t l) {
        edge_map->erase(l);
    }

    void temporal_right_vertex::remove_edge_set(uint32_t l,  const shared_ptr<unordered_set<shared_ptr<temporal_bipartite_edge>>> &edge_set) {
        if(edge_map->count(l)){
            for(const auto&e:*edge_set){
                edge_map->at(l)->erase(e);
            }
        }
    }
}

