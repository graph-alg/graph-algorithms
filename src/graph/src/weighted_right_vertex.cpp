
#include "graph/weighted_right_vertex.h"

namespace scnu {
    /**
     * @details construct a new right vertex with the given right vertex id
     * @param right_vertex_id
     */
    weighted_right_vertex::weighted_right_vertex(uint32_t other_right_vertex_id) : right_vertex_id(
            other_right_vertex_id), edge_map(
            make_shared<unordered_map<uint32_t, shared_ptr<weighted_bipartite_edge>>>()) {

    }

    /**
     * @details construct a new right vertex with the given right vertex
     * @param other_right_vertex
     */
    weighted_right_vertex::weighted_right_vertex(const shared_ptr<weighted_right_vertex> &other_right_vertex) :
            weighted_right_vertex(other_right_vertex->right_vertex_id) {
        auto other_edge_map = other_right_vertex->get_edge_map();
        copy(other_edge_map->begin(), other_edge_map->end(), inserter(*edge_map,edge_map->end()));
    }

    /**
     * @details get the adjacent edge with the given left vertex id
     * @param right_vertex_id
     * @return
     */
    shared_ptr<weighted_bipartite_edge> weighted_right_vertex::get_edge(uint32_t other_right_vertex_id) {
        if (edge_map->count(other_right_vertex_id)) {
            return edge_map->at(other_right_vertex_id);
        }
        return shared_ptr<weighted_bipartite_edge>();
    }

    /**
     * @details get the adjacent edge map of this right vertex
     * @return
     */
    shared_ptr<unordered_map<uint32_t, shared_ptr<weighted_bipartite_edge>>> weighted_right_vertex::get_edge_map() {
        return edge_map;
    }


    /**
     * @details get the degree of this right vertex
     * @return
     */
    uint32_t weighted_right_vertex::get_degree() {
        return edge_map->size();
    }

    /**
     * @details get the vertex id of this vertex
     * @return
     */
    uint32_t weighted_right_vertex::get_right_vertex_id() const {
        return right_vertex_id;
    }

    /**
     * @details insert an adjacent edge
     * @param other_left_vertex_id
     * @param edge
     */
    void weighted_right_vertex::insert_edge(uint32_t other_left_vertex_id,
                                            const shared_ptr<weighted_bipartite_edge> &edge) {
        if (!edge_map->count(other_left_vertex_id)) {
            edge_map->insert({other_left_vertex_id, edge});
        }
    }

    /**
     * @details remove an adjacent edge
     * @param leftVertexId
     */
    void weighted_right_vertex::remove_edge(uint32_t other_left_vertex_id) {
        edge_map->erase(other_left_vertex_id);
    }
}

