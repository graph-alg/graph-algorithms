
#include "graph/abstract_left_vertex.h"

namespace scnu {
    /**
     * @details construct a vertex with the given left vertex id
     * @param other_left_vertex_id
     */
    abstract_left_vertex::abstract_left_vertex(uint32_t other_left_vertex_id) :
            left_vertex_id(other_left_vertex_id),
            edge_map(make_shared<unordered_map<uint32_t, shared_ptr<abstract_bipartite_edge>>>()) {

    }

    /**
     * @details construct a vertex with the given left vertex
     * @param left_vertex
     */
    abstract_left_vertex::abstract_left_vertex(const shared_ptr<abstract_left_vertex> &left_vertex) :
            abstract_left_vertex(left_vertex->left_vertex_id) {
        for (const auto&[right_vertex_id, e]:*left_vertex->edge_map) {
            edge_map->insert({right_vertex_id, e});
        }
    }

    /**
     * @details get the id of this left vertex
     * @return
     */
    uint32_t abstract_left_vertex::get_left_vertex_id() const {
        return left_vertex_id;
    }

    /**
     * @details get the degree of this left vertex
     * @return
     */
    uint32_t abstract_left_vertex::get_degree() {
        return edge_map->size();
    }

    /**
     * @details get the adjacent edge with the given right vertex id
     * @param right_vertex_id
     * @return
     */
    shared_ptr<abstract_bipartite_edge> abstract_left_vertex::get_edge(uint32_t right_vertex_id) {
        if (edge_map->count(right_vertex_id)) {
            return edge_map->at(right_vertex_id);
        }
        return shared_ptr<abstract_bipartite_edge>();
    }

    /**
     * @details get the adjacent edge map
     * @return
     */
    shared_ptr<unordered_map<uint32_t, shared_ptr<abstract_bipartite_edge>>> abstract_left_vertex::get_edge_map() {
        return edge_map;
    }

    /**
     * @details insert an edge into this vertex
     * @param right_vertex_id
     * @param edge
     */
    void abstract_left_vertex::insert_edge(uint32_t right_vertex_id, const shared_ptr<abstract_bipartite_edge> &edge) {
        if (!edge_map->count(right_vertex_id)) {
            edge_map->insert({right_vertex_id, edge});
        }
    }

    /**
     * @details remove an edge with the given right vertex id
     * @param right_vertex_id
     */
    void abstract_left_vertex::remove_edge(uint32_t right_vertex_id) {
        edge_map->erase(right_vertex_id);
    }
}



