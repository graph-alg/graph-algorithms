
#include "graph/abstract_vertex.h"

namespace scnu {
    /**
     * @details construct a vertex with the given vertex id
     * @param other_vertex_id
     */
    abstract_vertex::abstract_vertex(uint32_t other_vertex_id) :
            vertex_id(other_vertex_id), edge_map(make_shared<unordered_map<uint32_t, shared_ptr<abstract_edge>>>()) {

    }

    /**
     * @details construct a vertex with the given vertex
     * @param other_vertex
     */
    abstract_vertex::abstract_vertex(const shared_ptr<abstract_vertex> &other_vertex)
            : abstract_vertex(other_vertex->get_vertex_id()) {
        for (const auto &edge_pair:*other_vertex->get_edge_map()) {
            edge_map->insert({edge_pair.first, edge_pair.second});
        }
    }

    /**
     * @details get the degree of this vertex
     * @return
     */
    uint32_t abstract_vertex::get_degree() {
        return edge_map->size();
    }

    /**
     * @details get the degree of this vertex
     * @return
     */
    uint32_t abstract_vertex::get_degree(const shared_ptr<unordered_set<shared_ptr<abstract_edge>>>& removed_edge_set) {
        auto degree = 0;
        for(const auto &[v,e]:*get_edge_map())
        {
            if(removed_edge_set->count(e))
            {
                continue;
            }
            ++degree;
        }
        return degree;
    }

    /**
     * @details get the edge map of this vertex
     * @return
     */
    shared_ptr<unordered_map<uint32_t, shared_ptr<abstract_edge>>> abstract_vertex::get_edge_map() {
        return edge_map;
    }

    /**
     * @details get an edge with the given adjacent vertex id
     * @param adjacent_vertex_id
     * @return
     */
    shared_ptr<abstract_edge> abstract_vertex::get_edge(uint32_t adjacent_vertex_id) {
        return edge_map->count(adjacent_vertex_id) ? edge_map->at(adjacent_vertex_id): shared_ptr<abstract_edge>();
    }

    /**
     * @details get the vertex id of this vertex
     * @return
     */
    uint32_t abstract_vertex::get_vertex_id() const {
        return vertex_id;
    }

    /**
     * @details insert an adjacent adjacent_edge
     * @param adjacent_vertex_id
     * @param adjacent_edge
     */
    void abstract_vertex::insert_edge(uint32_t adjacent_vertex_id,
                                      const shared_ptr<abstract_edge> &adjacent_edge) {
        if (!edge_map->count(adjacent_vertex_id)) {
            edge_map->insert({adjacent_vertex_id, adjacent_edge});
        }
    }

    /**
     * @details remove an adjacent edge
     * @param adjacent_vertex_id
     */
    void abstract_vertex::remove_edge(uint32_t adjacent_vertex_id) {
        if (edge_map->count(adjacent_vertex_id)) {
            edge_map->erase(adjacent_vertex_id);
        }
    }
}
