
#include "graph/abstract_graph.h"

namespace scnu {
    /**
     * @details construct an empty graph
     */
    abstract_graph::abstract_graph() : vertex_map(make_shared<unordered_map<uint32_t, shared_ptr<abstract_vertex>>>()) {

    }

    /**
     * @details construct a graph with the given vertex map
     * @param other_vertex_map
     */
    abstract_graph::abstract_graph(const shared_ptr<unordered_map<uint32_t, shared_ptr<abstract_vertex>>> &other_vertex_map)
            : abstract_graph() {
        for (const auto &[vertex_id, vertex]:*other_vertex_map) {
            vertex_map->insert({vertex_id, make_shared<abstract_vertex>(vertex)});
        }
    }

    /**
     * @details construct a graph with a given graph
     * @param other_graph
     */
    abstract_graph::abstract_graph(const shared_ptr<abstract_graph> &other_graph) :
            abstract_graph(other_graph->get_vertex_map()) {

    }

    /**
     * @details check this graph is empty or not
     * @return
     */
    bool abstract_graph::empty() {
        return vertex_map->empty();
    }

    /**
     * @details get an edge with the given source vertex id and the given destination vertex id
     * @param source_vertex_id
     * @param destination_vertex_id
     * @return
     */
    shared_ptr<abstract_edge>
    abstract_graph::get_edge(uint32_t source_vertex_id, uint32_t destination_vertex_id) {
        return vertex_map->count(source_vertex_id) ? vertex_map->at(source_vertex_id)->get_edge(destination_vertex_id):
        shared_ptr<abstract_edge>();
    }

    /**
     * @details get the maximal degree of vertices in this graph
     * @return
     */
    uint32_t abstract_graph::get_maximal_degree() {
        uint32_t value = 0;
        for (const auto &[vertex_id, vertex]:*vertex_map) {
            if (vertex->get_degree() > value) {
                value = vertex->get_degree();
            }
        }
        return value;
    }

    /**
     * @details get the total number of edges of this graph
     * @remarks this is undirected graph, so there is only one edge between two vertices
     * @return
     */
    uint32_t abstract_graph::get_edge_number() {
        uint32_t value = 0;
        for (const auto &[vertex_id, vertex]:*vertex_map) {
            value += vertex->get_degree();
        }
        return value / 2;
    }

    /**
     * @details get the total number of vertices of this graph
     * @return
     */
    uint32_t abstract_graph::get_vertex_number() {
        return vertex_map->size();
    }

    /**
     * @details get the edge set of this graph
     * @return
     */
    shared_ptr<unordered_set<shared_ptr<abstract_edge>>> abstract_graph::get_edge_set() {
        auto edge_set = make_shared<unordered_set<shared_ptr<abstract_edge>>>();
        for (const auto &[source_vertex_id, source_vertex] : *vertex_map) {
            for (const auto &[destination_vertex_id, e] : *source_vertex->get_edge_map()) {
                edge_set->insert(e);
            }
        }
        return edge_set;
    }

    /**
     * @details get a vertex with the given vertex id
     * @param vertex_id
     * @return
     */
    shared_ptr<abstract_vertex> abstract_graph::get_vertex(uint32_t vertex_id) {
        shared_ptr<abstract_vertex> vertex;
        if(vertex_map->count(vertex_id)){
            vertex = vertex_map->at(vertex_id);
        }
        return vertex;
    }

    /**
     * @details get the vertex map of this graph
     * @return
     */
    shared_ptr<unordered_map<uint32_t, shared_ptr<abstract_vertex>>> abstract_graph::get_vertex_map() {
        return vertex_map;
    }

    /**
     * @details insert an edge into this graph
     * @param edge
     */
    void abstract_graph::insert_edge(const shared_ptr<abstract_edge> &edge) {
        auto source_vertex_id = edge->get_source_vertex_id();
        auto destination_vertex_id = edge->get_destination_vertex_id();

        if (!get_edge(source_vertex_id, destination_vertex_id)) {
            auto source_vertex = get_vertex(source_vertex_id);
            if (!source_vertex) {
                source_vertex = make_shared<abstract_vertex>(source_vertex_id);
                vertex_map->insert({source_vertex_id, source_vertex});
            }
            source_vertex->insert_edge(destination_vertex_id, edge);

            auto destination_vertex = get_vertex(destination_vertex_id);
            if (!destination_vertex) {
                destination_vertex = make_shared<abstract_vertex>(destination_vertex_id);
                vertex_map->insert({destination_vertex_id, destination_vertex});
            }
            destination_vertex->insert_edge(source_vertex_id, edge);
        }
    }


    /**
     * @brief insert a new vertex with the vertex id
     * @param vertex_id
     */
    void abstract_graph::insert_vertex(uint32_t vertex_id) {
        if (!vertex_map->count(vertex_id)) {
            vertex_map->insert({vertex_id, make_shared<abstract_vertex>(vertex_id)});
        }
    }

    /**
     * @details remove an edge from this graph
     * @param edge
     */
    shared_ptr<unordered_set<uint32_t>> abstract_graph::remove_edge(const shared_ptr<abstract_edge> &edge) {
        auto source_vertex_id = edge->get_source_vertex_id();
        auto destination_vertex_id = edge->get_destination_vertex_id();

        auto isolated_vertex_set = make_shared<unordered_set<uint32_t>>();
        if (get_edge(source_vertex_id, destination_vertex_id)) {
            auto source_vertex = get_vertex(source_vertex_id);
            source_vertex->remove_edge(destination_vertex_id);

            if (source_vertex->get_degree() == 0) {
                remove_vertex(source_vertex_id);
                isolated_vertex_set->insert(source_vertex_id);
            }

            auto destination_vertex = get_vertex(destination_vertex_id);
            destination_vertex->remove_edge(source_vertex_id);

            if (destination_vertex->get_degree() == 0) {
                remove_vertex(destination_vertex_id);
                isolated_vertex_set->insert(destination_vertex_id);
            }
        }
        return isolated_vertex_set;
    }

    shared_ptr<unordered_set<uint32_t>> abstract_graph::remove_edge(const shared_ptr<abstract_edge> &edge,
                                                                    const shared_ptr<unordered_set<uint32_t>>& isolated_vertex_set) {
        auto source_vertex_id = edge->get_source_vertex_id();
        auto destination_vertex_id = edge->get_destination_vertex_id();

        if (get_edge(source_vertex_id, destination_vertex_id)) {
            auto source_vertex = get_vertex(source_vertex_id);
            source_vertex->remove_edge(destination_vertex_id);

            if (source_vertex->get_degree() == 0) {
                remove_vertex(source_vertex_id);
                isolated_vertex_set->insert(source_vertex_id);
            }

            auto destination_vertex = get_vertex(destination_vertex_id);
            destination_vertex->remove_edge(source_vertex_id);

            if (destination_vertex->get_degree() == 0) {
                remove_vertex(destination_vertex_id);
                isolated_vertex_set->insert(destination_vertex_id);
            }
        }
        return isolated_vertex_set;
    }

    /**
     * @details remove a vertex with the given vertex id from this graph
     * @param vertex_id
     */
    void abstract_graph::remove_vertex(uint32_t vertex_id) {
        if (vertex_map->count(vertex_id)) {
            vertex_map->erase(vertex_id);
        }
    }
}
