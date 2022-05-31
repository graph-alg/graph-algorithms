
#include "graph/abstract_bipartite_graph.h"

namespace scnu {
    /**
     * @details construct an empty bipartite graph
     */
    abstract_bipartite_graph::abstract_bipartite_graph() :
            left_vertex_map(make_shared<unordered_map<uint32_t, shared_ptr<abstract_left_vertex>>>()),
            right_vertex_map(make_shared<unordered_map<uint32_t, shared_ptr<abstract_right_vertex>>>()) {
    }

    /**
     * @details construct a bipartite graph with given left and right vertex maps
     * @param left_vertex_map
     * @param right_vertex_map
     */
    abstract_bipartite_graph::abstract_bipartite_graph(
            const shared_ptr<unordered_map<uint32_t, shared_ptr<abstract_left_vertex>>> &other_left_vertex_map,
            const shared_ptr<unordered_map<uint32_t, shared_ptr<abstract_right_vertex>>> &other_right_vertex_map)
            : abstract_bipartite_graph() {
        for (const auto &[left_vertex_id, left_vertex]:*other_left_vertex_map) {
            left_vertex_map->insert({left_vertex_id, make_shared<abstract_left_vertex>(left_vertex)});
        }

        for (const auto &[right_vertex_id, right_vertex]:*other_right_vertex_map) {
            right_vertex_map->insert({right_vertex_id, make_shared<abstract_right_vertex>(right_vertex)});
        }
    }


    /**
     * @details construct a bipartite graph with a given graph
     * @param graph
     */
    abstract_bipartite_graph::abstract_bipartite_graph(const shared_ptr<abstract_bipartite_graph> &graph)
            : abstract_bipartite_graph(graph->get_left_vertex_map(), graph->get_right_vertex_map()) {

    }

    /**
     * @details compute the density of this bipartite graph (density = |E|/(|L||R|))
     * @return
     */
    double abstract_bipartite_graph::get_density() {
        return double(get_edge_number()) / double(get_left_vertex_number() * get_right_vertex_number());
    }

    /**
     * @details get the edge with the given left vertex id and the right vertex id
     * @param left_vertex_id
     * @param right_vertex_id
     * @return
     */
    shared_ptr<abstract_bipartite_edge> abstract_bipartite_graph::get_edge(uint32_t left_vertex_id,
                                                                           uint32_t right_vertex_id) {
        auto left_vertex = get_left_vertex(left_vertex_id);
        return left_vertex ? left_vertex->get_edge(right_vertex_id) : shared_ptr<abstract_bipartite_edge>();
    }

    /**
     * @details get the  bipartite edge set of this graph
     * @return
     */
    shared_ptr<unordered_set<shared_ptr<abstract_bipartite_edge>>> abstract_bipartite_graph::get_edge_set() {
        auto edge_set = make_shared<unordered_set<shared_ptr<abstract_bipartite_edge>>>();
        for (const auto &[left_vertex_id, left_vertex]:*left_vertex_map) {
            for (const auto &[right_vertex_id, e]:*left_vertex->get_edge_map()) {
                edge_set->insert(e);
            }
        }
        return edge_set;
    }

    /**
     * @details get a left vertex with the given left vertex id
     * @param left_vertex_id
     * @return
     */
    shared_ptr<abstract_left_vertex> abstract_bipartite_graph::get_left_vertex(uint32_t left_vertex_id) {
        return left_vertex_map->count(left_vertex_id) ? left_vertex_map->at(left_vertex_id)
                                                      : shared_ptr<abstract_left_vertex>();
    }

    /**
      * @details get the left vertex map of this graph
      * @return
      */
    shared_ptr<unordered_map<uint32_t, shared_ptr<abstract_left_vertex>>>
    abstract_bipartite_graph::get_left_vertex_map() {
        return left_vertex_map;
    }

    /**
     * @details get the maximal degree of left vertices
     * @return
     */
    uint32_t abstract_bipartite_graph::get_maximal_left_vertex_degree() {
        uint32_t max_degree = 0;
        for (const auto &[left_vertex_id, left_vertex]:*left_vertex_map) {
            if (left_vertex->get_degree() > max_degree) {
                max_degree = left_vertex->get_degree();
            }
        }
        return max_degree;
    }

    /**
     * @details get the maximal degree of right vertices
     * @return
     */
    uint32_t abstract_bipartite_graph::get_maximal_right_vertex_degree() {
        uint32_t max_degree = 0;
        for (const auto &[right_vertex_id, right_vertex]:*right_vertex_map) {
            if (right_vertex->get_degree() > max_degree) {
                max_degree = right_vertex->get_degree();
            }
        }
        return max_degree;
    }

    /**
     * @details get the left vertex with the given right vertex id
     * @param right_vertex_id
     * @return
     */
    shared_ptr<abstract_right_vertex> abstract_bipartite_graph::get_right_vertex(uint32_t right_vertex_id) {
        return right_vertex_map->count(right_vertex_id) ? right_vertex_map->at(right_vertex_id)
                                                        : shared_ptr<abstract_right_vertex>();
    }


    /**
     * @details get the number of edges of this bipartite graph
     * @return
     */
    uint32_t abstract_bipartite_graph::get_edge_number() {
        uint32_t edge_number = 0;
        for (const auto &[left_vertex_id, left_vertex]:*left_vertex_map) {
            edge_number += left_vertex->get_degree();
        }
        return edge_number;
    }

    /**
     * @details get the number of left vertices
     * @return
     */
    uint32_t abstract_bipartite_graph::get_left_vertex_number() {
        return left_vertex_map->size();
    }

    /**
     * @details get the number of right vertices of this bipartite graph
     * @return
     */
    uint32_t abstract_bipartite_graph::get_right_vertex_number() {
        return right_vertex_map->size();
    }

    /**
     * @details get the right vertex map of this bipartite graph
     * @return
     */
    shared_ptr<unordered_map<uint32_t, shared_ptr<abstract_right_vertex>>>
    abstract_bipartite_graph::get_right_vertex_map() {
        return right_vertex_map;
    }

    /**
     * @details insert an edge into this bipartite graph
     * @param e
     */
    void abstract_bipartite_graph::insert_edge(const shared_ptr<abstract_bipartite_edge> &e) {
        auto left_vertex_id = e->get_left_vertex_id();
        auto right_vertex_id = e->get_right_vertex_id();

        auto left_vertex = get_left_vertex(left_vertex_id);
        if (!left_vertex) {
            insert_left_vertex(left_vertex_id);
            left_vertex = get_left_vertex(left_vertex_id);
        }
        left_vertex->insert_edge(right_vertex_id, e);

        auto right_vertex = get_right_vertex(right_vertex_id);
        if (!right_vertex) {
            insert_right_vertex(right_vertex_id);
            right_vertex = get_right_vertex(right_vertex_id);
        }
        right_vertex->insert_edge(e->get_left_vertex_id(), e);
    }

    /**
     * @details insert a new left vertex in this bipartite graph with the given left vertex id
     * @param left_vertex_id
     */
    void abstract_bipartite_graph::insert_left_vertex(uint32_t left_vertex_id) {
        left_vertex_map->insert({left_vertex_id, make_shared<abstract_left_vertex>(left_vertex_id)});
    }

    /**
     * @details insert a new right vertex in this bipartite graph with the given right vertex Id
     * @param right_vertex_id
     */
    void abstract_bipartite_graph::insert_right_vertex(uint32_t right_vertex_id) {
        right_vertex_map->insert({right_vertex_id, make_shared<abstract_right_vertex>(right_vertex_id)});
    }

    /**
     * @details check whether this bipartite graph is empty or not
     * @return
     */
    bool abstract_bipartite_graph::empty() {
        return get_left_vertex_number() == 0 || get_right_vertex_number() == 0;
    }

    /**
      * @details remove the edge from this bipartite graph
      * @param edge
      */
    void abstract_bipartite_graph::remove_edge(const shared_ptr<abstract_bipartite_edge> &edge) {
        auto left_vertex_id = edge->get_left_vertex_id();
        auto right_vertex_id = edge->get_right_vertex_id();

        if (get_edge(left_vertex_id, right_vertex_id)) {
            auto left_vertex = get_left_vertex(left_vertex_id);
            left_vertex->remove_edge(right_vertex_id);

            if (left_vertex->get_degree() == 0) {
                remove_left_vertex(left_vertex_id);
            }

            auto right_vertex = get_right_vertex(right_vertex_id);
            right_vertex->remove_edge(left_vertex_id);
            if (right_vertex->get_degree() == 0) {
                remove_right_vertex(right_vertex_id);
            }
        }
    }

    /**
     * @details remove the edge from this bipartite graph
     * @param edge
     * @param isolated_vertex_set
    */
    void abstract_bipartite_graph::remove_edge(const shared_ptr<abstract_bipartite_edge> &edge,
                                               const shared_ptr<unordered_set<uint32_t>>& isolated_left_vertex_set,
                                               const shared_ptr<unordered_set<uint32_t>>& isolated_right_vertex_set) {
        auto left_vertex_id = edge->get_left_vertex_id();
        auto right_vertex_id = edge->get_right_vertex_id();

        if (get_edge(left_vertex_id, right_vertex_id)) {
            auto left_vertex = get_left_vertex(left_vertex_id);
            left_vertex->remove_edge(right_vertex_id);

            if (left_vertex->get_degree() == 0) {
                remove_left_vertex(left_vertex_id);
                isolated_left_vertex_set->insert(left_vertex_id);
            }

            auto right_vertex = get_right_vertex(right_vertex_id);
            right_vertex->remove_edge(left_vertex_id);
            if (right_vertex->get_degree() == 0) {
                remove_right_vertex(right_vertex_id);
                isolated_right_vertex_set->insert(right_vertex_id);
            }
        }
    }

    /**
     * @details remove the left vertex from this bipartite graph with the given left vertex id
     * @param left_vertex_id
     */
    void abstract_bipartite_graph::remove_left_vertex(uint32_t left_vertex_id) {
        left_vertex_map->erase(left_vertex_id);
    }

    /**
     * @details remove the right vertex from this bipartite graph with the given right vertex id
     * @param right_vertex_id
     */
    void abstract_bipartite_graph::remove_right_vertex(uint32_t right_vertex_id) {
        right_vertex_map->erase(right_vertex_id);
    }
}