

#include "graph/weighted_bipartite_edge.h"

namespace scnu
{
    weighted_bipartite_edge::weighted_bipartite_edge(uint32_t other_left_vertex_id, uint32_t other_right_vertex_id) :
            weighted_bipartite_edge(other_left_vertex_id, other_right_vertex_id, 0){

    }

    weighted_bipartite_edge::weighted_bipartite_edge(uint32_t other_left_vertex_id, uint32_t other_right_vertex_id,
                                                     double other_weight) :
            abstract_bipartite_edge(other_left_vertex_id, other_right_vertex_id), weight(other_weight) {

    }

    double weighted_bipartite_edge::get_weight() const {
        return weight;
    }

    void weighted_bipartite_edge::update_weight(double other_weight) {
        weight = other_weight;
    }

    bool weighted_bipartite_edge::operator==(const shared_ptr<weighted_bipartite_edge> &other_edge) const {
        return get_left_vertex_id() == other_edge->get_left_vertex_id()
               && get_right_vertex_id() == other_edge->get_right_vertex_id()
               && weight == other_edge->get_weight();
    }
}