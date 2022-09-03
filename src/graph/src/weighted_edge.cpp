
#include "graph/weighted_edge.h"

namespace scnu{

    weighted_edge::weighted_edge(uint32_t other_source_vertex_id, uint32_t other_destination_vertex_id) :
            weighted_edge(other_source_vertex_id, other_destination_vertex_id,0){

    }

    weighted_edge::weighted_edge(uint32_t other_source_vertex_id, uint32_t other_destination_vertex_id,
                                 uint32_t other_weight) :
            abstract_edge(other_source_vertex_id, other_destination_vertex_id),weight(other_weight) {

    }

    /**
     * @details get the timestamp of this edge
     * @return
     */
    double weighted_edge::get_weight() const {
        return weight;
    }

    /**
     * @details update the timestamp of this edge
     * @param timestamp
     */
    void weighted_edge::set_weight(uint32_t other_weight) {
        weight = other_weight;
    }
}