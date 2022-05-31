
#include "graph/temporal_bipartite_edge.h"

namespace scnu
{
    /**
     * @details construct a new edge with timestamp = 0
     * @remarks default case constructs an abstract edge
     * @param other_left_vertex_id
     * @param other_right_vertex_id
     */
    temporal_bipartite_edge::temporal_bipartite_edge(uint32_t other_left_vertex_id, uint32_t other_right_vertex_id) :
            temporal_bipartite_edge(other_left_vertex_id, other_right_vertex_id, 0){

    }

    /**
     * @details construct a new edge with the given timestamp
     * @param other_left_vertex_id
     * @param other_right_vertex_id
     * @param other_timestamp
     */
    temporal_bipartite_edge::temporal_bipartite_edge(uint32_t other_left_vertex_id, uint32_t other_right_vertex_id,
                                 uint32_t other_timestamp) :
            abstract_bipartite_edge(other_left_vertex_id, other_right_vertex_id), timestamp(other_timestamp) {

    }

    /**
     * @details get the timestamp of this edge
     * @return
     */
    uint32_t temporal_bipartite_edge::get_timestamp() const {
        return timestamp;
    }


    /**
     * @details update the timestamp of this edge
     * @param timestamp
     */
    void temporal_bipartite_edge::update_timestamp(uint32_t other_timestamp) {
        timestamp = other_timestamp;
    }

    /**
     * @details overload operator < to implement compare between two edges
     * @param other_edge
     * @return
     */
    bool temporal_bipartite_edge::operator<(const shared_ptr<temporal_bipartite_edge> &other_edge) const {
        return timestamp < other_edge->timestamp;
    }

    /**
     * @details overload operator == to remove repeated edges
     * @param other_edge
     * @return
     */
    bool temporal_bipartite_edge::operator==(const shared_ptr<temporal_bipartite_edge> &other_edge) const {
        return get_left_vertex_id() == other_edge->get_left_vertex_id()
               && get_right_vertex_id() == other_edge->get_right_vertex_id()
                && timestamp == other_edge->timestamp;
    }
}