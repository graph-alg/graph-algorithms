
#include "graph/temporal_edge.h"

namespace scnu
{
    /**
     * @details construct a new edge with timestamp = 0
     * @remarks default case constructs an abstract edge
     * @param other_source_vertex_id
     * @param other_destination_vertex_id
     */
    temporal_edge::temporal_edge(uint32_t other_source_vertex_id, uint32_t other_destination_vertex_id) :
            temporal_edge(other_source_vertex_id, other_destination_vertex_id,0){

    }

    /**
     * @details construct a new edge with the given timestamp
     * @param other_source_vertex_id
     * @param other_destination_vertex_id
     * @param other_timestamp
     */
    temporal_edge::temporal_edge(uint32_t other_source_vertex_id, uint32_t other_destination_vertex_id,
                                 uint32_t other_timestamp) :
            abstract_edge(other_source_vertex_id, other_destination_vertex_id),timestamp(other_timestamp) {

    }

    /**
     * @details get the timestamp of this edge
     * @return
     */
    uint32_t temporal_edge::get_timestamp() const {
        return timestamp;
    }

    /**
     * @details update the timestamp of this edge
     * @param timestamp
     */
    void temporal_edge::update_timestamp(uint32_t other_timestamp) {
        timestamp = other_timestamp;
    }

    /**
     * @details overload operator < to implement compare between two edges
     * @param other_edge
     * @return
     */
    bool temporal_edge::operator<(const shared_ptr<temporal_edge> &other_edge) const {
        return timestamp < other_edge->timestamp;
    }

    /**
     * @details overload operator == to remove repeated edges
     * @param other_edge
     * @return
     */
    bool temporal_edge::operator==(const shared_ptr<temporal_edge> &other_edge) const {
        return get_source_vertex_id() == other_edge->get_source_vertex_id()
               && get_destination_vertex_id() == other_edge->get_destination_vertex_id()
        && timestamp == other_edge->timestamp;
    }
}