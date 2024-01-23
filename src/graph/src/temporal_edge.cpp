
#include "graph/temporal_edge.h"

namespace scnu
{
    temporal_edge::temporal_edge(uint32_t other_source_vertex_id, uint32_t other_destination_vertex_id) :
            temporal_edge(other_source_vertex_id, other_destination_vertex_id,0,0){

    }

    temporal_edge::temporal_edge(uint32_t other_source_vertex_id, uint32_t other_destination_vertex_id,
                                 double other_weight) :
            temporal_edge(other_source_vertex_id, other_destination_vertex_id, other_weight, 0) {

    }


    temporal_edge::temporal_edge(uint32_t other_source_vertex_id, uint32_t other_destination_vertex_id,
                                 double other_weight,
                                 uint32_t other_timestamp) :
            weighted_edge(other_source_vertex_id, other_destination_vertex_id, other_weight),timestamp(other_timestamp){

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
    void temporal_edge::set_timestamp(uint32_t other_timestamp) {
        timestamp = other_timestamp;
    }
}