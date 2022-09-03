/******************************************************************************************************************
* @copyright  : Boost Software License
* @file       : temporal_edge.h
* @brief      : A simple edge with extra time information
* @version    : 1.1
* @date       : 2020/10/15
******************************************************************************************************************/

#pragma once

#include "graph/graph_utility.h"
#include "graph/weighted_edge.h"

namespace scnu
{
    /**
     * @details a simple edge with a comparable timestamp
     */
    class temporal_edge : public weighted_edge{
    public:
        temporal_edge(uint32_t other_source_vertex_id, uint32_t other_destination_vertex_id);

        temporal_edge(uint32_t other_source_vertex_id, uint32_t other_destination_vertex_id,
                      double weight);

        temporal_edge(uint32_t other_source_vertex_id, uint32_t other_destination_vertex_id,
                      double weight,
                      uint32_t other_timestamp);

        [[nodiscard]] uint32_t get_timestamp() const;

        void set_timestamp(uint32_t other_timestamp);

    private:
        uint32_t timestamp;
    };
    /**
     * @details the hash struct for temporal edge
     */
    struct hash_temporal_edge {
        size_t operator()(const shared_ptr<temporal_edge> &e) const {
            stringstream input_stream;
            input_stream<<e->get_source_vertex_id()<< "," <<e->get_destination_vertex_id()<<"," << e->get_weight() <<","<<e->get_timestamp();
            return hash<std::string>()(input_stream.str());
        }
    };

    struct equal_temporal_edge {
        bool operator()(const shared_ptr<temporal_edge> &e1, const shared_ptr<temporal_edge> &e2) const {
            return e1->get_source_vertex_id() == e2->get_source_vertex_id()
                   && e1->get_destination_vertex_id() == e2->get_destination_vertex_id()
                   && e1->get_weight() == e2->get_weight()
                   && e1->get_timestamp() == e2->get_timestamp();
        }
    };
}



