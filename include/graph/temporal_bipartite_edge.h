/******************************************************************************************************************
* @copyright  : Boost Software License
* @file       : temporal_bipartite_edge.h
* @brief      : A simple bipartite edge with extra time information
* @version    : 1.1
* @date       : 2020/02/19
******************************************************************************************************************/

#pragma once
#include "graph/weighted_bipartite_edge.h"

namespace scnu
{
    /**
     * @details a simple edge with a comparable timestamp
     */
    class temporal_bipartite_edge : public weighted_bipartite_edge{
    public:
        temporal_bipartite_edge(uint32_t other_left_vertex_id, uint32_t other_right_vertex_id);

        temporal_bipartite_edge(uint32_t other_left_vertex_id, uint32_t other_right_vertex_id,
                                double other_weight);

        temporal_bipartite_edge(uint32_t other_left_vertex_id, uint32_t other_right_vertex_id,
                                double other_weight,
                                uint32_t other_timestamp);

        [[nodiscard]] uint32_t get_timestamp() const;

        void update_timestamp(uint32_t other_timestamp);

        bool operator<(const shared_ptr<temporal_bipartite_edge> &other_edge) const;

        bool operator==(const shared_ptr<temporal_bipartite_edge> &other_edge) const;

    private:
        uint32_t timestamp;
    };

    struct hash_temporal_bipartite_edge {
        size_t operator()(const shared_ptr<temporal_bipartite_edge> &e) const {
            stringstream input_stream;
            input_stream<<e->get_left_vertex_id()<< "," <<e->get_right_vertex_id()<<","<<e->get_weight()<< ","<<e->get_timestamp();
            return hash<std::string>()(input_stream.str());
        }
    };

    struct equal_temporal_bipartite_edge {
        bool operator()(const shared_ptr<temporal_bipartite_edge> &e1, const shared_ptr<temporal_bipartite_edge> &e2) const {
            return e1->get_left_vertex_id() == e2->get_left_vertex_id()
                   && e1->get_right_vertex_id() == e2->get_right_vertex_id()
                   && e1->get_weight() == e2->get_weight()
                   && e1->get_timestamp() == e2->get_timestamp();
        }
    };
}



