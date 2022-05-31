/******************************************************************************************************************
* @copyright  : Boost Software License
* @file       : temporal_bipartite_edge.h
* @brief      : A simple bipartite edge with extra time information
* @version    : 1.1
* @date       : 2020/02/19
******************************************************************************************************************/

#pragma once
#include "graph/abstract_bipartite_edge.h"

namespace scnu
{
    /**
     * @details a simple edge with a comparable timestamp
     */
    class temporal_bipartite_edge : public abstract_bipartite_edge{
    public:
        temporal_bipartite_edge(uint32_t other_left_vertex_id, uint32_t other_right_vertex_id);

        temporal_bipartite_edge(uint32_t other_left_vertex_id, uint32_t other_right_vertex_id,
                      uint32_t other_timestamp);

        [[nodiscard]] uint32_t get_timestamp() const;

        void update_timestamp(uint32_t other_timestamp);

        bool operator<(const shared_ptr<temporal_bipartite_edge> &other_edge) const;

        bool operator==(const shared_ptr<temporal_bipartite_edge> &other_edge) const;

    private:
        uint32_t timestamp;
    };
}



