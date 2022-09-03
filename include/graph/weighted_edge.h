/******************************************************************************************************************
* @copyright  : Boost Software License
* @file       : weighted_edge.h
* @brief      : A simple edge with extra weight
* @version    : 1.1
* @date       : 2022/6/15
******************************************************************************************************************/

#pragma once

#include "graph/graph_utility.h"
#include "graph/abstract_edge.h"

#pragma once
namespace scnu{
    /**
     * @details a simple edge with a comparable timestamp
     */
    class weighted_edge : public abstract_edge{
    public:
        weighted_edge(uint32_t other_left_vertex_id, uint32_t other_right_vertex_id);

        weighted_edge(uint32_t other_left_vertex_id, uint32_t other_right_vertex_id,
                      uint32_t other_weight);

        double get_weight() const;

        void set_weight(uint32_t other_weight);

    private:
        double weight;
    };

    /**
     * @details the hash struct for weighted edge
     */
    struct hash_weighted_edge {
        size_t operator()(const shared_ptr<weighted_edge> &e) const {
            stringstream input_stream;
            input_stream<<e->get_source_vertex_id()<< "," <<e->get_destination_vertex_id()<<","<<e->get_weight();
            return hash<std::string>()(input_stream.str());
        }
    };

    struct equal_weighted_edge {
        bool operator()(const shared_ptr<weighted_edge> &e1, const shared_ptr<weighted_edge> &e2) const {
            return e1->get_source_vertex_id() == e2->get_source_vertex_id()
                   && e1->get_destination_vertex_id() == e2->get_destination_vertex_id()
                   && e1->get_weight() == e2->get_weight();
        }
    };
}

