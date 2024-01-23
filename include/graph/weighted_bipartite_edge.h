/******************************************************************************************************************
* @copyright  : Boost Software License
* @file       : weighted_bipartite_edge.h
* @brief      : A bipartite edge with extra weight
* @version    : 1.1
* @date       : 2022/06/17
******************************************************************************************************************/
#pragma once
#include "graph/abstract_bipartite_edge.h"

namespace scnu
{
    class weighted_bipartite_edge : public abstract_bipartite_edge {
    public:
        weighted_bipartite_edge(uint32_t other_left_vertex_id, uint32_t other_right_vertex_id);

        weighted_bipartite_edge(uint32_t other_left_vertex_id, uint32_t other_right_vertex_id,
                                double other_weight);

        double get_weight() const;

        void update_weight(double other_weight);


        bool operator==(const shared_ptr<weighted_bipartite_edge> &other_edge) const;

    private:
        double weight;
    };

    struct hash_weighted_bipartite_edge {
        size_t operator()(const shared_ptr<weighted_bipartite_edge> &e) const {
            stringstream input_stream;
            input_stream<<e->get_left_vertex_id()<< "," <<e->get_right_vertex_id() << "," <<e->get_weight();
            return hash<std::string>()(input_stream.str());
        }
    };

    struct equal_weighted_bipartite_edge {
        bool operator()(const shared_ptr<weighted_bipartite_edge> &e1, const shared_ptr<weighted_bipartite_edge> &e2) const {
            return e1->get_left_vertex_id() == e2->get_left_vertex_id()
                   && e1->get_right_vertex_id() == e2->get_right_vertex_id()
                   && e1->get_weight() == e2->get_weight();
        }
    };
}

