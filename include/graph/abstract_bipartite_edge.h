/******************************************************************************************************************
* @copyright  : Boost Software License
* @file       : abstract_bipartite_edge.h
* @brief      : An abstract edge structure for bipartite graph
* @version    : 1.1
* @date       : 2020/02/10
******************************************************************************************************************/

#pragma once
#include "graph/graph_utility.h"

namespace scnu {
    /**
     * @details an abstract edge class, without direction and weight
     */
    class abstract_bipartite_edge {
    public:
        abstract_bipartite_edge(uint32_t other_left_vertex_id, uint32_t other_right_vertex_id);

        explicit abstract_bipartite_edge(const shared_ptr<abstract_bipartite_edge> &e);

        [[nodiscard]] uint32_t get_left_vertex_id() const;

        [[nodiscard]] uint32_t get_right_vertex_id() const;

    private:
        uint32_t left_vertex_id;
        uint32_t right_vertex_id;
    };

    struct hash_abstract_bipartite_edge {
        size_t operator()(const shared_ptr<abstract_bipartite_edge> &e) const {
            stringstream input_stream;
            input_stream<<e->get_left_vertex_id()<< "," <<e->get_right_vertex_id();
            return hash<std::string>()(input_stream.str());
        }
    };

    struct equal_abstract_bipartite_edge {
        bool operator()(const shared_ptr<abstract_bipartite_edge> &e1, const shared_ptr<abstract_bipartite_edge> &e2) const {
            return e1->get_left_vertex_id() == e2->get_left_vertex_id()
                   && e1->get_right_vertex_id() == e2->get_right_vertex_id();
        }
    };
}



