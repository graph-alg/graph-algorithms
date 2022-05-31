
#include "graph/abstract_bipartite_edge.h"

namespace scnu {
    /**
     * @details construct an edge with the given left vertex id and the right vertex id
     * @param other_left_vertex_id
     * @param other_right_vertex_id
     */
    abstract_bipartite_edge::abstract_bipartite_edge(uint32_t other_left_vertex_id, uint32_t other_right_vertex_id) {
        left_vertex_id = other_left_vertex_id;
        right_vertex_id = other_right_vertex_id;
    }

    /**
     * @details  construct an edge with the given edge
     * @param e
     */
    abstract_bipartite_edge::abstract_bipartite_edge(const shared_ptr<abstract_bipartite_edge> &e)
            : abstract_bipartite_edge(e->left_vertex_id, e->right_vertex_id) {

    }

    /**
     * @details get the left vertex id of this edge
     * @return
     */
    uint32_t abstract_bipartite_edge::get_left_vertex_id() const {
        return left_vertex_id;
    }

    /**
     * @details get the right vertex id of this edge
     * @return
     */
    uint32_t abstract_bipartite_edge::get_right_vertex_id() const {
        return right_vertex_id;
    }
}





