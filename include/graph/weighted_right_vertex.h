/******************************************************************************************************************
* @copyright  : Boost Software License
* @file       : weighted_right_vertex.h
* @brief      : a right vertex with weighted edges
* @version    : 1.0
* @date       : 2022/06/17
******************************************************************************************************************/

#include "graph/weighted_bipartite_edge.h"

namespace scnu
{
    /**
     * @brief An abstract right vertex class for a abstract bipartite graph
     */
    class weighted_right_vertex
    {
    public:
        explicit weighted_right_vertex(uint32_t other_right_vertex_id);

        explicit weighted_right_vertex(const shared_ptr<weighted_right_vertex>& other_right_vertex);

        virtual ~weighted_right_vertex() = default;

        shared_ptr<weighted_bipartite_edge> get_edge(uint32_t other_right_vertex_id);

        shared_ptr<unordered_map<uint32_t,shared_ptr<weighted_bipartite_edge>>> get_edge_map();

        uint32_t get_degree();

        [[nodiscard]] uint32_t get_right_vertex_id() const;

        void insert_edge(uint32_t other_left_vertex_id, const shared_ptr<weighted_bipartite_edge>& edge);

        void remove_edge(uint32_t other_left_vertex_id);

    private:
        uint32_t right_vertex_id;
        shared_ptr<unordered_map<uint32_t,shared_ptr<weighted_bipartite_edge>>> edge_map;
    };
}



