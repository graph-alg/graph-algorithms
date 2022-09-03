/******************************************************************************************************************
* @copyright  : Boost Software License
* @file       : weighted_left_vertex.h
* @brief      : a right vertex with weighted edges
* @version    : 1.0
* @date       : 2022/06/17
******************************************************************************************************************/

#include "graph/temporal_bipartite_edge.h"

namespace scnu
{
    /**
     * @brief An abstract right vertex class for a abstract bipartite graph
     */
    class temporal_right_vertex
    {
    public:
        explicit temporal_right_vertex(uint32_t other_r);

        explicit temporal_right_vertex(const shared_ptr<temporal_right_vertex>& other_r_vertex);

        virtual ~temporal_right_vertex() = default;

        shared_ptr<unordered_set<shared_ptr<temporal_bipartite_edge>>> get_edge_set(uint32_t other_r);

        shared_ptr<unordered_map<uint32_t,shared_ptr<unordered_set<shared_ptr<temporal_bipartite_edge>>>>> get_edge_map();

        uint32_t get_neighbor_size();

        [[nodiscard]] uint32_t get_right_vertex_id() const;

        void insert_edge(uint32_t l, const shared_ptr<temporal_bipartite_edge> &e);

        void insert_edge_set(uint32_t l, const shared_ptr<unordered_set<shared_ptr<temporal_bipartite_edge>>>& edge_set);

        void remove_edge(uint32_t l, const shared_ptr<temporal_bipartite_edge> &e);

        void remove_edge_set(uint32_t l);

        void remove_edge_set(uint32_t l, const shared_ptr<unordered_set<shared_ptr<temporal_bipartite_edge>>>& edge_set);

    private:
        uint32_t right_vertex_id;
        shared_ptr<unordered_map<uint32_t,shared_ptr<unordered_set<shared_ptr<temporal_bipartite_edge>>>>> edge_map;
    };
}

