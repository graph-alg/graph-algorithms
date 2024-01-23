
#include "wing/bipartite_edge_index.h"

namespace scnu{
    bipartite_edge_index::bipartite_edge_index(const shared_ptr<abstract_bipartite_edge>& other_edge):edge(other_edge),left_vertex_id_set(
            make_shared<unordered_set<uint32_t>>()),right_vertex_id_set(make_shared<unordered_set<uint32_t>>())
    {

    }

    shared_ptr<abstract_bipartite_edge> bipartite_edge_index::get_edge()
    {
        return edge;
    }



    shared_ptr<unordered_set<uint32_t>> bipartite_edge_index::get_left_vertex_id_set(){
        return left_vertex_id_set;
    }

    shared_ptr<unordered_set<uint32_t>> bipartite_edge_index::get_right_vertex_id_set(){
        return right_vertex_id_set;
    }



    void bipartite_edge_index::insert_dual_edge(const shared_ptr<abstract_bipartite_edge>& dual_edge){
        auto l = dual_edge->get_left_vertex_id();
        auto r = dual_edge->get_right_vertex_id();
        left_vertex_id_set->insert(l);
        right_vertex_id_set->insert(r);
    }

    void bipartite_edge_index::remove_dual_edge(const shared_ptr<abstract_bipartite_edge>& dual_edge){
        auto l = dual_edge->get_left_vertex_id();
        auto r = dual_edge->get_right_vertex_id();
        left_vertex_id_set->erase(l);
        right_vertex_id_set->erase(r);
    }
}