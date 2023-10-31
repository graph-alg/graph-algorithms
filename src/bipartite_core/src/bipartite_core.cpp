
#include "bipartite_core/bipartite_core.h"

namespace scnu{
    bipartite_core::bipartite_core(): left_vertex_set(make_shared<unordered_set<uint32_t>>()),
                      right_vertex_set(make_shared<unordered_set<uint32_t>>()){
    }

    bipartite_core::bipartite_core(const shared_ptr<bipartite_core>& other_bipartite_core):
            bipartite_core(other_bipartite_core->get_left_vertex_set(),
                           other_bipartite_core->get_right_vertex_set())
    {

    }

    bipartite_core::bipartite_core(const shared_ptr<unordered_set<uint32_t>>& other_left_vertex_set,
                   const shared_ptr<unordered_set<uint32_t>>& other_right_vertex_set):
            bipartite_core(){
        copy(other_left_vertex_set->begin(), other_left_vertex_set->end(),
             inserter(*left_vertex_set, left_vertex_set->begin()));
        copy(other_right_vertex_set->begin(), other_right_vertex_set->end(),
             inserter(*right_vertex_set, right_vertex_set->begin()));
    }

    bool bipartite_core::count_left_vertex(uint32_t l)
    {
        return left_vertex_set->count(l);
    }


    shared_ptr<unordered_set<uint32_t>> bipartite_core::get_left_vertex_set()
    {
        return left_vertex_set;
    }

    shared_ptr<unordered_set<uint32_t>> bipartite_core::get_right_vertex_set()
    {
        return right_vertex_set;
    }

    void bipartite_core::insert_left_vertex(uint32_t l)
    {
        left_vertex_set->insert(l);
    }

    void bipartite_core::insert_right_vertex(uint32_t r)
    {
        right_vertex_set->insert(r);
    }


    void bipartite_core::remove_left_vertex(uint32_t l)
    {
        left_vertex_set->erase(l);
    }

    void bipartite_core::remove_right_vertex(uint32_t r)
    {
        right_vertex_set->erase(r);
    }


    bool bipartite_core::count_right_vertex(uint32_t r)
    {
        return right_vertex_set->count(r);
    }


    bool bipartite_core::empty()
    {
        return left_vertex_set->empty() || right_vertex_set->empty();
    }
}