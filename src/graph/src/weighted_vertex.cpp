
#include "graph/weighted_vertex.h"

namespace scnu{
    weighted_vertex::weighted_vertex(uint32_t other_vertex_id) {
        this->vertex_id = other_vertex_id;
        this->neighbor_vertex_map = make_shared<unordered_map<uint32_t, shared_ptr<weighted_edge>>>();
    }

    weighted_vertex::weighted_vertex(const shared_ptr<weighted_vertex> &other_vertex) {
        this->vertex_id = other_vertex->get_vertex_id();
        this->neighbor_vertex_map = make_shared<unordered_map<uint32_t, shared_ptr<weighted_edge>>>();
        for (const auto&[v,e]:*other_vertex->get_neighbor_map()) {
            this->neighbor_vertex_map->insert({v,e});
        }
    }

    weighted_vertex::weighted_vertex(uint32_t other_vertex_id,
                                     const shared_ptr<unordered_map<uint32_t, shared_ptr<weighted_edge>>> &other_neighbor_vertex_map) {
        this->vertex_id = other_vertex_id;
        this->neighbor_vertex_map = other_neighbor_vertex_map;
    }

    bool weighted_vertex::exist_edge(uint32_t neighbor_vertex_id){
        return get_edge(neighbor_vertex_id) != shared_ptr<weighted_edge>();
    }

    uint32_t weighted_vertex::get_vertex_id() const {
        return this->vertex_id;
    }

    shared_ptr<weighted_edge> weighted_vertex::get_edge(uint32_t neighbor_vertex_id) {
        return neighbor_vertex_map->count(neighbor_vertex_id)? neighbor_vertex_map->at(neighbor_vertex_id): shared_ptr<weighted_edge>();
    }

    shared_ptr<unordered_map<uint32_t, shared_ptr<weighted_edge>>> weighted_vertex::get_neighbor_map() {
        return this->neighbor_vertex_map;
    }

    uint32_t weighted_vertex::get_degree(){
        return this->neighbor_vertex_map->size();
    }

    void weighted_vertex::insert_edge(uint32_t neighbor_vertex_id, const shared_ptr<weighted_edge> &edge)
    {
        if (!neighbor_vertex_map->count(neighbor_vertex_id))
        {
            neighbor_vertex_map->insert({neighbor_vertex_id, edge});
        }
    }

    void weighted_vertex::remove_edge(uint32_t neighbor_vertex_id) {
        this->neighbor_vertex_map->erase(neighbor_vertex_id);
    }
}
