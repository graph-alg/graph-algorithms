
#include "graph/temporal_vertex.h"


namespace scnu{
    temporal_vertex::temporal_vertex(uint32_t other_vertex_id) {
        this->vertex_id = other_vertex_id;
        this->neighbor_vertex_map = make_shared<unordered_map<uint32_t, shared_ptr<unordered_set<shared_ptr<temporal_edge>>>>>();
    }

    temporal_vertex::temporal_vertex(const shared_ptr<temporal_vertex> &other_vertex) {
        this->vertex_id = other_vertex->get_vertex_id();
        this->neighbor_vertex_map = make_shared<unordered_map<uint32_t, shared_ptr<unordered_set<shared_ptr<temporal_edge>>>>>();
        for (const auto&[v,v_temporal_edge_set]:*other_vertex->get_neighbor_map()) {
            auto new_temporal_edge_set = make_shared<unordered_set<shared_ptr<temporal_edge>>>();
            copy(v_temporal_edge_set->begin(),v_temporal_edge_set->end(), inserter(*new_temporal_edge_set,new_temporal_edge_set->end()));
            this->neighbor_vertex_map->insert({v,new_temporal_edge_set});
        }
    }

    temporal_vertex::temporal_vertex(uint32_t other_vertex_id,
                                     const shared_ptr<unordered_map<uint32_t, shared_ptr<unordered_set<shared_ptr<temporal_edge>>>>> &other_neighbor_vertex_map) {
        this->vertex_id = other_vertex_id;
        this->neighbor_vertex_map = other_neighbor_vertex_map;
    }


    bool temporal_vertex::exist_edge(uint32_t neighbor_vertex_id, const shared_ptr<temporal_edge> &edge) {
        if (neighbor_vertex_map->count(neighbor_vertex_id) == 0) {
            return false;
        }
        return neighbor_vertex_map->at(neighbor_vertex_id)->count(edge);
    }


    uint32_t temporal_vertex::get_vertex_id() const {
        return this->vertex_id;
    }

    shared_ptr<unordered_set<shared_ptr<temporal_edge>>> temporal_vertex::get_temporal_edge_set(uint32_t neighbor_vertex_id) {
        if (neighbor_vertex_map->count(neighbor_vertex_id)) {
            return neighbor_vertex_map->at(neighbor_vertex_id);
        }
        return shared_ptr<unordered_set<shared_ptr<temporal_edge>>>();
    }

    uint32_t temporal_vertex::get_temporal_edge_size(uint32_t neighbor_vertex_id) {
        return neighbor_vertex_map->count(neighbor_vertex_id)?neighbor_vertex_map->at(neighbor_vertex_id)->size():0;
    }

    shared_ptr<unordered_map<uint32_t, shared_ptr<unordered_set<shared_ptr<temporal_edge>>>>> temporal_vertex::get_neighbor_map() {
        return this->neighbor_vertex_map;
    }

    uint32_t temporal_vertex::get_neighbor_size(){
        return this->neighbor_vertex_map->size();
    }

    void temporal_vertex::insert_edge(uint32_t neighbor_vertex_id, const shared_ptr<temporal_edge> &edge)
    {
        if (!neighbor_vertex_map->count(neighbor_vertex_id))
        {
            neighbor_vertex_map->insert({neighbor_vertex_id, make_shared<unordered_set<shared_ptr<temporal_edge>>>()});
        }
        neighbor_vertex_map->at(neighbor_vertex_id)->insert(edge);
    }

    void temporal_vertex::insert_edge(uint32_t neighbor_vertex_id, const shared_ptr<unordered_set<shared_ptr<temporal_edge>>>& temporal_edge_set) {
        if (!neighbor_vertex_map->count(neighbor_vertex_id)) {
            auto new_temporal_edge_set = make_shared<unordered_set<shared_ptr<temporal_edge>>>();
            copy(temporal_edge_set->begin(),temporal_edge_set->end(), inserter(*new_temporal_edge_set,new_temporal_edge_set->end()));
            neighbor_vertex_map->insert({neighbor_vertex_id, new_temporal_edge_set});
        }else
        {
            for (const auto& edge:*temporal_edge_set) {
                neighbor_vertex_map->at(neighbor_vertex_id)->insert(edge);
            }
        }
    }

    void temporal_vertex::insert_neighbor_vertex(uint32_t neighbor_vertex_id) {
        neighbor_vertex_map->insert({neighbor_vertex_id, make_shared<unordered_set<shared_ptr<temporal_edge>>>()});
    }

    void temporal_vertex::insert_neighbor_vertex(uint32_t neighbor_vertex_id, const shared_ptr<unordered_set<shared_ptr<temporal_edge>>>& temporal_edge_set) {
        neighbor_vertex_map->insert({neighbor_vertex_id, temporal_edge_set});
    }

    void temporal_vertex::remove_edge(uint32_t neighbor_vertex_id, const shared_ptr<temporal_edge> &edge) {
        this->neighbor_vertex_map->at(neighbor_vertex_id)->erase(edge);
    }

    void temporal_vertex::remove_edge(uint32_t neighbor_vertex_id, const shared_ptr<unordered_set<shared_ptr<temporal_edge>>>& sub_temporal_edge_set) {
        for (const auto& edge:*sub_temporal_edge_set) {
            this->neighbor_vertex_map->at(neighbor_vertex_id)->erase(edge);
        }
    }

    void temporal_vertex::remove_neighbor_vertex(uint32_t neighbor_vertex_id) {
        this->neighbor_vertex_map->erase(neighbor_vertex_id);
    }
}