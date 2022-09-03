
#include "graph/temporal_graph.h"


namespace scnu{
    temporal_graph::temporal_graph() {
        vertex_map = make_shared<unordered_map<uint32_t,shared_ptr<temporal_vertex>>>();
    }

    temporal_graph::temporal_graph(const shared_ptr<unordered_map<uint32_t , shared_ptr<temporal_vertex>>> &other_vertex_map) {
        this->vertex_map = other_vertex_map;
    }

    temporal_graph::temporal_graph(const shared_ptr<temporal_graph> &graph):temporal_graph() {
        for (const auto&[v,v_vertex]: *graph->get_vertex_map()) {
            vertex_map->insert({v, make_shared<temporal_vertex>(v_vertex)});
        }
    }

    temporal_graph::~temporal_graph() = default;

    bool temporal_graph::empty(){
        return vertex_map->empty();
    }

    bool temporal_graph::exist_edge(const shared_ptr<temporal_edge>& edge) {
        auto source_vertex = get_vertex(edge->get_source_vertex_id());
        if (!source_vertex) {
            return false;
        }
        return source_vertex->exist_edge(edge->get_destination_vertex_id(), edge);
    }

    bool temporal_graph::exist_edge(uint32_t source_vertex_id, uint32_t destination_vertex_id) {
        auto source_vertex = get_vertex(source_vertex_id);
        if (!source_vertex) {
            return false;
        }
        return source_vertex->get_neighbor_map()->count(destination_vertex_id);
    }

    bool temporal_graph::exist_vertex(uint32_t vertex_id) {
        return vertex_map->count(vertex_id);
    }

    uint32_t temporal_graph::get_edge_size() {
        uint32_t edge_size = 0;
        for (const auto&[u,u_vertex]:*vertex_map) {
            for (const auto& [v, edge_set]:*u_vertex->get_neighbor_map()) {
                edge_size += edge_set->size();
            }
        }
        return edge_size / 2;
    }

    double temporal_graph::get_average_edge_size(){
        double edge_size = 0;
        uint32_t vertex_pair = 0;
        for (const auto&[u,u_vertex]:*vertex_map) {
            for (const auto& [v, edge_set]:*u_vertex->get_neighbor_map()) {
                vertex_pair ++;
                edge_size += edge_set->size();
            }
        }
        return edge_size / vertex_pair;
    }

    uint32_t temporal_graph::get_maximal_neighbor_vertex_size(){
        uint32_t max_number = 0;
        for (const auto&[v,v_vertex]:*vertex_map) {
            if (v_vertex->get_neighbor_map()->size() > max_number)
                max_number = v_vertex->get_neighbor_map()->size();
        }
        return max_number;
    }

    double temporal_graph::get_average_neighbor_vertex_size(){
        double max_number = 0;
        for (const auto&[v,v_vertex]:*vertex_map) {
            max_number += v_vertex->get_neighbor_map()->size();
        }
        return max_number/vertex_map->size();
    }

    uint32_t temporal_graph::get_maximal_parallel_edge_size(){
        uint32_t max_number = 0;
        for (const auto&[u,u_vertex]:*vertex_map) {
            for(const auto &[v,edge_set]:*u_vertex->get_neighbor_map()){
                if (edge_set->size() > max_number)
                    max_number = edge_set->size();
            }
        }
        return max_number;
    }

    shared_ptr<temporal_vertex> temporal_graph::get_vertex(uint32_t vertex_id) {
        return vertex_map->count(vertex_id) ? vertex_map->at(vertex_id) : nullptr;
    }

    shared_ptr<unordered_set<shared_ptr<temporal_edge>>> temporal_graph::get_edge_set() {
        auto edge_set = make_shared<unordered_set<shared_ptr<temporal_edge>>>();
        for (const auto&[v,v_vertex]:*vertex_map) {
            for (const auto& [neighbor_vertex_id,neighbor_edge_set]:*v_vertex->get_neighbor_map()) {
                for (const auto& edge:*neighbor_edge_set) {
                    edge_set->insert(edge);
                }
            }
        }
        return edge_set;
    }

    shared_ptr<unordered_set<shared_ptr<temporal_edge>>> temporal_graph::get_edge_set(uint32_t source_vertex_id, uint32_t destination_vertex_id) {
        auto source_vertex = get_vertex(source_vertex_id);
        return source_vertex ? source_vertex->get_temporal_edge_set(destination_vertex_id): shared_ptr<unordered_set<shared_ptr<temporal_edge>>>();
    }

    shared_ptr<unordered_map<uint32_t , shared_ptr<temporal_vertex>>> temporal_graph::get_vertex_map() {
        return this->vertex_map;
    }

    uint32_t temporal_graph::get_vertex_size()  {
        return vertex_map->size();
    }

    void temporal_graph::insert_edge(const shared_ptr<temporal_edge> &edge) {
        auto source_vertex_id = edge->get_source_vertex_id();
        auto destination_vertex_id = edge->get_destination_vertex_id();
        auto source_vertex = get_vertex(source_vertex_id);
        if (!source_vertex) {
            source_vertex = make_shared<temporal_vertex>(source_vertex_id);
            vertex_map->insert({source_vertex_id, source_vertex});
        }
        source_vertex->insert_edge(destination_vertex_id, edge);
        auto destination_vertex = get_vertex(destination_vertex_id);
        if (!destination_vertex) {
            destination_vertex = make_shared<temporal_vertex>(destination_vertex_id);
            vertex_map->insert({destination_vertex_id, destination_vertex});
        }
        destination_vertex->insert_edge(source_vertex_id, edge);
    }


    void temporal_graph::insert_edge(const shared_ptr<unordered_set<shared_ptr<temporal_edge>>> &edge_set) {
        for (const auto& edge:*edge_set) {
            insert_edge(edge);
        }
    }

    void temporal_graph::insert_vertex(uint32_t vertex_id) {
        vertex_map->insert({vertex_id, make_shared<temporal_vertex>(vertex_id)});
    }

    void temporal_graph::insert_vertex(const shared_ptr<temporal_vertex>& vertex) {
        vertex_map->insert({vertex->get_vertex_id(), vertex});
    }

    void temporal_graph::remove_vertex(uint32_t vertex_id) {
        vertex_map->erase(vertex_id);
    }

    void temporal_graph::remove_edge(const shared_ptr<temporal_edge> &edge) {
        auto source_vertex = get_vertex(edge->get_source_vertex_id());
        source_vertex->remove_edge(edge->get_destination_vertex_id(), edge);
        if (source_vertex->get_temporal_edge_size(edge->get_destination_vertex_id()) == 0) {
            source_vertex->remove_neighbor_vertex(edge->get_destination_vertex_id());
        }
        auto destination_vertex = get_vertex(edge->get_destination_vertex_id());
        destination_vertex->remove_edge(edge->get_source_vertex_id(), edge);
        if (destination_vertex->get_temporal_edge_size(edge->get_source_vertex_id()) == 0) {
            destination_vertex->remove_neighbor_vertex(edge->get_source_vertex_id());
        }
    }

    void temporal_graph::remove_edges(uint32_t source_vertex_id, uint32_t destination_vertex_id) {
        auto source_vertex = get_vertex(source_vertex_id);
        source_vertex->remove_neighbor_vertex(destination_vertex_id);

        auto destination_vertex = get_vertex(destination_vertex_id);
        destination_vertex->remove_neighbor_vertex(source_vertex_id);
    }
}