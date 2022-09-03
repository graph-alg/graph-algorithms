
#include "graph/weighted_graph.h"

namespace scnu{
    weighted_graph::weighted_graph() {
        vertex_map = make_shared<unordered_map<uint32_t,shared_ptr<weighted_vertex>>>();
    }

    weighted_graph::weighted_graph(const shared_ptr<unordered_map<uint32_t , shared_ptr<weighted_vertex>>> &other_vertex_map) {
        this->vertex_map = other_vertex_map;
    }

    weighted_graph::weighted_graph(const shared_ptr<weighted_graph> &graph):weighted_graph() {
        for (const auto&[v,v_vertex]: *graph->get_vertex_map()) {
            vertex_map->insert({v, make_shared<weighted_vertex>(v_vertex)});
        }
    }

    weighted_graph::~weighted_graph() = default;

    bool weighted_graph::empty(){
        return vertex_map->empty();
    }

    bool weighted_graph::exist_edge(const shared_ptr<weighted_edge>& edge) {
        auto source_vertex = get_vertex(edge->get_source_vertex_id());
        if (!source_vertex) {
            return false;
        }
        return source_vertex->exist_edge(edge->get_destination_vertex_id());
    }

    bool weighted_graph::exist_edge(uint32_t source_vertex_id, uint32_t destination_vertex_id) {
        auto source_vertex = get_vertex(source_vertex_id);
        if (!source_vertex) {
            return false;
        }
        return source_vertex->get_neighbor_map()->count(destination_vertex_id);
    }

    bool weighted_graph::exist_vertex(uint32_t vertex_id) {
        return vertex_map->count(vertex_id);
    }

    uint32_t weighted_graph::get_edge_size() {
        uint32_t edge_size = 0;
        for (const auto&[u,u_vertex]:*vertex_map) {
            edge_size += u_vertex->get_degree();
        }
        return edge_size / 2;
    }


    shared_ptr<weighted_vertex> weighted_graph::get_vertex(uint32_t vertex_id) {
        return vertex_map->count(vertex_id) ? vertex_map->at(vertex_id) : shared_ptr<weighted_vertex>();
    }

    shared_ptr<weighted_edge> weighted_graph::get_edge(uint32_t source_vertex_id, uint32_t destination_vertex_id){
        if(vertex_map->count(source_vertex_id)){
            return vertex_map->at(source_vertex_id)->get_edge(destination_vertex_id);
        }
        return shared_ptr<weighted_edge>();
    }


    shared_ptr<unordered_set<shared_ptr<weighted_edge>>> weighted_graph::get_edge_set() {
        auto edge_set = make_shared<unordered_set<shared_ptr<weighted_edge>>>();
        for (const auto&[v,v_vertex]:*vertex_map) {
            for (const auto& [neighbor_vertex_id,e]:*v_vertex->get_neighbor_map()) {
                edge_set->insert(e);
            }
        }
        return edge_set;
    }

    shared_ptr<unordered_map<uint32_t , shared_ptr<weighted_vertex>>> weighted_graph::get_vertex_map() {
        return this->vertex_map;
    }

    uint32_t weighted_graph::get_vertex_size()  {
        return vertex_map->size();
    }

    void weighted_graph::insert_edge(const shared_ptr<weighted_edge> &edge) {
        auto source_vertex_id = edge->get_source_vertex_id();
        auto destination_vertex_id = edge->get_destination_vertex_id();
        auto source_vertex = get_vertex(source_vertex_id);
        if (!source_vertex) {
            source_vertex = make_shared<weighted_vertex>(source_vertex_id);
            vertex_map->insert({source_vertex_id, source_vertex});
        }
        source_vertex->insert_edge(destination_vertex_id, edge);
        auto destination_vertex = get_vertex(destination_vertex_id);
        if (!destination_vertex) {
            destination_vertex = make_shared<weighted_vertex>(destination_vertex_id);
            vertex_map->insert({destination_vertex_id, destination_vertex});
        }
        destination_vertex->insert_edge(source_vertex_id, edge);
    }


    void weighted_graph::insert_edge(const shared_ptr<unordered_set<shared_ptr<weighted_edge>>> &edge_set) {
        for (const auto& edge:*edge_set) {
            insert_edge(edge);
        }
    }

    void weighted_graph::insert_vertex(uint32_t vertex_id) {
        vertex_map->insert({vertex_id, make_shared<weighted_vertex>(vertex_id)});
    }

    void weighted_graph::insert_vertex(const shared_ptr<weighted_vertex>& vertex) {
        vertex_map->insert({vertex->get_vertex_id(), vertex});
    }

    void weighted_graph::remove_vertex(uint32_t vertex_id) {
        vertex_map->erase(vertex_id);
    }

    void weighted_graph::remove_edge(const shared_ptr<weighted_edge> &edge) {

    }

    void weighted_graph::remove_edges(uint32_t source_vertex_id, uint32_t destination_vertex_id) {
        auto source_vertex = get_vertex(source_vertex_id);
        source_vertex->remove_edge(destination_vertex_id);
        if(source_vertex->get_degree() == 0){
            vertex_map->erase(source_vertex->get_vertex_id());
        }

        auto destination_vertex = get_vertex(destination_vertex_id);
        destination_vertex->remove_edge(source_vertex_id);
        if(destination_vertex->get_degree() == 0){
            vertex_map->erase(destination_vertex->get_vertex_id());
        }
    }
}