
#include "../../../include/io/weighted_graph_io.h"

namespace scnu{

    shared_ptr<vector<shared_ptr<weighted_edge>>> weighted_graph_io::get_edge_vector(const string &path,
                                                                                     const string &file_name) {
        ifstream input_stream(path + file_name);
        string line;
        auto edge_vector = make_shared<vector<shared_ptr<weighted_edge>>>();
        while (input_stream.good()) {
            getline(input_stream, line);
            if (line.empty()) {
                continue;
            }
            auto line_vector = string_algorithm::regex_split(line, ",");
            auto source_vertex_id = stoul(line_vector->at(0));
            auto destination_vertex_id = stoul(line_vector->at(1));
            auto weight = stoul(line_vector->at(2));
            auto edge = make_shared<weighted_edge>(source_vertex_id, destination_vertex_id,weight);
            edge_vector->push_back(edge);
        }
        input_stream.close();
        return edge_vector;
    }

    shared_ptr<weighted_graph>
    weighted_graph_io::load_graph(const shared_ptr<vector<shared_ptr<temporal_edge>>> &edge_vector) {
        auto graph = make_shared<weighted_graph>();
        for (const auto &edge: *edge_vector) {
            auto u = edge->get_source_vertex_id();
            auto v = edge->get_destination_vertex_id();
            auto weighted_e = graph->get_edge(u,v);
            if(weighted_e){
                weighted_e->set_weight(weighted_e->get_weight() + 1);
            }else
            {
                weighted_e = make_shared<weighted_edge>(u,v,1);
                graph->insert_edge(weighted_e);
            }
        }
        return graph;
    }

    shared_ptr<weighted_graph>
    weighted_graph_io::load_graph(const shared_ptr<vector<shared_ptr<weighted_edge>>> &edge_vector) {
        auto graph = make_shared<weighted_graph>();
        for (const auto &e: *edge_vector) {
            graph->insert_edge(e);
        }
        return graph;
    }

    shared_ptr<weighted_graph>
    weighted_graph_io::load_graph(const shared_ptr<vector<shared_ptr<temporal_edge>>> &edge_vector,
                                  uint32_t thread_number) {
        auto graph = make_shared<weighted_graph>();
        auto vertex_map = graph->get_vertex_map();
        auto vertex_mutex_map = make_shared<unordered_map<uint32_t, shared_ptr<mutex>>>();
        auto edge_map = make_shared<unordered_map<pair<uint32_t, uint32_t>, shared_ptr<weighted_edge>, hash_pair,equal_pair>>();

        for (const auto &edge: *edge_vector) {
            auto u = edge->get_source_vertex_id();
            if(!vertex_map->count(u)){
                vertex_map->insert({u, make_shared<weighted_vertex>(u)});
                vertex_mutex_map->insert({u, make_shared<mutex>()});
            }
            auto v = edge->get_destination_vertex_id();
            if(!vertex_map->count(v)){
                vertex_map->insert({v, make_shared<weighted_vertex>(v)});
                vertex_mutex_map->insert({v, make_shared<mutex>()});
            }

            if(edge_map->count({u,v})){
                auto e = edge_map->at({u,v});
                e->set_weight(e->get_weight() + 1);
            }else
            {
                edge_map->insert({{u,v}, make_shared<weighted_edge>(u,v,1)});
            }
        }
        auto pool = make_shared<thread_pool>(thread_number);
        for(const auto &p: *edge_map){
            pool->submit_task([=]{
                auto edge = p.second;
                auto u = edge->get_source_vertex_id();
                auto v = edge->get_destination_vertex_id();
                auto u_vertex = vertex_map->at(u);
                auto v_vertex = vertex_map->at(v);

                vertex_mutex_map->at(u)->lock();
                u_vertex->insert_edge(v,edge);
                vertex_mutex_map->at(u)->unlock();

                vertex_mutex_map->at(v)->lock();
                v_vertex->insert_edge(v,edge);
                vertex_mutex_map->at(v)->unlock();
            });
        }
        pool->barrier();

        return graph;
    }

    shared_ptr<weighted_graph>
    weighted_graph_io::load_graph(const shared_ptr<vector<shared_ptr<weighted_edge>>> &edge_vector,
                                  uint32_t thread_number) {
        auto graph = make_shared<weighted_graph>();
        auto vertex_map = graph->get_vertex_map();
        auto vertex_mutex_map = make_shared<unordered_map<uint32_t, shared_ptr<mutex>>>();
        for (const auto &edge: *edge_vector) {
            auto u = edge->get_source_vertex_id();
            if(!vertex_map->count(u)){
                vertex_map->insert({u, make_shared<weighted_vertex>(u)});
                vertex_mutex_map->insert({u, make_shared<mutex>()});
            }
            auto v = edge->get_destination_vertex_id();
            if(!vertex_map->count(v)){
                vertex_map->insert({v, make_shared<weighted_vertex>(v)});
                vertex_mutex_map->insert({v, make_shared<mutex>()});
            }
        }
        auto pool = make_shared<thread_pool>(thread_number);
        for(const auto &edge: *edge_vector){
            pool->submit_task([=]{
                auto u = edge->get_source_vertex_id();
                auto v = edge->get_destination_vertex_id();
                auto u_vertex = vertex_map->at(u);
                auto v_vertex = vertex_map->at(v);

                vertex_mutex_map->at(u)->lock();
                u_vertex->insert_edge(v,edge);
                vertex_mutex_map->at(u)->unlock();

                vertex_mutex_map->at(v)->lock();
                v_vertex->insert_edge(v,edge);
                vertex_mutex_map->at(v)->unlock();
            });
        }
        pool->barrier();

        return graph;
    }

    /**
     * @details store graph in csv format
     * @param input_path
     * @param output_path
     */
    void weighted_graph_io::store_graph(const string &input_path, const string &output_path, uint32_t thread_number) {
        auto directory = path(input_path);

        thread_pool pool(thread_number);
        for (auto &file_iter:std::filesystem::directory_iterator(input_path)) {
            if(!std::filesystem::is_regular_file(file_iter)){
                continue;
            }
            pool.submit_task([=] {
                auto edge_set = make_shared<unordered_set<shared_ptr<weighted_edge>>>();

                unordered_map<uint32_t, uint32_t> vertex_id_map;
                uint32_t vertex_id = 1;

                const auto &file = file_iter.path();
                auto file_name = file.filename().string();
                ifstream input_stream(input_path + file_name);

                string line;
                while (getline(input_stream, line).good()) {
                    if (line.empty() || line[0] == '%' || line[0] == '#') {
                        continue;
                    }

                    line = string_algorithm::replace_all(line, "[ |\\t|\\r]+", ",");
                    auto line_list = string_algorithm::regex_split(line, ",");

                    auto u = stoull(line_list->at(0));
                    auto v = stoull(line_list->at(1));
                    auto w = std::stod(line_list->at(2));

                    /**
                     * @brief remove self-loop edge
                     */
                    if (u == v) {
                        continue;
                    }

                    if (!vertex_id_map.count(u)) {
                        vertex_id_map.insert({u, vertex_id++});
                    }

                    if (!vertex_id_map.count(v)) {
                        vertex_id_map.insert({v, vertex_id++});
                    }

                    u = vertex_id_map[u];
                    v = vertex_id_map[v];

                    /**
                     * @brief swap u and v
                     */
                    if (u > v) {
                        swap(u, v);
                    }

                    auto e = make_shared<weighted_edge>(u,v,w);
                    edge_set->insert(e);
                }
                input_stream.close();

                auto begin_index = file_name.find_last_of('.');
                file_name = file_name.substr(begin_index + 1);
                ofstream output_stream(output_path + file_name);
                for(const auto& e:*edge_set){
                    output_stream << e->get_source_vertex_id() << ',' << e->get_destination_vertex_id() << ',' << e->get_weight() << '\n';
                }
                output_stream<< '\n';
                output_stream.close();
            });
        }
        pool.barrier();
    }
}
