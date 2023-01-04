
#include "io/temporal_graph_io.h"

namespace scnu{
    /**
     * @details get a vector of edges from the given file
     * @param path
     * @param file_name
     * @return
    */
    shared_ptr<vector<shared_ptr<temporal_edge>>> temporal_graph_io::get_edge_vector(const string &path,
                                                                                     const string &file_name) {
        ifstream input_stream(path + file_name);
        string line;
        auto edge_vector = make_shared<vector<shared_ptr<temporal_edge>>>();
        while (input_stream.good()) {
            getline(input_stream, line);
            if (line.empty()) {
                continue;
            }
            auto line_vector = string_algorithm::regex_split(line, ",");
            auto source_vertex_id = stoul(line_vector->at(0));
            auto destination_vertex_id = stoul(line_vector->at(1));
            auto weight =  stod(line_vector->at(2));
            auto edge_timestamp = stoul(line_vector->at(3));
            auto edge = make_shared<temporal_edge>(source_vertex_id, destination_vertex_id, weight,edge_timestamp);
            edge_vector->push_back(edge);
        }
        input_stream.close();
        return edge_vector;
    }

    shared_ptr<temporal_graph>
    temporal_graph_io::load_graph(const shared_ptr<vector<shared_ptr<temporal_edge>>> &edge_vector) {
        auto graph = make_shared<temporal_graph>();
        for (const auto &edge: *edge_vector) {
            graph->insert_edge(edge);
        }
        return graph;
    }

    shared_ptr<temporal_graph>
    temporal_graph_io::load_graph(const shared_ptr<vector<shared_ptr<temporal_edge>>> &edge_vector,
                                  const shared_ptr<thread_pool>& pool) {
        auto graph = make_shared<temporal_graph>();
        auto vertex_map = graph->get_vertex_map();
        auto vertex_mutex_map = make_shared<unordered_map<uint32_t, shared_ptr<mutex>>>();
        for (const auto &edge: *edge_vector) {
            auto u = edge->get_source_vertex_id();
            if(!vertex_map->count(u)){
                vertex_map->insert({u, make_shared<temporal_vertex>(u)});
                vertex_mutex_map->insert({u, make_shared<mutex>()});
            }
            auto v = edge->get_destination_vertex_id();
            if(!vertex_map->count(v)){
                vertex_map->insert({v, make_shared<temporal_vertex>(v)});
                vertex_mutex_map->insert({v, make_shared<mutex>()});
            }
        }

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
                v_vertex->insert_edge(u,edge);
                vertex_mutex_map->at(v)->unlock();
            });
        }
        pool->barrier();

        return graph;
    }

    /**
     * @details convert original line to csv format
     * @param input_path
     * @param output_path
     */
    void temporal_graph_io::store_graph(const string &input_path, const string &output_path,
                                        const shared_ptr<thread_pool>& pool) {
        auto directory = path(input_path);

        for (auto &file_iter:std::filesystem::directory_iterator(input_path)) {
            if(!std::filesystem::is_regular_file(file_iter)){
                continue;
            }
            pool->submit_task([=] {
                auto edge_map = make_shared<map<uint32_t, shared_ptr<unordered_set<shared_ptr<temporal_edge>, hash_temporal_edge,equal_temporal_edge>>>>();

                auto vertex_id_map = make_shared<unordered_map<uint32_t, uint32_t>>();
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

                    auto u = stoul(line_list->at(0));
                    auto v = stoul(line_list->at(1));

                    /**
                     * @brief remove self-loop edge
                     */
                    if (u == v) {
                        continue;
                    }

                    /**
                     * @brief by default, the timestamp is located in forth column
                     */
                    double w = std::stod(line_list->at(2));
                    uint32_t t = std::stoul(line_list->at(3));

                    if (!vertex_id_map->count(u)) {
                        vertex_id_map->insert({u, vertex_id++});
                    }

                    if (!vertex_id_map->count(v)) {
                        vertex_id_map->insert({v, vertex_id++});
                    }

                    u = vertex_id_map->at(u);
                    v = vertex_id_map->at(v);

                    /**
                     * @brief swap u and v
                     */
                    if (u > v) {
                        swap(u, v);
                    }

                    if(!edge_map->count(t)){
                        edge_map->insert({t, make_shared<unordered_set<shared_ptr<temporal_edge>,hash_temporal_edge, equal_temporal_edge>>()});
                    }

                    auto e = make_shared<temporal_edge>(u, v, w,t);
                    edge_map->at(t)->insert(e);
                }
                input_stream.close();

                auto begin_index = file_name.find_last_of('.');
                file_name = file_name.substr(begin_index + 1);
                ofstream output_stream(output_path + file_name);
                for(const auto& [t, t_set]:*edge_map){
                    for (const auto &e:*t_set) {
                        output_stream << e->get_source_vertex_id() << ','
                                      << e->get_destination_vertex_id() << ','
                                      << e->get_weight()<< ','
                                      << t << '\n';
                    }
                }
                output_stream<< '\n';
                output_stream.close();
            });
        }
        pool->barrier();
    }

    void temporal_graph_io::output_unique_graph(const string &input_path, const string &output_path,
                                                const shared_ptr<thread_pool>& pool) {
        auto directory = path(input_path);

        for (auto &file_iter:std::filesystem::directory_iterator(input_path)) {
            if(!std::filesystem::is_regular_file(file_iter)){
                continue;
            }
            pool->submit_task([=] {
                auto edge_set =  make_shared<unordered_set<shared_ptr<abstract_edge>, hash_abstract_edge, equal_abstract_edge>>();
                auto edge_map = make_shared<map<uint32_t, shared_ptr<unordered_set<shared_ptr<abstract_edge>>>>>();

                auto vertex_id_map = make_shared<unordered_map<uint32_t, uint32_t> >();
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

                    /**
                     * @brief remove self-loop edge
                     */
                    if (u == v) {
                        continue;
                    }

                    /**
                     * @brief by default, the timestamp is located in forth column
                     */
                    uint32_t t =  std::stoull(line_list->at(3));

                    if (!vertex_id_map->count(u)) {
                        vertex_id_map->insert({u, vertex_id++});
                    }

                    if (!vertex_id_map->count(v)) {
                        vertex_id_map->insert({v, vertex_id++});
                    }

                    u = vertex_id_map->at(u);
                    v = vertex_id_map->at(v);

                    /**
                     * @brief swap u and v
                     */
                    if (u > v) {
                        swap(u, v);
                    }


                    auto e  = make_shared<abstract_edge>(u, v);
                    if(!edge_set->count(e)){
                        edge_set->insert(e);
                        if(!edge_map->count(t)){
                            edge_map->insert({t, make_shared<unordered_set<shared_ptr<abstract_edge>>>()});
                        }
                        edge_map->at(t)->insert(e);
                    }
                }
                input_stream.close();
                edge_set->clear();

                auto begin_index = file_name.find_last_of('.');
                file_name = file_name.substr(begin_index+1);
                ofstream output_stream(output_path + file_name);
                for(const auto& [t, e_set]:*edge_map){
                    for(const auto &e :*e_set){
                        output_stream <<  e->get_source_vertex_id() << ','
                                      << e->get_destination_vertex_id() << ','
                                      << t << '\n';
                    }
                }
                output_stream<< '\n';
                output_stream.close();
            });
        }
        pool->barrier();
    }
}