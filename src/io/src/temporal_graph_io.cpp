
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
            auto edge_timestamp = stoul(line_vector->at(2));
            auto edge = make_shared<temporal_edge>(source_vertex_id, destination_vertex_id,edge_timestamp);
            edge_vector->push_back(edge);
        }
        input_stream.close();
        return edge_vector;
    }

    /**
     * @details load a graph from a vector of edges
     * @param edge_vector
     * @return
     */
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
                                  uint32_t thread_number) {
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
    void temporal_graph_io::output_csv_file(const string &input_path, const string &output_path, uint32_t thread_number) {
        auto directory = path(input_path);

        thread_pool pool(thread_number);
        for (auto &file_iter:std::filesystem::directory_iterator(input_path)) {
//            if(!std::filesystem::is_regular_file(file_iter)){
//                continue;
//            }
            pool.submit_task([=] {
                unordered_set<shared_ptr<temporal_edge>,hash_temporal_edge,equal_temporal_edge> line_set;

                multiset<shared_ptr<temporal_edge>, temporal_edge_compare> edge_set;

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

                    auto source_vertex_id = stoull(line_list->at(0));
                    auto destination_vertex_id = stoull(line_list->at(1));

                    /**
                     * @brief remove self-loop edge
                     */
                    if (source_vertex_id == destination_vertex_id) {
                        continue;
                    }

                    /**
                     * @brief by default, the timestamp is located in forth column
                     */
                    uint32_t edge_timestamp = 0;
                    if (line_list->size() > 2 && !line_list->at(2).empty()) {
                        edge_timestamp = std::stoull(line_list->at(3));
                    }

                    /**
                     * @brief skip the edge lack of time information
                     */
                    if(edge_timestamp == 0){
                        continue;
                    }

                    if (!vertex_id_map.count(source_vertex_id)) {
                        vertex_id_map.insert({source_vertex_id, vertex_id++});
                    }

                    if (!vertex_id_map.count(destination_vertex_id)) {
                        vertex_id_map.insert({destination_vertex_id, vertex_id++});
                    }

                    source_vertex_id = vertex_id_map[source_vertex_id];
                    destination_vertex_id = vertex_id_map[destination_vertex_id];

                    /**
                     * @brief swap source_vertex_id and destination_vertex_id
                     */
                    if (source_vertex_id > destination_vertex_id) {
                        swap(source_vertex_id, destination_vertex_id);
                    }

                    auto e = make_shared<temporal_edge>(source_vertex_id, destination_vertex_id, edge_timestamp);
                    if (!line_set.count(e)) {
                        line_set.insert(e);
                        edge_set.insert(e);
                    }
                }
                input_stream.close();

                auto begin_index = file_name.find_last_of('.');
                file_name = file_name.substr(begin_index+1);
                ofstream output_stream(output_path + file_name);
                for (const auto &e:edge_set) {
                    output_stream << e->get_source_vertex_id() << ','
                                  << e->get_destination_vertex_id() << ','
                                  <<e->get_timestamp() << '\n';
                }
                output_stream<< '\n';
                output_stream.close();
            });
        }
        pool.barrier();
    }

    void temporal_graph_io::unique_output_csv_file(const string &input_path, const string &output_path, uint32_t thread_number) {
        auto directory = path(input_path);

        thread_pool pool(thread_number);
        for (auto &file_iter:std::filesystem::directory_iterator(input_path)) {
//            if(!std::filesystem::is_regular_file(file_iter)){
//                continue;
//            }
            pool.submit_task([=] {
                unordered_set<pair<uint32_t,uint32_t>,hash_pair,equal_pair> line_set;

                set<shared_ptr<temporal_edge>, temporal_edge_compare> edge_set;

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

                    auto source_vertex_id = stoull(line_list->at(0));
                    auto destination_vertex_id = stoull(line_list->at(1));

                    /**
                     * @brief remove self-loop edge
                     */
                    if (source_vertex_id == destination_vertex_id) {
                        continue;
                    }

                    /**
                     * @brief by default, the timestamp is located in forth column
                     */
                    uint32_t edge_timestamp = 0;
                    if (line_list->size() > 2 && !line_list->at(2).empty()) {
                        edge_timestamp = std::stoull(line_list->at(3));
                    }

                    /**
                     * @brief skip the edge lack of time information
                     */
                    if(edge_timestamp == 0){
                        continue;
                    }

                    if (!vertex_id_map.count(source_vertex_id)) {
                        vertex_id_map.insert({source_vertex_id, vertex_id++});
                    }

                    if (!vertex_id_map.count(destination_vertex_id)) {
                        vertex_id_map.insert({destination_vertex_id, vertex_id++});
                    }

                    source_vertex_id = vertex_id_map[source_vertex_id];
                    destination_vertex_id = vertex_id_map[destination_vertex_id];

                    /**
                     * @brief swap source_vertex_id and destination_vertex_id
                     */
                    if (source_vertex_id > destination_vertex_id) {
                        swap(source_vertex_id, destination_vertex_id);
                    }

                    auto p = make_pair(source_vertex_id,destination_vertex_id);
                    if (!line_set.count(p)) {
                        line_set.insert(p);
                        auto e = make_shared<temporal_edge>(source_vertex_id, destination_vertex_id, edge_timestamp);
                        edge_set.insert(e);
                    }
                }
                input_stream.close();

                auto begin_index = file_name.find_last_of('.');
                file_name = file_name.substr(begin_index+1);
                ofstream output_stream(output_path + file_name);
                for (const auto &e:edge_set) {
                    output_stream << e->get_source_vertex_id() << ','
                                  << e->get_destination_vertex_id() << ','
                                  <<e->get_timestamp() << '\n';
                }
                output_stream<< '\n';
                output_stream.close();
            });
        }
        pool.barrier();
    }
}