
#include "io/abstract_graph_io.h"

namespace scnu {
    /**
     * @details get a vector of edges from the given file
     * @param path
     * @param file_name
     * @return
     */
    shared_ptr<vector<shared_ptr<abstract_edge>>> abstract_graph_io::get_edge_vector(const string &path,
                                                                                     const string &file_name) {
        ifstream input_stream(path + file_name);
        string line;
        auto edge_vector = make_shared<vector<shared_ptr<abstract_edge>>>();
        while (input_stream.good()) {
            getline(input_stream, line);
            if (line.empty()) {
                continue;
            }
            auto line_vector = string_algorithm::regex_split(line, ",");
            auto source_vertex_id = stoul(line_vector->at(0));
            auto destination_vertex_id = stoul(line_vector->at(1));
            auto edge = make_shared<abstract_edge>(source_vertex_id, destination_vertex_id);
            edge_vector->push_back(edge);
        }
        input_stream.close();
        return edge_vector;
    }

    shared_ptr<vector<shared_ptr<abstract_edge>>> abstract_graph_io::get_shrink_edge_vector(const string &path,
                                                                                            const string &file_name,
                                                                                            double rate) {
        ifstream input_stream(path + file_name);
        string line;
        auto edge_vector = make_shared<vector<shared_ptr<abstract_edge>>>();
        while (input_stream.good()) {
            getline(input_stream, line);
            if (line.empty()) {
                continue;
            }
            auto line_vector = string_algorithm::regex_split(line, ",");
            auto source_vertex_id = stoul(line_vector->at(0));
            auto destination_vertex_id = stoul(line_vector->at(1));
            auto edge = make_shared<abstract_edge>(source_vertex_id, destination_vertex_id);
            edge_vector->push_back(edge);
        }
        input_stream.close();
        auto size = uint32_t (edge_vector->size() * rate) + 1;
        auto sub_edge_vector = make_shared<vector<shared_ptr<abstract_edge>>>(size, shared_ptr<abstract_edge>());
        for(uint32_t i = 0; i < size; ++i){
            sub_edge_vector->at(i) = edge_vector->at(i);
        }
        return sub_edge_vector;
    }

    /**
     * @details load a graph from a vector of edges
     * @param edge_vector
     * @return
     */
    shared_ptr<abstract_graph>
    abstract_graph_io::load_graph(const shared_ptr<vector<shared_ptr<abstract_edge>>> &edge_vector) {
        auto graph = make_shared<abstract_graph>();
        for (const auto &edge: *edge_vector) {
            graph->insert_edge(edge);
        }
        return graph;
    }

    shared_ptr<abstract_graph>
    abstract_graph_io::load_graph(const shared_ptr<vector<shared_ptr<abstract_edge>>> &edge_vector, uint32_t thread_number) {
        auto graph = make_shared<abstract_graph>();
        auto vertex_map = graph->get_vertex_map();
        auto vertex_mutex_map = make_shared<unordered_map<uint32_t, shared_ptr<mutex>>>();
        for (const auto &edge: *edge_vector) {
            auto u = edge->get_source_vertex_id();
            if(!vertex_map->count(u)){
                vertex_map->insert({u, make_shared<abstract_vertex>(u)});
                vertex_mutex_map->insert({u, make_shared<mutex>()});
            }
            auto v = edge->get_destination_vertex_id();
            if(!vertex_map->count(v)){
                vertex_map->insert({v, make_shared<abstract_vertex>(v)});
                vertex_mutex_map->insert({v, make_shared<mutex>()});
            }
        }

        auto pool = make_shared<thread_pool>();
        uint32_t task_number = edge_vector->size() % thread_number == 0? edge_vector->size()/ thread_number: edge_vector->size()/ thread_number+1;
        for(uint32_t i = 0;i< thread_number;++i){
            pool->submit_task([=]{
                auto end_point = std::min((i+1)*task_number,uint32_t(edge_vector->size()));
                for(uint32_t j = i * task_number;j < end_point;++j){
                    auto edge = edge_vector->at(j);
                    auto u = edge->get_source_vertex_id();
                    auto v = edge->get_destination_vertex_id();

                    auto u_vertex = vertex_map->at(u);
                    vertex_mutex_map->at(u)->lock();
                    u_vertex->insert_edge(v,edge);
                    vertex_mutex_map->at(u)->unlock();

                    auto v_vertex = vertex_map->at(v);
                    vertex_mutex_map->at(v)->lock();
                    v_vertex->insert_edge(u,edge);
                    vertex_mutex_map->at(v)->unlock();
                }
            });
        }
        pool->barrier();
        return graph;
    }

    shared_ptr<abstract_graph>
    abstract_graph_io::load_graph(const shared_ptr<abstract_graph> &other_graph, uint32_t thread_number) {
        auto graph = make_shared<abstract_graph>();
        auto vertex_map = graph->get_vertex_map();
        for (const auto &[u,u_vertex]: *other_graph->get_vertex_map()) {
            vertex_map->insert({u, shared_ptr<abstract_vertex>()});
        }
        auto pool = make_shared<thread_pool>(thread_number);
        for (const auto &p: *other_graph->get_vertex_map()) {
            pool->submit_task([=]{
                auto [u,other_u_vertex] = p;
                vertex_map->at(u) = make_shared<abstract_vertex>(other_u_vertex);
            });
        }
        pool->barrier();
        return graph;
    }

    /**
     * @details output the core number of each vertex
     * @param path
     * @param file_name
     * @param vertex_core_map
     */
    void abstract_graph_io::output_core_number(const string &path, const string &file_name,
                                               const shared_ptr<unordered_map<uint32_t, uint32_t>> &vertex_core_map) {
        ofstream output_stream(path + file_name);
        output_stream << "vertex_id,core_number\n";
        for (const auto&[v, core_number]:*vertex_core_map) {
            output_stream << v << "," << core_number << "\n";
        }
        output_stream.close();
    }

    /**
     * @details convert original line to csv format
     * @remarks the method will remove multiple edges
     * @param input_path
     * @param output_path
     */
    void abstract_graph_io::output_csv_file(const string &input_path, const string &output_path, uint32_t thread_number) {
        auto directory = path(input_path);

        thread_pool pool(thread_number);
        for (auto &file_iter:std::filesystem::directory_iterator(input_path)) {
//            if(!std::filesystem::is_regular_file(file_iter)){
//                continue;
//            }
            pool.submit_task([=] {
                unordered_set<pair<uint32_t, uint32_t>, hash_pair, equal_pair> edge_set;

                unordered_map<uint32_t, uint32_t> vertex_id_map;
                uint32_t vertex_id = 1;

                const auto &file = file_iter.path();
                auto file_name = file.filename().string();
                ifstream input_stream(input_path + file_name);

                auto begin_index = file_name.find_last_of('.');
                file_name = file_name.substr(begin_index+1);
                ofstream output_stream(output_path + file_name);

                string line;
                while (getline(input_stream, line).good()) {
                    if (line.empty() || line[0] == '%' || line[0] == '#') {
                        continue;
                    }

                    line = string_algorithm::replace_all(line, "[ |\\t|\\r]+", ",");
                    auto line_list = string_algorithm::regex_split(line, ",");

                    auto source_vertex_id = stoul(line_list->at(0));
                    auto destination_vertex_id = stoul(line_list->at(1));

                    /**
                     * @brief remove self-loop edge
                     */
                    if (source_vertex_id == destination_vertex_id) {
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

                    auto edge_pair = make_pair(source_vertex_id, destination_vertex_id);
                    if (!edge_set.count(edge_pair)) {
                        edge_set.insert(edge_pair);
                        output_stream << edge_pair.first << ',' << edge_pair.second<< '\n';
                    }
                }
                input_stream.close();

                output_stream.close();
            });
        }
        pool.barrier();
    }
}



