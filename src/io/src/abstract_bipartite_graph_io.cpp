
#include "io/abstract_bipartite_graph_io.h"

namespace scnu {
    /**
     * @details load all edges from a given file
     * @param path
     * @param file_name
     * @return
     */
    shared_ptr<vector<shared_ptr<abstract_bipartite_edge>>>
    abstract_bipartite_graph_io::get_edge_vector(const string &path,
                                                 const string &file_name) {
        ifstream input_stream(path + file_name);
        string line;
        auto edge_vector = make_shared<vector<shared_ptr<abstract_bipartite_edge>>>();
        while (getline(input_stream, line).good()) {
            if (line.empty()) {
                continue;
            }
            auto line_vector = string_algorithm::regex_split(line, ",");
            auto left_vertex_id = stoul(line_vector->at(0));
            auto right_vertex_id = stoul(line_vector->at(1));
            auto edge = make_shared<abstract_bipartite_edge>(left_vertex_id, right_vertex_id);
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
    shared_ptr<abstract_bipartite_graph>
    abstract_bipartite_graph_io::load_graph(const shared_ptr<vector<shared_ptr<abstract_bipartite_edge>>> &edge_vector) {
        auto graph = make_shared<abstract_bipartite_graph>();
        for (const auto &edge: *edge_vector) {
            graph->insert_edge(edge);
        }
        return graph;
    }

    shared_ptr<abstract_bipartite_graph>
    abstract_bipartite_graph_io::load_graph(const shared_ptr<vector<shared_ptr<abstract_bipartite_edge>>> &edge_vector, uint32_t thread_number) {
        auto graph = make_shared<abstract_bipartite_graph>();
        auto left_vertex_map = graph->get_left_vertex_map();
        auto right_vertex_map = graph->get_right_vertex_map();

        auto left_vertex_mutex_map = make_shared<unordered_map<uint32_t, shared_ptr<mutex>>>();
        auto right_vertex_mutex_map = make_shared<unordered_map<uint32_t, shared_ptr<mutex>>>();

        auto pool = make_shared<thread_pool>(thread_number);
        pool->submit_task([=]{
            for (const auto &edge: *edge_vector) {
                auto l = edge->get_left_vertex_id();
                if(!left_vertex_map->count(l)){
                    left_vertex_map->insert({l, make_shared<abstract_left_vertex>(l)});
                    left_vertex_mutex_map->insert({l, make_shared<mutex>()});
                }
            }
        });
        pool->submit_task([=]{
            for (const auto &edge: *edge_vector) {
                auto r = edge->get_right_vertex_id();
                if(!right_vertex_map->count(r)){
                    right_vertex_map->insert({r, make_shared<abstract_right_vertex>(r)});
                    right_vertex_mutex_map->insert({r, make_shared<mutex>()});
                }
            }
        });
        pool->barrier();

        for(const auto&edge:*edge_vector){
            pool->submit_task([=]{
                auto l = edge->get_left_vertex_id();
                auto r = edge->get_right_vertex_id();

                auto l_vertex = graph->get_left_vertex(l);
                left_vertex_mutex_map->at(l)->lock();
                l_vertex->insert_edge(r, edge);
                left_vertex_mutex_map->at(l)->unlock();

                auto r_vertex = graph->get_right_vertex(r);
                right_vertex_mutex_map->at(r)->lock();
                r_vertex->insert_edge(l, edge);
                right_vertex_mutex_map->at(r)->unlock();
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
    void abstract_bipartite_graph_io::output_csv_file(const string &input_path, const string &output_path, uint32_t thread_number) {
        auto directory = path(input_path);

        thread_pool pool(thread_number);
        for (const auto &file_iter:std::filesystem::directory_iterator(input_path)) {
            if(!std::filesystem::is_regular_file(file_iter)){
                continue;
            }
            pool.submit_task([=] {
                unordered_set<pair<uint32_t, uint32_t>, hash_pair, equal_pair> line_set;


                auto file_name = file_iter.path().filename().string();
                ifstream input_stream(input_path + file_name);

                uint32_t new_left_vertex_id = 0;
                uint32_t new_right_vertex_id = 0;
                auto left_vertex_id_map = make_shared<unordered_map<uint32_t, uint32_t>>();
                auto right_vertex_id_map = make_shared<unordered_map<uint32_t, uint32_t>>();
                string line;
                while (getline(input_stream, line).good()) {
                    if (line.empty() || line[0] == '%' || line[0] == '#') {
                        continue;
                    }

                    line = string_algorithm::replace_all(line, "[ |\\t|\\r]+", ",");
                    auto line_list = string_algorithm::regex_split(line, ",");

                    auto left_vertex_id = stoul(line_list->at(0));
                    auto right_vertex_id = stoul(line_list->at(1));

                    /**
                     * @brief renumber the left vertex id
                     */
                    if (!left_vertex_id_map->count(left_vertex_id)) {
                        left_vertex_id_map->insert({left_vertex_id, ++new_left_vertex_id});
                    }

                    /**
                     * @brief renumber the right vertex id
                     */
                    if (!right_vertex_id_map->count(right_vertex_id)) {
                        right_vertex_id_map->insert({right_vertex_id, ++new_right_vertex_id});
                    }

                    if (!line_set.count({left_vertex_id, right_vertex_id})) {
                        line_set.insert({left_vertex_id, right_vertex_id});
                    }
                }
                input_stream.close();

                multiset<shared_ptr<abstract_bipartite_edge>, abstract_bipartite_edge_compare> edge_set;
                for (const auto&[l, r]:line_set) {
                    auto left_vertex_id = left_vertex_id_map->at(l);
                    auto right_vertex_id = left_vertex_id_map->size() + right_vertex_id_map->at(r);
                    auto edge = make_shared<abstract_bipartite_edge>(left_vertex_id, right_vertex_id);
                    edge_set.insert(edge);
                }

                auto begin_index = file_name.find_first_of('.');
                file_name = file_name.substr(begin_index + 1);
                ofstream output_stream(output_path + file_name);
                for (const auto &e:edge_set) {
                    output_stream << e->get_left_vertex_id() << ',' << e->get_right_vertex_id() << '\n';
                }
                output_stream.close();
            });
        }
        pool.barrier();
    }
}


