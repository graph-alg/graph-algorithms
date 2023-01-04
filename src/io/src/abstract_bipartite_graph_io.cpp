
#include "io/abstract_bipartite_graph_io.h"

namespace scnu {
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
    abstract_bipartite_graph_io::load_graph(const shared_ptr<vector<shared_ptr<abstract_bipartite_edge>>> &edge_vector,
                                            const shared_ptr<thread_pool>& pool) {
        auto graph = make_shared<abstract_bipartite_graph>();
        auto left_vertex_map = graph->get_left_vertex_map();
        auto right_vertex_map = graph->get_right_vertex_map();

        auto left_vertex_mutex_map = make_shared<unordered_map<uint32_t, shared_ptr<mutex>>>();
        auto right_vertex_mutex_map = make_shared<unordered_map<uint32_t, shared_ptr<mutex>>>();

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

                left_vertex_mutex_map->at(l)->lock();
                auto l_vertex = graph->get_left_vertex(l);
                l_vertex->insert_edge(r, edge);
                left_vertex_mutex_map->at(l)->unlock();

                right_vertex_mutex_map->at(r)->lock();
                auto r_vertex = graph->get_right_vertex(r);
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
    void abstract_bipartite_graph_io::store_graph(const string &input_path, const string &output_path,
                                                  const shared_ptr<thread_pool>& pool) {
        auto directory = path(input_path);

        for (const auto &file_iter:std::filesystem::directory_iterator(input_path)) {
            if(!std::filesystem::is_regular_file(file_iter)){
                continue;
            }
            pool->submit_task([=] {

                auto file_name = file_iter.path().filename().string();
                ifstream input_stream(input_path + file_name);

                auto edge_set = make_shared<unordered_set<shared_ptr<abstract_bipartite_edge>, hash_abstract_bipartite_edge, equal_abstract_bipartite_edge>>();

                uint32_t new_left_vertex_id = 0;
                uint32_t new_right_vertex_id = 0;
                auto l_map = make_shared<unordered_map<uint32_t, uint32_t>>();
                auto r_map = make_shared<unordered_map<uint32_t, uint32_t>>();
                string line;
                while (getline(input_stream, line).good()) {
                    if (line.empty() || line[0] == '%' || line[0] == '#') {
                        continue;
                    }

                    line = string_algorithm::replace_all(line, "[ |\\t|\\r]+", ",");
                    auto line_list = string_algorithm::regex_split(line, ",");

                    auto l = stoul(line_list->at(0));
                    auto r = stoul(line_list->at(1));

                    /**
                     * @brief renumber the left vertex id
                     */
                    if (!l_map->count(l)) {
                        l_map->insert({l, ++new_left_vertex_id});
                    }
                    l = l_map->at(l);

                    /**
                     * @brief renumber the right vertex id
                     */
                    if (!r_map->count(r)) {
                        r_map->insert({r, ++new_right_vertex_id});
                    }
                    r = r_map->at(r);

                    auto e = make_shared<abstract_bipartite_edge>(l,r);
                    edge_set->insert(e);
                }
                input_stream.close();

                auto begin_index = file_name.find_first_of('.');
                file_name = file_name.substr(begin_index + 1);
                ofstream output_stream(output_path + file_name);
                for (const auto&e:*edge_set) {
                    output_stream << e->get_left_vertex_id()  << ',' <<(l_map->size() + e->get_right_vertex_id()) << '\n';
                }
                output_stream.close();
            });
        }
        pool->barrier();
    }
}


