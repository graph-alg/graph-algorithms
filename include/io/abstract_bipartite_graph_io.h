/******************************************************************************************************************
* @copyright  : Boost Software License
* @file       : abstract_bipartite_graph_io.h
* @brief      : An I/O class for bipartite graph input and output
* @version    : 1.0
* @date       : 2020/09/06
******************************************************************************************************************/

#pragma once
#include "graph/abstract_bipartite_graph.h"
#include "io/io_utility.h"

namespace scnu
{
    class abstract_bipartite_graph_io
    {
    public:
        static shared_ptr<vector<shared_ptr<abstract_bipartite_edge>>> get_edge_vector(const string &path,
                                                                        const string &file_name);

        static shared_ptr<abstract_bipartite_graph> load_graph(const shared_ptr<vector<shared_ptr<abstract_bipartite_edge>>>& edge_vector);

        static shared_ptr<abstract_bipartite_graph> load_graph(const shared_ptr<vector<shared_ptr<abstract_bipartite_edge>>> &edge_vector, uint32_t thread_number);

        static void store_graph(const string &input_path, const string &output_path, uint32_t thread_number);


        /**
         * @details convert edge collection to csv format
         * @param edge_vector
         * @param output_file
         */
        template <typename T>
        static void output_csv_file(const T edge_container, const string &output_file) {
            ofstream output_stream(output_file);
            for (const auto &e:*edge_container) {
                output_stream << e->get_left_vertex_id() << ',' << e->get_right_vertex_id() << '\n';
            }
            output_stream.close();
        }
    };
}


