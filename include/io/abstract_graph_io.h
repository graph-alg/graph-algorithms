/******************************************************************************************************************
* @copyright  : Boost Software License
* @file       : abstract_graph_io.h
* @brief      : an implementation of input/output for abstract graphs
* @version    : 1.1
* @date       : 2020/10/16
******************************************************************************************************************/

#pragma once

#include "io/io_utility.h"

namespace scnu {
    /**
     * @details A graph I/O class contains some input and output functions
     */
    class abstract_graph_io {
    public:
        static shared_ptr<vector<shared_ptr<abstract_edge>>>
        get_edge_vector(const string &path, const string &file_name);

        static shared_ptr<vector<shared_ptr<abstract_edge>>>
        get_shrink_edge_vector(const string &path, const string &file_name, double rate);

        static shared_ptr<abstract_graph> load_graph(const shared_ptr<vector<shared_ptr<abstract_edge>>> &edge_vector);

        static shared_ptr<abstract_graph> load_graph(const shared_ptr<vector<shared_ptr<abstract_edge>>> &edge_vector, uint32_t thread_number);

        static shared_ptr<abstract_graph> load_graph(const shared_ptr<abstract_graph> &other_graph, uint32_t thread_number);

        static void output_core_number(const string &path, const string &file_name,
                                       const shared_ptr<unordered_map<uint32_t , uint32_t>> &vertex_core_map);

        static void output_csv_file(const string &input_path, const string &output_path,uint32_t thread_number);



        /**
         * @details convert edge collection to csv format
         * @param edge_vector
         * @param output_file
         */
        template <typename T>
        static void output_csv_file(const T edge_container, const string &output_file) {
            ofstream output_stream(output_file);
            for (const auto &e:*edge_container) {
                output_stream << e->get_source_vertex_id() << ',' << e->get_destination_vertex_id() << '\n';
            }
            output_stream.close();
        }
    };
}


