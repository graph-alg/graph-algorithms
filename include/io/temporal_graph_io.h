/******************************************************************************************************************
* @copyright  : Boost Software License
* @file       : io_utility.h
* @brief      : An implementation of functions for temporal_graph input/output
* @version    : 1.0
* @date       : 2021/09/07
******************************************************************************************************************/

#pragma once
#include "io/io_utility.h"

namespace scnu{
    class temporal_graph_io {
    public:
        static shared_ptr<vector<shared_ptr<temporal_edge>>>
        get_edge_vector(const string &path, const string &file_name);

        static shared_ptr<temporal_graph> load_graph(const shared_ptr<vector<shared_ptr<temporal_edge>>> &edge_vector);

        static shared_ptr<temporal_graph> load_graph(const shared_ptr<vector<shared_ptr<temporal_edge>>> &edge_vector,
                                                     const shared_ptr<thread_pool>& pool);

        static void store_graph(const string &input_path, const string &output_path,
                                const shared_ptr<thread_pool>& pool);

        static void output_unique_graph(const string &input_path, const string &output_path,
                                        const shared_ptr<thread_pool>& pool);

        /**
        * @details convert edge collection to csv format
        * @param edge_vector
        * @param output_file
        */
        template <typename T>
        static void output_csv_file(const T edge_container, const string &output_file) {
            ofstream output_stream(output_file);
            for (const auto &e:*edge_container) {
                output_stream << e->get_source_vertex_id() << ',' << e->get_destination_vertex_id() << ',' <<e->get_weight() <<',' << e->get_timestamp() << '\n';
            }
            output_stream.close();
        }
    };
}



