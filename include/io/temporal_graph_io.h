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
                                                     uint32_t thread_number);

        static void output_csv_file(const string &input_path, const string &output_path, uint32_t thread_number);

        static void unique_output_csv_file(const string &input_path, const string &output_path, uint32_t thread_number);
    };
}



