/******************************************************************************************************************
* @copyright  : Boost Software License
* @file       : weighted_graph_io.h
* @brief      : An implementation of functions for weighted_graph input/output
* @version    : 1.0
* @date       : 2021/09/07
******************************************************************************************************************/
#pragma once
#include "io/io_utility.h"

namespace scnu{
    class weighted_graph_io {
    public:
        static shared_ptr<vector<shared_ptr<weighted_edge>>>
        get_edge_vector(const string &path, const string &file_name);

        static shared_ptr<weighted_graph> load_graph(const shared_ptr<vector<shared_ptr<temporal_edge>>> &edge_vector);

        static shared_ptr<weighted_graph> load_graph(const shared_ptr<vector<shared_ptr<weighted_edge>>> &edge_vector);

        static shared_ptr<weighted_graph> load_graph(const shared_ptr<vector<shared_ptr<temporal_edge>>> &edge_vector,
                                                     const shared_ptr<thread_pool>& pool);

        static shared_ptr<weighted_graph> load_graph(const shared_ptr<vector<shared_ptr<weighted_edge>>> &edge_vector,
                                                     const shared_ptr<thread_pool>& pool);

        static void store_graph(const string &input_path, const string &output_path, const shared_ptr<thread_pool>& pool);
    };
}
