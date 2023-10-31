
#include "bipartite_core_test/bipartite_core_test.h"

int main(int argc, char **argv) {
    if (argc < 4) {
        std::cout << "Usage: please input path, file name and thread_number!";
    }
    string path = argv[1];
    string input_file_name = argv[2];
    auto thread_number = std::stoul(argv[3]);
    /**
     *@brief prepare log file
     **/
    string log_path = path + "log/";
    auto directory = std::filesystem::path(log_path);
    if (!exists(directory)) {
        create_directory(directory);
    }

    string log_file_name = log_path + "decomposition_test.log";
    scnu::simple_logger logger(log_file_name);

    auto edge_vector = abstract_bipartite_graph_io::get_edge_vector(path, input_file_name);

    auto B = shared_ptr<abstract_bipartite_graph>();
    {
        auto pool = make_shared<thread_pool>(thread_number);
        B = abstract_bipartite_graph_io::load_graph(edge_vector, pool);
    }

    LOG(logger, LOG_RANK::INFO) << input_file_name << "\n";

    auto index_left_index_map = make_shared<unordered_map<uint32_t, shared_ptr<bipartite_core_left_store_index>>>();
    auto index_right_index_map = make_shared<unordered_map<uint32_t, shared_ptr<bipartite_core_right_store_index>>>();
    uint32_t delta = 0;
    {
        simple_timer decomposition_timer;
        auto left_mutex_map = make_shared<unordered_map<uint32_t, shared_ptr<mutex>>>();
        auto right_mutex_map = make_shared<unordered_map<uint32_t, shared_ptr<mutex>>>();

        auto pool = make_shared<thread_pool>(thread_number);

        share_bipartite_core_decomposition::init(B, left_mutex_map, right_mutex_map,
                                                 index_left_index_map, index_right_index_map, pool);
        delta = share_bipartite_core_decomposition::decompose(B, left_mutex_map, right_mutex_map,
                                                              index_left_index_map,
                                                              index_right_index_map, pool);
        auto decomposition_time = decomposition_timer.get_elapse_second();

        LOG(logger, LOG_RANK::INFO) << "Share Decomposition," << decomposition_time << "\n";
    }


//    {
//        auto index_left_index_map2 = make_shared<unordered_map<uint32_t, shared_ptr<bipartite_core_left_store_index>>>();
//        auto index_right_index_map2 = make_shared<unordered_map<uint32_t, shared_ptr<bipartite_core_right_store_index>>>();
//        simple_timer decomposition_timer;
//        auto left_mutex_map = make_shared<unordered_map<uint32_t, shared_ptr<mutex>>>();
//        auto right_mutex_map = make_shared<unordered_map<uint32_t, shared_ptr<mutex>>>();
//
//        auto pool = make_shared<thread_pool>(thread_number);
//
//        share_bipartite_core_decomposition::init(B, left_mutex_map, right_mutex_map,
//                                             index_left_index_map, index_right_index_map, pool);
//
//        share_bipartite_core_decomposition::decompose2(B, index_left_index_map2,
//                                                       index_right_index_map2, pool);
//        auto decomposition_time = decomposition_timer.get_elapse_second();
//
//        if (bipartite_core_compare::same(index_left_index_map, index_right_index_map,
//                                         index_left_index_map2, index_right_index_map2)) {
//            LOG(logger, LOG_RANK::INFO) << "Share Decomposition2," << decomposition_time << "\n";
//        }
//    }
//
//    {
//        auto branch_left_index_map = make_shared<unordered_map<uint32_t, shared_ptr<bipartite_core_left_store_index>>>();
//        auto branch_right_index_map = make_shared<unordered_map<uint32_t, shared_ptr<bipartite_core_right_store_index>>>();
//        simple_timer decomposition_timer;
//
//        auto left_mutex_map = make_shared<unordered_map<uint32_t, shared_ptr<mutex>>>();
//        auto right_mutex_map = make_shared<unordered_map<uint32_t, shared_ptr<mutex>>>();
//
//        auto branch_core_order_map = make_shared<bipartite_core_order_index>();
//        auto branch_core_rem_degree_map = make_shared<bipartite_core_rem_degree_index>();
//        auto branch_core_degree_map = make_shared<bipartite_core_degree_index>();
//
//        auto pool = make_shared<thread_pool>(thread_number);
//        branch_bipartite_core_decomposition::init(B, left_mutex_map, right_mutex_map, branch_left_index_map, branch_right_index_map, pool);
//        branch_bipartite_core_decomposition::decompose(B,left_mutex_map, right_mutex_map, branch_left_index_map,
//                                                       branch_right_index_map, branch_core_order_map,
//                                                       branch_core_rem_degree_map, branch_core_degree_map,
//                                                       pool);
//        auto decomposition_time = decomposition_timer.get_elapse_second();
//        if (bipartite_core_compare::same(branch_left_index_map, branch_right_index_map,
//                                         index_left_index_map, index_right_index_map)) {
//            LOG(logger, LOG_RANK::INFO) << "Branch Decomposition," << decomposition_time << "\n";
//        }
//    }
//
//    {
//        auto branch_left_index_map = make_shared<unordered_map<uint32_t, shared_ptr<bipartite_core_left_store_index>>>();
//        auto branch_right_index_map = make_shared<unordered_map<uint32_t, shared_ptr<bipartite_core_right_store_index>>>();
//        simple_timer decomposition_timer;
//
//        auto left_mutex_map = make_shared<unordered_map<uint32_t, shared_ptr<mutex>>>();
//        auto right_mutex_map = make_shared<unordered_map<uint32_t, shared_ptr<mutex>>>();
//
//        auto branch_core_order_map = make_shared<bipartite_core_order_index>();
//        auto branch_core_rem_degree_map = make_shared<bipartite_core_rem_degree_index>();
//        auto branch_core_degree_map = make_shared<bipartite_core_degree_index>();
//
//        auto pool = make_shared<thread_pool>(thread_number);
//        branch_bipartite_core_decomposition::init(B, left_mutex_map, right_mutex_map, branch_left_index_map, branch_right_index_map, pool);
//        branch_bipartite_core_decomposition::decompose2(B,left_mutex_map, right_mutex_map, branch_left_index_map,
//                                                       branch_right_index_map, branch_core_order_map,
//                                                       branch_core_rem_degree_map, branch_core_degree_map,
//                                                       pool);
//        auto decomposition_time = decomposition_timer.get_elapse_second();
//        if (bipartite_core_compare::same(branch_left_index_map, branch_right_index_map,
//                                         index_left_index_map, index_right_index_map)) {
//            LOG(logger, LOG_RANK::INFO) << "Branch Decomposition2," << decomposition_time << "\n";
//        }
//    }

    LOG(logger, LOG_RANK::INFO) << "NumberOfLeftVertices," << B->get_left_vertex_number() << "\n";
    LOG(logger, LOG_RANK::INFO) << "NumberOfRightVertices," << B->get_right_vertex_number() << "\n";
    LOG(logger, LOG_RANK::INFO) << "NumberOfEdges," << B->get_edge_number() << "\n";
    LOG(logger, LOG_RANK::INFO) << "MaximalLeftVertexDegree," << B->get_maximal_left_vertex_degree() << "\n";
    LOG(logger, LOG_RANK::INFO) << "MaximalRightVertexDegree," << B->get_maximal_right_vertex_degree() << "\n";
    LOG(logger, LOG_RANK::INFO) << "AverageLeftVertexDegree," << B->get_average_left_vertex_degree() << "\n";
    LOG(logger, LOG_RANK::INFO) << "AverageRightVertexDegree," << B->get_average_right_vertex_degree() << "\n";
    LOG(logger, LOG_RANK::INFO) << "Delta," << delta << "\n";

    return 0;
}

