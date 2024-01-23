
#include "truss_test/truss_test.h"

int main(int argc, char **argv) {
    if(argc < 4){
        std::cout<<"Usage: input path, file name, and thread number!";
        exit(0);
    }

    string path = argv[1];
    string input_file_name = argv[2];
    uint32_t  thread_number = std::stoul(argv[3]);

    /**
     *@brief prepare log file
     **/
    string log_path = path + "log/";
    auto directory = std::filesystem::path(log_path);
    if (!exists(directory)) {
        create_directory(directory);
    }

    string log_file_name = log_path  + "decomposition_test.log";
    simple_logger logger(log_file_name);

    auto edge_vector = abstract_graph_io::get_edge_vector(path, input_file_name);
    LOG(logger, LOG_RANK::INFO)  << input_file_name << "\n";

    uint32_t  k_max = 0;
    auto G = shared_ptr<abstract_graph>();
    {
        auto pool = make_shared<thread_pool>(thread_number);
        G = abstract_graph_io::load_graph(edge_vector, pool);
    }


    auto edge_truss_map1 = make_shared<unordered_map<shared_ptr<abstract_edge>,uint32_t>>();
    {
        simple_timer decomposition_timer;
        auto edge_mutex_map = make_shared<unordered_map<shared_ptr<abstract_edge>, shared_ptr<mutex>>>();
        auto edge_rank_map = make_shared<unordered_map<shared_ptr<abstract_edge>, uint32_t>>();
        auto edge_support_map = make_shared<unordered_map<shared_ptr<abstract_edge>, uint32_t>>();

        auto edge_truss_support_map = make_shared<unordered_map<shared_ptr<abstract_edge>, uint32_t>>();

        auto pool = make_shared<thread_pool>(thread_number);
        basic_truss_decomposition::init(G, edge_mutex_map, edge_rank_map, edge_support_map, edge_truss_map1, pool);
        k_max = basic_truss_decomposition::decompose(G, edge_mutex_map, edge_rank_map, edge_support_map,
                                                     edge_truss_map1, pool);
        auto decomposition_time = decomposition_timer.get_elapse_second();

        LOG(logger, LOG_RANK::INFO) << "Base decomposition," << decomposition_time << "\n";
    }

    /**
     * @brief test!
     */
//    uint32_t k_max2 = 0;
//    auto edge_truss_map2 = make_shared<unordered_map<shared_ptr<abstract_edge>,uint32_t>>();
//    {
//        simple_timer t2;
//        auto edge_support_map = make_shared<unordered_map<shared_ptr<abstract_edge>, uint32_t>>();
//        bin_sort_truss_decomposition::init(G, edge_support_map, edge_truss_map2);
//        k_max2 = bin_sort_truss_decomposition::decompose(G, edge_support_map, edge_truss_map2);
//        auto bin_decomposition_time = t2.get_elapse_second();
//        LOG(logger,LOG_RANK::INFO)<<"Bin decomposition,"<<bin_decomposition_time<<"\n";
//    }
//
//    if (truss_compare::same_associative_map(edge_truss_map1, edge_truss_map2)) {
//        LOG(logger,LOG_RANK::INFO)<< "OK!"<<"\n";
//    }

    LOG(logger,LOG_RANK::INFO) << "NumberOfVertices," << G->get_vertex_number() << "\n";
    LOG(logger,LOG_RANK::INFO) << "NumberOfEdges," << G->get_edge_number() << "\n";
    LOG(logger,LOG_RANK::INFO) << "MaximalDegree," << G->get_maximal_degree() << "\n";
    LOG(logger,LOG_RANK::INFO) << "kMax," << k_max<<"\n";
}
