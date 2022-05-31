
#include "core_test/core_test.h"

int main(int argc, char **argv) {
    if(argc < 4){
        std::cout<<"Usage: input path, file name, and thread number!";
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
    auto graph = abstract_graph_io::load_graph(edge_vector,thread_number);
    auto vertex_core_map1 = make_shared<unordered_map<uint32_t,uint32_t>>();
    {
        simple_timer decomposition_timer;
        k_max = basic_core_decomposition::decompose(graph, vertex_core_map1,thread_number);
        auto decomposition_time = decomposition_timer.get_elapse_second();

        LOG(logger,LOG_RANK::INFO) << "Base decomposition," << decomposition_time << "\n";
    }

    /**
     * @brief test!
     */
//    auto vertex_core_map2 = make_shared<unordered_map<uint32_t,uint32_t>>();{
//        simple_timer t2;
//        bin_sort_core_decomposition::decompose(graph, vertex_core_map2);
//        auto bin_decomposition_time = t2.get_elapse_second();
//        LOG(logger,LOG_RANK::INFO)<<"Bin decomposition,"<<bin_decomposition_time<<"\n";
//    }
//
//    if (container_compare::same_associative_map(vertex_core_map1, vertex_core_map2)) {
//        LOG(logger,LOG_RANK::INFO)<< "OK!"<<"\n";
//    }

    LOG(logger,LOG_RANK::INFO) << "NumberOfVertices," << graph->get_vertex_number()<<"\n";
    LOG(logger,LOG_RANK::INFO) << "NumberOfEdges," << graph->get_edge_number()<<"\n";
    LOG(logger,LOG_RANK::INFO) << "MaximalDegree," << graph->get_maximal_degree()<<"\n";
    LOG(logger,LOG_RANK::INFO) << "kMax," << k_max<<"\n";
}