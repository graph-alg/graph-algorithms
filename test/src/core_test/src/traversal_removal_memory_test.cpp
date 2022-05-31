
#include "core_test/core_test.h"
#include <fstream>

int main(int argc, char **argv) {
    if(argc < 4){
        std::cout<<"Usage: please input path, file name, and thread number!";
    }

    string path = argv[1];
    string input_file_name = argv[2];
    uint32_t thread_number = std::stoul(argv[3]);
    /**
     * @brief prepare log file
     */
    string log_path = path + "log/";
    string log_file_name = log_path + "removal_memory_test.log";
    simple_logger logger(log_file_name);

    auto edge_vector = scnu::abstract_graph_io::get_edge_vector(path, input_file_name);
    /**
     * @brief shuffle the edges, test!
     */
//    {
//        auto rd = make_shared<random_device>();
//        shuffle(edge_vector->begin(), edge_vector->end(), *random_generator::get_default_engine(rd));
//    }


    /**
     * @brief the number of removed edges
     */
    uint32_t m = 10000;
    LOG(logger, LOG_RANK::INFO)  << input_file_name << "," << m << "\n";

    auto removal_edge_vector = make_shared<vector<shared_ptr<abstract_edge>>>();
    for (uint32_t i = 0; i < m; ++i) {
        removal_edge_vector->emplace_back(edge_vector->at(i));
    }

    /**
     * @details prepare core number of vertices in previous graph
     */
    auto previous_graph = abstract_graph_io::load_graph(edge_vector, thread_number);
    auto previous_vertex_core_map = make_shared<unordered_map<uint32_t, uint32_t>>();
    {
        basic_core_decomposition::decompose(previous_graph, previous_vertex_core_map,thread_number);
    }

    auto contrastive_vertex_core_map = make_shared<unordered_map<uint32_t, uint32_t>>();
    {
        previous_graph->remove_edge_collection(removal_edge_vector);

        simple_timer decomposition_timer;
        basic_core_decomposition::decompose(previous_graph, contrastive_vertex_core_map, thread_number);
        auto decomposition_time = decomposition_timer.get_elapse_second();

        LOG(logger, LOG_RANK::INFO) << "Decomposition," << decomposition_time << "\n";

        /**
         * @brief restore the graph
         */
        previous_graph->insert_edge_collection(removal_edge_vector);
    }



    /**
     * @brief traversal methods
     */
    {
        auto previous_vertex_core_map1 = container_copy::to_unordered_map<uint32_t,uint32_t>(previous_vertex_core_map);
        auto rcd1 = make_shared<unordered_map<uint32_t, shared_ptr<unordered_map<uint32_t, uint32_t>>>>();
        uint32_t n1 = 2;

        simple_timer maintenance_timer;
        traversal_core_maintenance::init_rcd(previous_graph, previous_vertex_core_map1, rcd1, n1);
        for (const auto &edge:*removal_edge_vector) {
            traversal_core_maintenance::remove(previous_graph, previous_vertex_core_map1, rcd1, n1,
                                               edge);
        }
        auto maintenance_time = maintenance_timer.get_elapse_second();

        if (core_compare::same_associative_map(contrastive_vertex_core_map, previous_vertex_core_map1)) {
            auto memory_size = process_information::get_memory(getpid());
            LOG(logger, LOG_RANK::INFO) << "Traversal Removal," << maintenance_time<<","<<double (memory_size)/1024/1024  << "\n";
        }

        /**
         * @brief restore the graph
         */
        previous_graph->insert_edge_collection(removal_edge_vector);
    }
}


