
#include "core_test/core_test.h"


int main(int argc, char **argv) {
    if(argc < 3){
        std::cout<<"Usage: please input path, file name, and thread number!";
    }
    string path = argv[1];
    string input_file_name = argv[2];
    uint32_t thread_number = std::stoul(argv[3]);
    /**
     * @brief prepare log file
     */
    string log_path = path + "log/";
    string log_file_name = log_path + "removal_thread_test.log";
    simple_logger logger(log_file_name);

    auto edge_vector = abstract_graph_io::get_edge_vector(path, input_file_name);
    /**
     * @brief shuffle the edges, test!
     */
//    {
//        auto rd = make_shared<random_device>();
//        shuffle(edge_vector->begin(),edge_vector->end(),*random_generator::get_default_engine(rd));
//    }

    /**
     * @brief the number of removed edges
     */
    uint32_t m =  1000000;

    auto removal_edge_vector = make_shared<vector<shared_ptr<abstract_edge>>>();
    for(uint32_t i =0;i < m; ++i) {
        removal_edge_vector->emplace_back(edge_vector->at(i));
    }

    /**
      * @prepare core number of vertices in previous graph
      */
    auto previous_graph = shared_ptr<abstract_graph>();
    auto previous_vertex_core_map = make_shared<unordered_map<uint32_t,uint32_t>>();
    uint32_t previous_k_max_value = 0;
    {
        auto pool = make_shared<thread_pool>(thread_number);
        previous_graph = abstract_graph_io::load_graph(edge_vector,pool);
        previous_k_max_value = basic_core_decomposition::decompose(previous_graph, previous_vertex_core_map, thread_number);
    }

    auto contrastive_vertex_core_map = make_shared<unordered_map<uint32_t,uint32_t>>();
    double decomposition_time = 0;
    {

        previous_graph->remove_edge_collection(removal_edge_vector);

        simple_timer decomposition_timer;
        basic_core_decomposition::decompose(previous_graph, contrastive_vertex_core_map, thread_number);
        decomposition_time = decomposition_timer.get_elapse_second();

        /*
         * @brief restore the graph
         */
        previous_graph->insert_edge_collection(removal_edge_vector);
    }



    vector<uint32_t> thread_array{1, 2, 4, 6, 10, 14};

    for (auto sub_thread_number:thread_array) {
        LOG(logger,LOG_RANK::INFO) << input_file_name << "," << m <<  "," << sub_thread_number << "\n";

        LOG(logger,LOG_RANK::INFO) << "Decomposition," << decomposition_time << "\n";
        /**
         * @brief jes remove
         */
        {
            auto previous_vertex_core_map1 = container_copy::to_unordered_map<uint32_t,uint32_t>(previous_vertex_core_map);
            auto ER1 = container_copy::to_unordered_set<shared_ptr<abstract_edge>>(removal_edge_vector);

            simple_timer maintenance_timer;
            jes_core_maintenance::joint_delete(previous_graph, ER1, previous_vertex_core_map1, sub_thread_number);
            auto maintenance_time = maintenance_timer.get_elapse_second();

            if (core_compare::same_associative_map(contrastive_vertex_core_map, previous_vertex_core_map1)) {
                LOG(logger,LOG_RANK::INFO) << "JES Removal," << maintenance_time << "\n";
            }

            /*
             * @brief restore the graph
             */
            previous_graph->insert_edge_collection(removal_edge_vector);
        }

        /**
         * @brief parallel remove
         */
        {
            auto previous_vertex_core_map2 = container_copy::to_unordered_map<uint32_t,uint32_t>(previous_vertex_core_map);
            auto removal_edge_set2 = container_convert::to_unordered_set<shared_ptr<abstract_edge>>(removal_edge_vector);
            auto previous_k_max2 = make_shared<uint32_t>(previous_k_max_value);

            simple_timer maintenance_timer;
            quasi_core_maintenance::remove(previous_graph, removal_edge_set2, previous_vertex_core_map2,
                                           previous_k_max2, sub_thread_number);
            auto maintenance_time = maintenance_timer.get_elapse_second();

            if (core_compare::same_associative_map(contrastive_vertex_core_map, previous_vertex_core_map2)) {
                LOG(logger,LOG_RANK::INFO) << "Parallel Removal," << maintenance_time << "\n";
            }

            /*
             * @brief restore the graph
             */
            previous_graph->insert_edge_collection(removal_edge_vector);
        }
    }
}

