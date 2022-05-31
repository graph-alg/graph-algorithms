
#include "core_test/core_test.h"

int main(int argc, char **argv) {
    if(argc < 4){
        std::cout<<"Usage: please input path, file name, and thread number!";
    }
    string path = argv[1];
    string input_file_name = argv[2];
    uint32_t thread_number = std::stoul(argv[3]);
    /**
     * @brief prepare the log file
     */
    string log_path = path + "log/";
    string log_file_name = log_path + "insertion_thread_test.log";
    simple_logger logger(log_file_name);


    auto edge_vector = abstract_graph_io::get_edge_vector(path, input_file_name);
    /**
    * @brief shuffle the edges
    */
//    {
//        auto rd = make_shared<random_device>();
//        shuffle(edge_vector->begin(),edge_vector->end(),*random_generator::get_default_engine(rd));
//    }


    /**
     * @brief the number of inserted edges
     */
    uint32_t m =  1000000;

    auto previous_edge_vector = make_shared<vector<shared_ptr<abstract_edge>>>();
    auto insertion_edge_vector = make_shared<vector<shared_ptr<abstract_edge>>>();

    for (auto i = 0; i < edge_vector->size(); ++i) {
        if (i < edge_vector->size() - m) {
            previous_edge_vector->emplace_back(edge_vector->at(i));
        } else {
            insertion_edge_vector->emplace_back(edge_vector->at(i));
        }
    }

    /**
     * @prepare core number of vertices in previous graph
     */
    auto previous_graph = abstract_graph_io::load_graph(previous_edge_vector, thread_number);
    auto previous_vertex_core_map = make_shared<unordered_map<uint32_t,uint32_t>>();
    uint32_t previous_k_max_value = 0;
    {
        previous_k_max_value = basic_core_decomposition::decompose(previous_graph, previous_vertex_core_map, thread_number);
    }

    /**
     * @brief decomposition for checking results
     */
    auto contrastive_vertex_core_map = make_shared<unordered_map<uint32_t,uint32_t>>();
    double decomposition_time = 0;
    {
        previous_graph->insert_edge_collection(insertion_edge_vector);

        simple_timer decomposition_timer;
        basic_core_decomposition::decompose(previous_graph, contrastive_vertex_core_map, thread_number);
        decomposition_time = decomposition_timer.get_elapse_second();

        /*
         * @brief restore the graph
         */
        previous_graph->remove_edge_collection(insertion_edge_vector);
    }


    vector<uint32_t> thread_array{1, 2, 4, 6, 10, 14};

    for (const auto& sub_thread_number :thread_array) {
        LOG(logger,LOG_RANK::INFO) << input_file_name <<  "," << m << "," << sub_thread_number << '\n';

        LOG(logger,LOG_RANK::INFO) << "Decomposition("<< thread_number <<")," << decomposition_time << "\n";
        /**
         * @brief joint insert
         */
        {
            auto previous_vertex_core_map1 = container_copy::to_unordered_map<uint32_t,uint32_t>(previous_vertex_core_map);
            auto EI1 = container_copy::to_unordered_set<shared_ptr<abstract_edge>>(insertion_edge_vector);

            simple_timer maintenance_timer;
            jes_core_maintenance::joint_insert(previous_graph, EI1, previous_vertex_core_map1, sub_thread_number);
            auto maintenance_time = maintenance_timer.get_elapse_second();

            if (core_compare::same_associative_map(contrastive_vertex_core_map, previous_vertex_core_map1)) {
                LOG(logger,LOG_RANK::INFO) << "Joint Insertion," << maintenance_time << "\n";
            }

            /*
             * @brief restore the graph
             */
            previous_graph->remove_edge_collection(insertion_edge_vector);
        }

        /**
         * @brief parallel insert
         */
        {
            auto previous_vertex_core_map2 = container_copy::to_unordered_map<uint32_t,uint32_t>(previous_vertex_core_map);
            auto insertion_edge_set2 = container_copy::to_unordered_set<shared_ptr<abstract_edge>>(insertion_edge_vector);
            auto previous_k_max2 = make_shared<uint32_t>(previous_k_max_value);

            simple_timer maintenance_timer;
            quasi_core_maintenance::insert(previous_graph, insertion_edge_set2, previous_vertex_core_map2,
                                           previous_k_max2, sub_thread_number);
            auto maintenance_time = maintenance_timer.get_elapse_second();

            if(core_compare::same_associative_map(previous_vertex_core_map2, contrastive_vertex_core_map))
            {
                LOG(logger,LOG_RANK::INFO) << "Parallel Insertion," << maintenance_time << "\n";
            }

            /*
             * @brief restore the graph
             */
            previous_graph->remove_edge_collection(insertion_edge_vector);
        }
    }
}

