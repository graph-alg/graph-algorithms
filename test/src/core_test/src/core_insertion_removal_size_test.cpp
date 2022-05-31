
#include "core_test/core_test.h"


int main(int argc, char **argv) {
    if(argc < 4){
        std::cout<<"Usage: please input path, file name, and thread number!";
    }

    string path = argv[1];
    string input_file_name = argv[2];
    auto thread_number = std::stoul(argv[3]);
    /**
    * @brief prepare the log file
    */
    string log_path = path + "log/";
    string log_file_name = log_path  + "insertion_removal_size_test.log";
    simple_logger logger(log_file_name);

    auto total_edge_vector = abstract_graph_io::get_edge_vector(path, input_file_name);
    /**
     * @brief shuffle the edges, test!
     */
//    {
//        auto rd = make_shared<random_device>();
//        shuffle(edge_vector->begin(),edge_vector->end(),*random_generator::get_default_engine(rd));
//    }


    vector<double> edge_ratio{0.2, 0.4, 0.6, 0.8, 1.0};
    for (auto r :edge_ratio) {
        auto size = uint32_t (total_edge_vector->size()*r)+1;
        auto edge_vector = make_shared<vector<shared_ptr<abstract_edge>>>();
        for(uint32_t i=0;i< min(size, uint32_t (total_edge_vector->size()));++i){
            edge_vector->push_back(total_edge_vector->at(i));
        }

        LOG(logger,LOG_RANK::INFO) << input_file_name << "," << r << '\n';

        auto previous_edge_vector = make_shared<vector<shared_ptr<abstract_edge>>>();
        auto insertion_edge_vector = make_shared<vector<shared_ptr<abstract_edge>>>();
        auto removal_edge_vector = make_shared<vector<shared_ptr<abstract_edge>>>();

        auto m = 5000;
        for (auto i = 0; i < edge_vector->size(); ++i) {
            if (i < m) {
                removal_edge_vector->emplace_back(edge_vector->at(i));
            }
            if (i < edge_vector->size() - m) {
                previous_edge_vector->emplace_back(edge_vector->at(i));
            } else {
                insertion_edge_vector->emplace_back(edge_vector->at(i));
            }
        }

        /**
         * @details prepare core number of vertices in previous graph
         */
        auto previous_graph = abstract_graph_io::load_graph(previous_edge_vector,thread_number);
        auto previous_vertex_core_map = make_shared<unordered_map<uint32_t, uint32_t>>();
        uint32_t previous_k_max_value = 0;
        {
            previous_k_max_value = basic_core_decomposition::decompose(previous_graph, previous_vertex_core_map, thread_number);
        }

        /**
         * @brief decomposition method for checking results
         */
        auto contrastive_vertex_core_map = make_shared<unordered_map<uint32_t, uint32_t>>();
        {

            previous_graph->insert_edge_collection(insertion_edge_vector);
            previous_graph->remove_edge_collection(removal_edge_vector);

            simple_timer decomposition_timer;
            basic_core_decomposition::decompose(previous_graph, contrastive_vertex_core_map, thread_number);
            double decomposition_time = decomposition_timer.get_elapse_second();

            LOG(logger, LOG_RANK::INFO) << "Decomposition," << decomposition_time << "\n";

            /**
             * @brief restore the graph
             */
            previous_graph->remove_edge_collection(insertion_edge_vector);
            previous_graph->insert_edge_collection(removal_edge_vector);
        }


        /**
        * @brief joint maintenance
        */
        {
            auto previous_vertex_core_map1 = container_copy::to_unordered_map<uint32_t, uint32_t>(
                    previous_vertex_core_map);
            auto EI1 = container_convert::to_unordered_set<shared_ptr<abstract_edge>>(insertion_edge_vector);
            auto ER1 = container_convert::to_unordered_set<shared_ptr<abstract_edge>>(removal_edge_vector);

            simple_timer maintenance_timer;
            jes_core_maintenance::joint_insert(previous_graph, EI1, previous_vertex_core_map1, thread_number);
            jes_core_maintenance::joint_delete(previous_graph, ER1, previous_vertex_core_map1, thread_number);
            auto maintenance_time = maintenance_timer.get_elapse_second();

            if (core_compare::same_associative_map(contrastive_vertex_core_map, previous_vertex_core_map1)) {
                LOG(logger, LOG_RANK::INFO) << "Joint Maintenance," << maintenance_time << "\n";
            }

            /**
             * @brief restore the graph
             */
            previous_graph->remove_edge_collection(insertion_edge_vector);
            previous_graph->insert_edge_collection(removal_edge_vector);
        }

        /**
         * @brief parallel maintenance
         */
        {
            auto previous_vertex_core_map2 = container_copy::to_unordered_map<uint32_t,uint32_t>(previous_vertex_core_map);
            auto insertion_edge_set2 = container_convert::to_unordered_set<shared_ptr<abstract_edge>>(
                    insertion_edge_vector);
            auto removal_edge_set2 = container_convert::to_unordered_set<shared_ptr<abstract_edge>>(
                    removal_edge_vector);
            auto previous_k_max2 = make_shared<uint32_t>(previous_k_max_value);

            simple_timer maintenance_timer;
            quasi_core_maintenance::insert(previous_graph, insertion_edge_set2,
                                           previous_vertex_core_map2, previous_k_max2, thread_number);
            quasi_core_maintenance::remove(previous_graph, removal_edge_set2,
                                           previous_vertex_core_map2, previous_k_max2, thread_number);
            auto maintenance_time = maintenance_timer.get_elapse_second();

            if (core_compare::same_associative_map(previous_vertex_core_map2, contrastive_vertex_core_map)) {
                LOG(logger, LOG_RANK::INFO) << "Parallel Maintenance," << maintenance_time << "\n";
            }

            /**
             * @brief restore the graph
             */
            previous_graph->remove_edge_collection(insertion_edge_vector);
            previous_graph->insert_edge_collection(removal_edge_vector);
        }
    }
}



