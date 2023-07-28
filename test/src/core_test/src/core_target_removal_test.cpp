
#include "core_test/core_test.h"

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
    string log_file_name = log_path + "target_removal_test.log";
    simple_logger logger(log_file_name);

    auto edge_vector = abstract_graph_io::get_edge_vector(path, input_file_name);
    /**
     * @brief shuffle the edges, test!
     */
//    {
//        auto rd = make_shared<random_device>();
//        shuffle(total_edge_vector->begin(),total_edge_vector->end(),*random_generator::get_default_engine(rd));
//    }

    auto global_vertex_core_map = make_shared<unordered_map<uint32_t , uint32_t>>();
    {
        auto pool = make_shared<thread_pool>(thread_number);
        auto G = abstract_graph_io::load_graph(edge_vector, pool);
        basic_core_decomposition::decompose(G,global_vertex_core_map, thread_number);
    }

    uint32_t m = 1000000;

    vector<uint32_t> edge_ratio{10, 20, 30, 40, 50};
    for (const auto& r:edge_ratio) {
        LOG(logger, LOG_RANK::INFO) <<input_file_name <<"," << r << '\n';

        auto candidate_edge_vector = make_shared<vector<shared_ptr<abstract_edge>>>();
        for(const auto &e:*edge_vector){
            auto u = e->get_source_vertex_id();
            auto v = e->get_destination_vertex_id();
            auto min_k = std::min(global_vertex_core_map->at(u), global_vertex_core_map->at(v));
            if(min_k == r){
                candidate_edge_vector->emplace_back(e);
            }
        }
        auto removal_edge_vector = make_shared<vector<shared_ptr<abstract_edge>>>();
        for(uint32_t i = 0;i<candidate_edge_vector->size();++i){
            if(i < m){
                removal_edge_vector->emplace_back(candidate_edge_vector->at(i));
            }
        }


        auto previous_graph = shared_ptr<abstract_graph>(0);
        auto previous_vertex_core_map = make_shared<unordered_map<uint32_t,uint32_t>>();
        uint32_t previous_k_max_value = 0;
        {
            auto pool = make_shared<thread_pool>(thread_number);
            previous_graph = abstract_graph_io::load_graph(edge_vector,pool);
            previous_k_max_value = basic_core_decomposition::decompose(previous_graph, previous_vertex_core_map, thread_number);
        }


        /**
        * @brief decomposition for checking results
        */
        auto contrastive_vertex_core_map = make_shared<unordered_map<uint32_t, uint32_t>>();
        {
            previous_graph->remove_edge_collection(removal_edge_vector);

            simple_timer decomposition_timer;
            basic_core_decomposition::decompose(previous_graph, contrastive_vertex_core_map, thread_number);
            double decomposition_time = decomposition_timer.get_elapse_second();

            LOG(logger,LOG_RANK::INFO) << "Decomposition," << decomposition_time << "\n";

            /**
             * @brief restore the graph
             */
            previous_graph->insert_edge_collection(removal_edge_vector);
        }

        /**
         * @brief joint insert
         */
        {
            auto previous_vertex_core_map1 = container_copy::to_unordered_map<uint32_t,uint32_t>(previous_vertex_core_map);
            auto ER1 = container_copy::to_unordered_set<shared_ptr<abstract_edge>>(removal_edge_vector);

            simple_timer maintenance_timer;
            jes_core_maintenance::joint_delete(previous_graph, ER1, previous_vertex_core_map1, thread_number);
            auto maintenance_time = maintenance_timer.get_elapse_second();

            if (core_compare::same_associative_map(contrastive_vertex_core_map, previous_vertex_core_map1)) {
                LOG(logger,LOG_RANK::INFO) << "JES Removal," << maintenance_time << "\n";
            }

            /**
             * @brief restore the graph
             */
            previous_graph->insert_edge_collection(removal_edge_vector);
        }


        /**
         * @brief parallel remove methods
         */
        {
            auto previous_vertex_core_map2 = container_copy::to_unordered_map<uint32_t,uint32_t>(previous_vertex_core_map);
            auto removal_edge_set2 = container_convert::to_unordered_set<shared_ptr<abstract_edge>>(removal_edge_vector);
            auto previous_k_max2 = make_shared<uint32_t>(previous_k_max_value);

            simple_timer maintenance_timer2;
            quasi_core_maintenance::remove(previous_graph, removal_edge_set2, previous_vertex_core_map2,
                                           previous_k_max2, thread_number);
            auto maintenance_time2 = maintenance_timer2.get_elapse_second();

            if (core_compare::same_associative_map(contrastive_vertex_core_map, previous_vertex_core_map2)) {
                LOG(logger,LOG_RANK::INFO) << "Parallel Removal," << maintenance_time2 << "\n";
            }

            /**
             * @brief restore the graph
             */
            previous_graph->insert_edge_collection(removal_edge_vector);
        }
    }
}

