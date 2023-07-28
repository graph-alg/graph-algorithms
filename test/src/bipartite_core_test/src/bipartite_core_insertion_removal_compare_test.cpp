
#include "bipartite_core_test/bipartite_core_test.h"


int main(int argc, char **argv){
    if(argc < 4){
        std::cout<<"Usage: please input path, file name, and thread number!";
    }

    string path = argv[1];
    string input_file_name = argv[2];
    uint32_t thread_number = std::stoul(argv[3]);
    /**
     *@brief prepare log file
     **/
    string log_path = path + "log/";
    auto directory = std::filesystem::path(log_path);

    string log_file_name = log_path  + "insertion_removal_compare_test.log";
    scnu::simple_logger logger(log_file_name);

    auto total_edge_vector = abstract_bipartite_graph_io::get_edge_vector(path, input_file_name);
//    {
//        auto rd = make_shared<random_device>();
//        shuffle(total_edge_vector->begin(), total_edge_vector->end(), *random_generator::get_default_engine(rd));
//    }


    auto total_previous_edge_vector = make_shared<vector<shared_ptr<abstract_bipartite_edge>>>();
    auto total_insertion_edge_vector = make_shared<vector<shared_ptr<abstract_bipartite_edge>>>();
    auto total_removal_edge_vector = make_shared<vector<shared_ptr<abstract_bipartite_edge>>>();

    for (uint32_t i = 0; i < total_edge_vector->size(); i++)
    {
        if(i < 1250)
        {
            total_removal_edge_vector->push_back(total_edge_vector->at(i));
        }

        if(i < uint32_t (total_edge_vector->size()) - 1250)
        {
            total_previous_edge_vector->push_back(total_edge_vector->at(i));
        } else
        {
            total_insertion_edge_vector->push_back(total_edge_vector->at(i));
        }
    }

    auto B = shared_ptr<abstract_bipartite_graph>();

    auto previous_left_index_map = make_shared<unordered_map<uint32_t, shared_ptr<bipartite_core_left_store_index>>>();
    auto previous_right_index_map = make_shared<unordered_map<uint32_t, shared_ptr<bipartite_core_right_store_index>>>();
    auto previous_delta = make_shared<uint32_t>(0);
    {


        auto pool = make_shared<thread_pool>(thread_number);
        B = abstract_bipartite_graph_io::load_graph(total_previous_edge_vector, pool);

        *previous_delta = branch_bipartite_core_decomposition::decompose(B, previous_left_index_map,
                                                                         previous_right_index_map, pool);
    }

    vector<uint32_t> rate_vector {250, 500, 750, 1000, 1250};
    for(const auto &m:rate_vector){
        LOG(logger, LOG_RANK::INFO) << input_file_name<< ","<<m<<"\n";

        auto insertion_edge_vector = make_shared<vector<shared_ptr<abstract_bipartite_edge>>>();
        for(uint32_t i = 0;i < m; ++i){
            insertion_edge_vector->push_back(total_insertion_edge_vector->at(i));
        }

        auto removal_edge_vector = make_shared<vector<shared_ptr<abstract_bipartite_edge>>>();
        for(uint32_t i = 0;i < total_removal_edge_vector->size(); ++i){
            if(i < m){
                removal_edge_vector->push_back(total_removal_edge_vector->at(i));
            }
        }

        /**
         * directly decompose the remain graph
         */
        auto contrastive_left_index_map = make_shared<unordered_map<uint32_t, shared_ptr<bipartite_core_left_store_index>>>();
        auto contrastive_right_index_map = make_shared<unordered_map<uint32_t, shared_ptr<bipartite_core_right_store_index>>>();
        {
            B->insert_edge_collection(insertion_edge_vector);
            B->remove_edge_collection(removal_edge_vector);

            auto pool = make_shared<thread_pool>(thread_number);
            branch_bipartite_core_decomposition::decompose(B, contrastive_left_index_map,contrastive_right_index_map, pool);

            B->insert_edge_collection(removal_edge_vector);
            B->remove_edge_collection(insertion_edge_vector);
        }

        /**
         * @brief branch maintenance
         */
        {
            auto branch_left_vertex_map = make_shared<unordered_map<uint32_t, shared_ptr<bipartite_core_left_store_index>>>();
            auto branch_right_vertex_map = make_shared<unordered_map<uint32_t, shared_ptr<bipartite_core_right_store_index>>>();
            auto inserted_edge_set = container_copy::to_unordered_set<shared_ptr<abstract_bipartite_edge>>(
                    insertion_edge_vector);
            auto removed_edge_set = container_copy::to_unordered_set<shared_ptr<abstract_bipartite_edge>>(
                    removal_edge_vector);
            for (const auto &[l, l_node]: *previous_left_index_map) {
                branch_left_vertex_map->insert({l, make_shared<bipartite_core_left_store_index>(l_node)});
            }
            for (const auto &[r, r_node]: *previous_right_index_map) {
                branch_right_vertex_map->insert({r, make_shared<bipartite_core_right_store_index>(r_node)});
            }

            simple_timer maintenance_timer;
            {
                auto new_left_index_map = make_shared<unordered_map<uint32_t, shared_ptr<bipartite_core_left_store_index>>>();
                auto new_right_index_map = make_shared<unordered_map<uint32_t, shared_ptr<bipartite_core_right_store_index>>>();

                auto pool = make_shared<thread_pool>(thread_number);

                branch_bipartite_core_maintenance::init(B, new_left_index_map, new_right_index_map, pool);
                branch_bipartite_core_maintenance::remove(B, removed_edge_set, branch_left_vertex_map,branch_right_vertex_map,
                                                          new_left_index_map, new_right_index_map, pool);
                branch_bipartite_core_maintenance::insert(B, inserted_edge_set, branch_left_vertex_map,branch_right_vertex_map,
                                                          new_left_index_map, new_right_index_map, pool);
            }
            auto maintenance_time = maintenance_timer.get_elapse_second();

            if (bipartite_core_compare::same(branch_left_vertex_map, branch_right_vertex_map, contrastive_left_index_map,contrastive_right_index_map)) {
                LOG(logger, LOG_RANK::INFO) << "Branch Maintenance," << maintenance_time << "\n";
            }

            B->insert_edge_collection(removal_edge_vector);
            B->remove_edge_collection(insertion_edge_vector);
        }

        /**
         * @brief branch maintenance
         */
        {
            auto branch_left_index_map = make_shared<unordered_map<uint32_t, shared_ptr<bipartite_core_left_store_index>>>();
            auto branch_right_index_map = make_shared<unordered_map<uint32_t, shared_ptr<bipartite_core_right_store_index>>>();
            auto inserted_edge_set = container_copy::to_unordered_set<shared_ptr<abstract_bipartite_edge>>(
                    insertion_edge_vector);
            auto removed_edge_set = container_copy::to_unordered_set<shared_ptr<abstract_bipartite_edge>>(
                    removal_edge_vector);
            for (const auto &[l, l_node]: *previous_left_index_map) {
                branch_left_index_map->insert({l, make_shared<bipartite_core_left_store_index>(l_node)});
            }
            for (const auto &[r, r_node]: *previous_right_index_map) {
                branch_right_index_map->insert({r, make_shared<bipartite_core_right_store_index>(r_node)});
            }

            simple_timer maintenance_timer;
            {
                auto left_core_degree_map = make_shared<unordered_map<uint32_t, shared_ptr<map<uint32_t, shared_ptr<unordered_set<uint32_t>>>>>>();
                auto right_core_degree_map = make_shared<unordered_map<uint32_t, shared_ptr<map<uint32_t, shared_ptr<unordered_set<uint32_t>>>>>>();

                auto new_left_index_map = make_shared<unordered_map<uint32_t, shared_ptr<bipartite_core_left_store_index>>>();
                auto new_right_index_map = make_shared<unordered_map<uint32_t, shared_ptr<bipartite_core_right_store_index>>>();

                auto pool = make_shared<thread_pool>(thread_number);

                branch_bipartite_core_maintenance::init(B, branch_left_index_map, branch_right_index_map, left_core_degree_map, right_core_degree_map,
                                                        new_left_index_map, new_right_index_map, pool);
                branch_bipartite_core_maintenance::remove(B,  left_core_degree_map, right_core_degree_map, removed_edge_set, branch_left_index_map, branch_right_index_map,
                                                          new_left_index_map, new_right_index_map, pool);
                branch_bipartite_core_maintenance::insert(B, left_core_degree_map, right_core_degree_map, inserted_edge_set, branch_left_index_map, branch_right_index_map,
                                                          new_left_index_map, new_right_index_map, pool);
            }
            auto maintenance_time = maintenance_timer.get_elapse_second();

            if (bipartite_core_compare::same(branch_left_index_map, branch_right_index_map, contrastive_left_index_map, contrastive_right_index_map)) {
                LOG(logger, LOG_RANK::INFO) << "Branch Maintenance*," << maintenance_time << "\n";
            }

            B->insert_edge_collection(removal_edge_vector);
            B->remove_edge_collection(insertion_edge_vector);
        }
    }
    return 0;
}


