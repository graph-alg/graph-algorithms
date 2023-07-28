#include "bipartite_core_test/bipartite_core_test.h"

int main(int argc, char **argv){
    if(argc < 2){
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

    string log_file_name = log_path + "insertion_test.log";
    simple_logger logger(log_file_name);

    LOG(logger, LOG_RANK::INFO) << input_file_name << "\n";

    auto edge_vector = abstract_bipartite_graph_io::get_edge_vector(path, input_file_name);

    uint32_t m = 2;

    auto rd = make_shared<random_device>();
    //shuffle(edge_vector->begin(), edge_vector->end(), *random_generator::get_default_engine(rd));

    while (true) {
//        auto e = *edge_vector->begin();
//        edge_vector->erase(edge_vector->begin());

        auto insertion_edge_vector = make_shared<vector<shared_ptr<abstract_bipartite_edge>>>();
        auto B = shared_ptr<abstract_bipartite_graph>();
        {

            {


//            for (const auto &e: *edge_vector) {
//                printf("%u,%u\n", e->get_left_vertex_id(), e->get_right_vertex_id());
//            }
            }

            auto previous_edge_vector = make_shared<vector<shared_ptr<abstract_bipartite_edge>>>();
            for (uint32_t i = 0; i < edge_vector->size(); i++) {
                if (i < uint32_t(edge_vector->size()) - m) {
                    previous_edge_vector->push_back(edge_vector->at(i));
                } else {
                    insertion_edge_vector->push_back(edge_vector->at(i));
                }
            }

            auto pool = make_shared<thread_pool>(thread_number);
            B = abstract_bipartite_graph_io::load_graph(previous_edge_vector, pool);
        }

        auto previous_left_index_map = make_shared<unordered_map<uint32_t, shared_ptr<bipartite_core_left_store_index>>>();
        auto previous_right_index_map = make_shared<unordered_map<uint32_t, shared_ptr<bipartite_core_right_store_index>>>();
        auto previous_core_order_index = make_shared<bipartite_core_order_index>();
        auto previous_core_rem_degree_index = make_shared<bipartite_core_rem_degree_index>();
        auto previous_core_degree_index = make_shared<bipartite_core_degree_index>();

        auto previous_delta = make_shared<uint32_t>(0);
        {
            auto pool = make_shared<thread_pool>(thread_number);
            *previous_delta = branch_bipartite_core_decomposition::decompose(B,
                                                                             previous_left_index_map,
                                                                             previous_right_index_map,
                                                                             previous_core_order_index,
                                                                             previous_core_rem_degree_index,
                                                                             previous_core_degree_index,
                                                                             pool);
        }

        /**
          * directly decompose the graph
        */
        auto contrastive_left_index_map = make_shared<unordered_map<uint32_t, shared_ptr<bipartite_core_left_store_index>>>();
        auto contrastive_right_index_map = make_shared<unordered_map<uint32_t, shared_ptr<bipartite_core_right_store_index>>>();
        {

            B->insert_edge_collection(insertion_edge_vector);
            simple_timer decomposition_timer;
            auto pool = make_shared<thread_pool>(thread_number);
            share_bipartite_core_decomposition::decompose(B, contrastive_left_index_map, contrastive_right_index_map,
                                                          pool);
            double decomposition_time = decomposition_timer.get_elapse_second();

            LOG(logger, LOG_RANK::INFO) << "Decomposition," << decomposition_time << "\n";

            B->remove_edge_collection(insertion_edge_vector);
        }

        /**
        * @brief edge insert
        */
        {
            auto edge_left_vertex_map = make_shared<unordered_map<uint32_t, shared_ptr<bipartite_core_left_store_index>>>();
            auto edge_right_vertex_map = make_shared<unordered_map<uint32_t, shared_ptr<bipartite_core_right_store_index>>>();

            for (const auto &[l, l_node]: *previous_left_index_map) {
                edge_left_vertex_map->insert({l, make_shared<bipartite_core_left_store_index>(l_node)});
            }
            for (const auto &[r, r_node]: *previous_right_index_map) {
                edge_right_vertex_map->insert({r, make_shared<bipartite_core_right_store_index>(r_node)});
            }
            auto delta = make_shared<uint32_t>(*previous_delta);


            simple_timer maintenance_timer;
            {

                auto new_left_index_map = make_shared<unordered_map<uint32_t, shared_ptr<bipartite_core_left_store_index>>>();
                auto new_right_index_map = make_shared<unordered_map<uint32_t, shared_ptr<bipartite_core_right_store_index>>>();

                auto pool = make_shared<thread_pool>(thread_number);
                edge_bipartite_core_maintenance::init(B, new_left_index_map, new_right_index_map, pool);
                for (const auto &e: *insertion_edge_vector) {
                    edge_bipartite_core_maintenance::insert(B, e, edge_left_vertex_map, edge_right_vertex_map,
                                                            new_left_index_map, new_right_index_map, delta,
                                                            pool);
                }
            }
            double maintenance_time = maintenance_timer.get_elapse_second();


            if (bipartite_core_compare::same(edge_left_vertex_map, edge_right_vertex_map, contrastive_left_index_map,
                                             contrastive_right_index_map)) {
                LOG(logger, LOG_RANK::INFO) << "Edge Insertion," << maintenance_time << "\n";
            }

            B->remove_edge_collection(insertion_edge_vector);
        }

        /**
         * @brief batch insert
         */
        {
            auto edge_left_vertex_map = make_shared<unordered_map<uint32_t, shared_ptr<bipartite_core_left_store_index>>>();
            auto edge_right_vertex_map = make_shared<unordered_map<uint32_t, shared_ptr<bipartite_core_right_store_index>>>();

            auto inserted_edge_set = container_copy::to_unordered_set<shared_ptr<abstract_bipartite_edge>>(
                    insertion_edge_vector);

            for (const auto &[l, l_node]: *previous_left_index_map) {
                edge_left_vertex_map->insert({l, make_shared<bipartite_core_left_store_index>(l_node)});
            }
            for (const auto &[r, r_node]: *previous_right_index_map) {
                edge_right_vertex_map->insert({r, make_shared<bipartite_core_right_store_index>(r_node)});
            }
            auto delta = make_shared<uint32_t>(*previous_delta);

            simple_timer maintenance_timer;
            {
                auto new_left_index_map = make_shared<unordered_map<uint32_t, shared_ptr<bipartite_core_left_store_index>>>();
                auto new_right_index_map = make_shared<unordered_map<uint32_t, shared_ptr<bipartite_core_right_store_index>>>();

                auto pool = make_shared<thread_pool>(thread_number);

                edge_bipartite_core_maintenance::init(B, new_left_index_map, new_right_index_map, pool);
                edge_bipartite_core_maintenance::batch_insert(B, inserted_edge_set, edge_left_vertex_map,
                                                              edge_right_vertex_map,
                                                              new_left_index_map, new_right_index_map, delta, pool);
            }
            double maintenance_time = maintenance_timer.get_elapse_second();

            if (bipartite_core_compare::same(edge_left_vertex_map, edge_right_vertex_map, contrastive_left_index_map,
                                             contrastive_right_index_map)) {
                LOG(logger, LOG_RANK::INFO) << "Batch Insertion," << maintenance_time << "\n";
            }

            B->remove_edge_collection(insertion_edge_vector);
        }

        {
            auto quasi_left_vertex_map = make_shared<unordered_map<uint32_t, shared_ptr<bipartite_core_left_store_index>>>();
            auto quasi_right_vertex_map = make_shared<unordered_map<uint32_t, shared_ptr<bipartite_core_right_store_index>>>();

            auto inserted_edge_set = container_copy::to_unordered_set<shared_ptr<abstract_bipartite_edge>>(
                    insertion_edge_vector);
            for (const auto &[l, l_node]: *previous_left_index_map) {
                quasi_left_vertex_map->insert({l, make_shared<bipartite_core_left_store_index>(l_node)});
            }
            for (const auto &[r, r_node]: *previous_right_index_map) {
                quasi_right_vertex_map->insert({r, make_shared<bipartite_core_right_store_index>(r_node)});
            }

            simple_timer maintenance_timer;
            {
                auto new_left_index_map = make_shared<unordered_map<uint32_t, shared_ptr<bipartite_core_left_store_index>>>();
                auto new_right_index_map = make_shared<unordered_map<uint32_t, shared_ptr<bipartite_core_right_store_index>>>();

                auto pool = make_shared<thread_pool>(thread_number);
                quasi_bipartite_core_maintenance::init(B, new_left_index_map, new_right_index_map, pool);
                quasi_bipartite_core_maintenance::insert(B, inserted_edge_set, quasi_left_vertex_map,
                                                         quasi_right_vertex_map, new_left_index_map,
                                                         new_right_index_map, pool);
            }

            auto maintenance_time = maintenance_timer.get_elapse_second();

            if (bipartite_core_compare::same(quasi_left_vertex_map, quasi_right_vertex_map, contrastive_left_index_map,
                                             contrastive_right_index_map)) {
                LOG(logger, LOG_RANK::INFO) << "Quasi Insertion," << maintenance_time << "\n";
            }

            B->remove_edge_collection(insertion_edge_vector);
        }

        /**
         * @brief branch maintenance without index
         */
        {
            auto branch_left_vertex_map = make_shared<unordered_map<uint32_t, shared_ptr<bipartite_core_left_store_index>>>();
            auto branch_right_vertex_map = make_shared<unordered_map<uint32_t, shared_ptr<bipartite_core_right_store_index>>>();
            auto inserted_edge_set = container_copy::to_unordered_set<shared_ptr<abstract_bipartite_edge>>(
                    insertion_edge_vector);
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
                branch_bipartite_core_maintenance::insert(B, inserted_edge_set, branch_left_vertex_map,
                                                          branch_right_vertex_map,
                                                          new_left_index_map, new_right_index_map,
                                                          pool);
            }
            auto maintenance_time = maintenance_timer.get_elapse_second();

            if (bipartite_core_compare::same(branch_left_vertex_map, branch_right_vertex_map,
                                             contrastive_left_index_map,
                                             contrastive_right_index_map)) {
                LOG(logger, LOG_RANK::INFO) << "Branch Insertion," << maintenance_time << "\n";
            }

            B->remove_edge_collection(insertion_edge_vector);
        }

        /**
         * @brief branch maintenance with index
         */
        {
            auto branch_left_vertex_map = make_shared<unordered_map<uint32_t, shared_ptr<bipartite_core_left_store_index>>>();
            auto branch_right_vertex_map = make_shared<unordered_map<uint32_t, shared_ptr<bipartite_core_right_store_index>>>();
            auto inserted_edge_set = container_copy::to_unordered_set<shared_ptr<abstract_bipartite_edge>>(
                    insertion_edge_vector);
            for (const auto &[l, l_node]: *previous_left_index_map) {
                branch_left_vertex_map->insert({l, make_shared<bipartite_core_left_store_index>(l_node)});
            }
            for (const auto &[r, r_node]: *previous_right_index_map) {
                branch_right_vertex_map->insert({r, make_shared<bipartite_core_right_store_index>(r_node)});
            }

            auto branch_core_order_index = make_shared<bipartite_core_order_index>(previous_core_order_index);
            auto branch_core_rem_degree_index = make_shared<bipartite_core_rem_degree_index>(
                    previous_core_rem_degree_index);
            auto branch_core_degree_index = make_shared<bipartite_core_degree_index>(previous_core_degree_index);

            simple_timer maintenance_timer;
            {
                auto new_left_index_map = make_shared<unordered_map<uint32_t, shared_ptr<bipartite_core_left_store_index>>>();
                auto new_right_index_map = make_shared<unordered_map<uint32_t, shared_ptr<bipartite_core_right_store_index>>>();

                auto pool = make_shared<thread_pool>(thread_number);
                branch_bipartite_core_maintenance::init(B, new_left_index_map, new_right_index_map, pool);
                branch_bipartite_core_maintenance::insert(B, inserted_edge_set,
                                                          branch_core_order_index,
                                                          branch_core_rem_degree_index,
                                                          branch_core_degree_index,
                                                          branch_left_vertex_map,
                                                          branch_right_vertex_map,
                                                          new_left_index_map, new_right_index_map,
                                                          pool);
            }
            auto maintenance_time = maintenance_timer.get_elapse_second();

            if (bipartite_core_compare::same(branch_left_vertex_map, branch_right_vertex_map,
                                             contrastive_left_index_map,
                                             contrastive_right_index_map)) {
                LOG(logger, LOG_RANK::INFO) << "Branch Insertion2," << maintenance_time << "\n";
                //edge_vector->push_back(e);
                shuffle(edge_vector->begin(), edge_vector->end(), *random_generator::get_default_engine(rd));
            } else {
                printf("%lu\n", edge_vector->size());
                if (edge_vector->size() < 21) {
                    for (const auto &e0: *edge_vector) {
                        printf("%u,%u\n", e0->get_left_vertex_id(), e0->get_right_vertex_id());
                    }
                    break;
                }

            }

            B->remove_edge_collection(insertion_edge_vector);
        }
    }


    return 0;
}
