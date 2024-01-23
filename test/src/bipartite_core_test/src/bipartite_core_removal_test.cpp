#include "bipartite_core_test/bipartite_core_test.h"

int main(int argc, char** argv){
    if(argc < 4){
        std::cout<<"Usage: please input path, file name, and thread number!";
    }
    string path = argv[1];
    string input_file_name = argv[2];
    uint32_t thread_number = std::stoul(argv[3]);

    string log_path = path + "log/";

    auto directory = std::filesystem::path(log_path);
    if (!exists(directory)) {
        create_directory(directory);
    }

    string log_file_name = log_path  + "removal_test.log";
    scnu::simple_logger logger(log_file_name);

    LOG(logger, LOG_RANK::INFO) << input_file_name << "\n";

    uint32_t m = 500;// the number of removed edges

    auto B = shared_ptr<abstract_bipartite_graph>();

    auto removal_edge_vector = make_shared<vector<shared_ptr<abstract_bipartite_edge>>>(m);
    {
        auto edge_vector = abstract_bipartite_graph_io::get_edge_vector(path, input_file_name);
        {
//            auto rd = make_shared<random_device>();
//            shuffle(edge_vector->begin(), edge_vector->end(), *random_generator::get_default_engine(rd));
//            for(const auto &e:*edge_vector){
//                printf("%u,%u\n",e->get_left_vertex_id(), e->get_right_vertex_id());
//            }
        }
        for (auto i = 0; i < m; i++) {
            removal_edge_vector->at(i) = edge_vector->at(i);
        }

        auto pool = make_shared<thread_pool>(thread_number);
        B = abstract_bipartite_graph_io::load_graph(edge_vector, pool);
    }

    auto previous_left_index_map = make_shared<unordered_map<uint32_t, shared_ptr<bipartite_core_left_store_index>>>();
    auto previous_right_index_map = make_shared<unordered_map<uint32_t, shared_ptr<bipartite_core_right_store_index>>>();

    auto previous_core_order_index = make_shared<bipartite_core_order_index>();
    auto previous_core_rem_degree_index = make_shared<bipartite_core_rem_degree_index>();
    auto previous_core_degree_index = make_shared<bipartite_core_degree_index>();

    auto previous_delta = make_shared<uint32_t>(0);
    {
        auto pool = make_shared<thread_pool>(thread_number);

        auto left_mutex_map = make_shared<unordered_map<uint32_t, shared_ptr<mutex>>>();
        auto right_mutex_map = make_shared<unordered_map<uint32_t, shared_ptr<mutex>>>();

        branch_bipartite_core_decomposition::init(B, left_mutex_map, right_mutex_map, previous_left_index_map,
                                                  previous_right_index_map, pool);
        *previous_delta = branch_bipartite_core_decomposition::decompose(B,
                                                                         left_mutex_map,
                                                                         right_mutex_map,
                                                                         previous_left_index_map,
                                                                         previous_right_index_map,
                                                                         previous_core_order_index,
                                                                         previous_core_rem_degree_index,
                                                                         previous_core_degree_index,
                                                                         pool);
    }

    /**
      * directly decompose the remain graph
      */
    auto contrastive_left_index_map = make_shared<unordered_map<uint32_t, shared_ptr<bipartite_core_left_store_index>>>();
    auto contrastive_right_index_map = make_shared<unordered_map<uint32_t, shared_ptr<bipartite_core_right_store_index>>>();
    {
        B->remove_edge_collection(removal_edge_vector);

        simple_timer decomposition_timer;
        auto pool = make_shared<thread_pool>(thread_number);

        auto left_mutex_map = make_shared<unordered_map<uint32_t, shared_ptr<mutex>>>();
        auto right_mutex_map = make_shared<unordered_map<uint32_t, shared_ptr<mutex>>>();

        share_bipartite_core_decomposition::init(B, left_mutex_map, right_mutex_map, contrastive_left_index_map,
                                                 contrastive_right_index_map, pool);
        share_bipartite_core_decomposition::decompose(B, left_mutex_map, right_mutex_map, contrastive_left_index_map,
                                                      contrastive_right_index_map, pool);
        auto decomposition_time = decomposition_timer.get_elapse_second();

        LOG(logger, LOG_RANK::INFO) << "Decomposition," << decomposition_time << ","
                                    << double(process_information::get_memory(getpid())) / 1024 / 1024 << "\n";

        B->insert_edge_collection(removal_edge_vector);
    }

    {
        auto edge_left_index_map = make_shared<unordered_map<uint32_t,shared_ptr<bipartite_core_left_store_index>>>();
        for(const auto& [l,l_node]:*previous_left_index_map){
            edge_left_index_map->insert({l,make_shared<bipartite_core_left_store_index>(l_node)});
        }
        auto edge_right_index_map = make_shared<unordered_map<uint32_t,shared_ptr<bipartite_core_right_store_index>>>();
        for(const auto& [r,r_node]:*previous_right_index_map)
        {
            edge_right_index_map->insert({r,make_shared<bipartite_core_right_store_index>(r_node)});
        }
        auto delta = make_shared<uint32_t>(*previous_delta);

        simple_timer maintenance_timer;
        {
            auto new_left_index_map = make_shared<unordered_map<uint32_t, shared_ptr<bipartite_core_left_store_index>>>();
            auto new_right_index_map = make_shared<unordered_map<uint32_t, shared_ptr<bipartite_core_right_store_index>>>();


            auto left_mutex_map = make_shared<unordered_map<uint32_t, shared_ptr<mutex>>>();
            auto right_mutex_map = make_shared<unordered_map<uint32_t, shared_ptr<mutex>>>();

            auto pool = make_shared<thread_pool>(thread_number);

            edge_bipartite_core_maintenance::init(B, left_mutex_map, right_mutex_map, new_left_index_map,
                                                  new_right_index_map, pool);
            for (const auto &e: *removal_edge_vector) {
                edge_bipartite_core_maintenance::remove(B, left_mutex_map, right_mutex_map, e, edge_left_index_map,
                                                        edge_right_index_map, new_left_index_map, new_right_index_map,
                                                        delta, pool);
            }
        }
        auto maintenance_time = maintenance_timer.get_elapse_second();

        if(bipartite_core_compare::same(edge_left_index_map,edge_right_index_map,contrastive_left_index_map,contrastive_right_index_map)){
            LOG(logger, LOG_RANK::INFO) << "Edge Removal," << maintenance_time  <<"," << double(process_information::get_memory(getpid())) / 1024 / 1024<< "\n";
        }

        B->insert_edge_collection(removal_edge_vector);
    }

    {
        auto batch_left_index_map = make_shared<unordered_map<uint32_t, shared_ptr<bipartite_core_left_store_index>>>();
        for (const auto &[l, l_node]: *previous_left_index_map) {
            batch_left_index_map->insert({l, make_shared<bipartite_core_left_store_index>(l_node)});
        }
        auto batch_right_index_map = make_shared<unordered_map<uint32_t, shared_ptr<bipartite_core_right_store_index>>>();
        for (const auto &[r, r_node]: *previous_right_index_map) {
            batch_right_index_map->insert({r, make_shared<bipartite_core_right_store_index>(r_node)});
        }
        auto removal_edge_set = container_copy::to_unordered_set<shared_ptr<abstract_bipartite_edge>>(
                removal_edge_vector);

        auto delta = make_shared<uint32_t>(*previous_delta);

        simple_timer maintenance_timer;
        {
            auto new_left_index_map = make_shared<unordered_map<uint32_t, shared_ptr<bipartite_core_left_store_index>>>();
            auto new_right_index_map = make_shared<unordered_map<uint32_t, shared_ptr<bipartite_core_right_store_index>>>();

            auto left_mutex_map = make_shared<unordered_map<uint32_t, shared_ptr<mutex>>>();
            auto right_mutex_map = make_shared<unordered_map<uint32_t, shared_ptr<mutex>>>();

            auto pool = make_shared<thread_pool>(thread_number);

            edge_bipartite_core_maintenance::init(B, left_mutex_map, right_mutex_map, new_left_index_map,
                                                  new_right_index_map, pool);
            edge_bipartite_core_maintenance::batch_remove(B, left_mutex_map, right_mutex_map, removal_edge_set,
                                                          batch_left_index_map,
                                                          batch_right_index_map, new_left_index_map,
                                                          new_right_index_map, delta, pool);
        }
        auto maintenance_time = maintenance_timer.get_elapse_second();

        if (bipartite_core_compare::same(batch_left_index_map, batch_right_index_map, contrastive_left_index_map,
                                         contrastive_right_index_map)) {
            LOG(logger, LOG_RANK::INFO) << "Batch Removal," << maintenance_time  <<"," << double(process_information::get_memory(getpid())) / 1024 / 1024 << "\n";
        }
        B->insert_edge_collection(removal_edge_vector);
    }

    /**
     * @brief quasi maintenance
     */
     {
        auto quasi_left_vertex_map = make_shared<unordered_map<uint32_t,shared_ptr<bipartite_core_left_store_index>>>();
        for(const auto& [l,l_node]:*previous_left_index_map){
            quasi_left_vertex_map->insert({l, make_shared<bipartite_core_left_store_index>(l_node)});
        }
        auto quasi_right_vertex_map = make_shared<unordered_map<uint32_t,shared_ptr<bipartite_core_right_store_index>>>();
        for(const auto& [r,r_node]:*previous_right_index_map)
        {
            quasi_right_vertex_map->insert({r, make_shared<bipartite_core_right_store_index>(r_node)});
        }
        auto removal_edge_set = container_copy::to_unordered_set<shared_ptr<abstract_bipartite_edge>>(removal_edge_vector);

        simple_timer maintenance_timer;
        {
            auto new_left_index_map = make_shared<unordered_map<uint32_t, shared_ptr<bipartite_core_left_store_index>>>();
            auto new_right_index_map = make_shared<unordered_map<uint32_t, shared_ptr<bipartite_core_right_store_index>>>();

            auto pool = make_shared<thread_pool>(thread_number);
            quasi_bipartite_core_maintenance::init(B, new_left_index_map, new_right_index_map, pool);
            quasi_bipartite_core_maintenance::remove(B, removal_edge_set, quasi_left_vertex_map,
                                                     quasi_right_vertex_map, new_left_index_map, new_right_index_map, pool);
        }
        auto maintenance_time = maintenance_timer.get_elapse_second();

        if(bipartite_core_compare::same(quasi_left_vertex_map,quasi_right_vertex_map,contrastive_left_index_map,contrastive_right_index_map)){
            LOG(logger, LOG_RANK::INFO) << "Quasi Removal," << maintenance_time  <<"," << double(process_information::get_memory(getpid())) / 1024 / 1024 << "\n";
        }

        B->insert_edge_collection(removal_edge_vector);
    }


    /**
     * @brief branch remove maintenance
     */
    {
        auto branch_left_vertex_map = make_shared<unordered_map<uint32_t, shared_ptr<bipartite_core_left_store_index>>>();
        for (const auto &[l, l_node]: *previous_left_index_map) {
            branch_left_vertex_map->insert({l, make_shared<bipartite_core_left_store_index>(l_node)});
        }
        auto branch_right_vertex_map = make_shared<unordered_map<uint32_t, shared_ptr<bipartite_core_right_store_index>>>();
        for (const auto &[r, r_node]: *previous_right_index_map) {
            branch_right_vertex_map->insert({r, make_shared<bipartite_core_right_store_index>(r_node)});
        }

        auto removal_edge_set = container_copy::to_unordered_set<shared_ptr<abstract_bipartite_edge>>(
                removal_edge_vector);

        auto branch_core_order_index = make_shared<bipartite_core_order_index>(previous_core_order_index);
        auto branch_core_rem_degree_index = make_shared<bipartite_core_rem_degree_index>(
                previous_core_rem_degree_index);
        auto branch_core_degree_index = make_shared<bipartite_core_degree_index>(previous_core_degree_index);

        simple_timer maintenance_timer;
        {

            auto new_left_index_map = make_shared<unordered_map<uint32_t, shared_ptr<bipartite_core_left_store_index>>>();
            auto new_right_index_map = make_shared<unordered_map<uint32_t, shared_ptr<bipartite_core_right_store_index>>>();

            auto pool = make_shared<thread_pool>(thread_number);

            auto left_mutex_map = make_shared<unordered_map<uint32_t, shared_ptr<mutex>>>();
            auto right_mutex_map = make_shared<unordered_map<uint32_t, shared_ptr<mutex>>>();

            branch_bipartite_core_maintenance::init(B, left_mutex_map, right_mutex_map, new_left_index_map,
                                                    new_right_index_map, pool);
            branch_bipartite_core_maintenance::remove(B, left_mutex_map, right_mutex_map, removal_edge_set,
                                                      branch_left_vertex_map,
                                                      branch_right_vertex_map,
                                                      new_left_index_map, new_right_index_map,
                                                      branch_core_order_index,
                                                      branch_core_rem_degree_index,
                                                      branch_core_degree_index,
                                                      pool);
        }
        auto maintenance_time = maintenance_timer.get_elapse_second();

        if (bipartite_core_compare::same(branch_left_vertex_map, branch_right_vertex_map, contrastive_left_index_map,
                                         contrastive_right_index_map)) {
            LOG(logger, LOG_RANK::INFO) << "Branch Removal," << maintenance_time << ","
                                        << double(process_information::get_memory(getpid())) / 1024 / 1024 << "\n";
        }

        B->insert_edge_collection(removal_edge_vector);
    }

    return 0;
}
