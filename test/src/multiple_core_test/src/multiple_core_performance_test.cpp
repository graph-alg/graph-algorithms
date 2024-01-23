
#include "multiple_core_test/multiple_core_test.h"

void decomposition_test(const shared_ptr<temporal_graph>& G,
                        uint32_t thread_number,
                        const string& log_path,
                        const string& input_file_name){
    const string log_file_name = log_path + "decomposition_test.log";
    scnu::simple_logger logger(log_file_name);

    LOG(logger,LOG_RANK::INFO) << input_file_name << "\n";

    auto vertex_base_map = make_shared<unordered_map<uint32_t,shared_ptr<multiple_core_pair_map_index>>>();
    {
        simple_timer decomposition_timer;
        auto pool = make_shared<thread_pool>(thread_number);
        auto vertex_mutex_map = make_shared<unordered_map<uint32_t, shared_ptr<mutex>>>();
        basic_multiple_core_decomposition::init(G, vertex_mutex_map, vertex_base_map, pool);
        basic_multiple_core_decomposition::decompose(G, vertex_mutex_map, vertex_base_map, pool);
        auto decomposition_time = decomposition_timer.get_elapse_second();

        LOG(logger, LOG_RANK::INFO) << "Based Decomposition2," << decomposition_time << "\n";
    }

    {
        auto vertex_hierarchical_map = make_shared<unordered_map<uint32_t,shared_ptr<multiple_core_pair_map_index>>>();

        simple_timer decomposition_timer;
        auto vertex_mutex_map = make_shared<unordered_map<uint32_t, shared_ptr<mutex>>>();
        auto pool = make_shared<thread_pool>(thread_number);
        hierarchy_multiple_core_decomposition::init(G, vertex_mutex_map, vertex_hierarchical_map, pool);
        hierarchy_multiple_core_decomposition::decompose(G, vertex_mutex_map, vertex_hierarchical_map, pool);
        auto decomposition_time = decomposition_timer.get_elapse_second();

        if(multiple_core_compare::same(vertex_base_map, vertex_hierarchical_map)) {
            LOG(logger, LOG_RANK::INFO) << "Hierarchical Decomposition2," << decomposition_time << "\n";
        }
    }

    {
        auto vertex_branch_map = make_shared<unordered_map<uint32_t,shared_ptr<multiple_core_pair_map_index>>>();

        simple_timer decomposition_timer;
        auto vertex_mutex_map = make_shared<unordered_map<uint32_t, shared_ptr<mutex>>>();
        auto pool = make_shared<thread_pool>(thread_number);
        branch_multiple_core_decomposition::init(G, vertex_mutex_map, vertex_branch_map, pool);
        branch_multiple_core_decomposition::decompose(G, vertex_mutex_map, vertex_branch_map, pool);
        auto decomposition_time = decomposition_timer.get_elapse_second();

        if(multiple_core_compare::same(vertex_base_map, vertex_branch_map)){
            LOG(logger,LOG_RANK::INFO) << "Branch Decomposition2," << decomposition_time << "\n";
        }
    }

    {
        auto index_vertex_branch_map = make_shared<unordered_map<uint32_t,shared_ptr<multiple_core_pair_map_index>>>();

        simple_timer decomposition_timer;
        auto vertex_edge_size_map = make_shared<unordered_map<uint32_t,
                shared_ptr<map<uint32_t ,shared_ptr<unordered_set<uint32_t>>>>>>();
        auto vertex_mutex_map = make_shared<unordered_map<uint32_t, shared_ptr<mutex>>>();
        auto pool = make_shared<thread_pool>(thread_number);
        branch_multiple_core_decomposition::init(G, vertex_mutex_map, vertex_edge_size_map,index_vertex_branch_map, pool);
        branch_multiple_core_decomposition::decompose(G, vertex_mutex_map, vertex_edge_size_map,
                                                      index_vertex_branch_map, pool);
        auto decomposition_time = decomposition_timer.get_elapse_second();

        if(multiple_core_compare::same(vertex_base_map, index_vertex_branch_map)){
            LOG(logger,LOG_RANK::INFO) << "Branch Decomposition2*," << decomposition_time << "\n";
        }
    }
}

void insertion_test(const shared_ptr<temporal_graph>& G,
                    const shared_ptr<vector<shared_ptr<temporal_edge>>>& edge_vector,
                    uint32_t thread_number,
                    const string& log_path,
                    const string& input_file_name){
    const string log_file_name = log_path + "insertion_test.log";
    scnu::simple_logger logger(log_file_name);


    LOG(logger, LOG_RANK::INFO) << input_file_name << "\n";

    const uint32_t m = 5000;
    auto previous_edge_vector = make_shared<vector<shared_ptr<temporal_edge>>>();
    auto insertion_edge_vector = make_shared<vector<shared_ptr<temporal_edge>>>();
    for (uint32_t i = 0; i < edge_vector->size(); ++i) {
        if (i < edge_vector->size() - m) {
            previous_edge_vector->push_back(edge_vector->at(i));
        } else {
            insertion_edge_vector->push_back(edge_vector->at(i));
        }
    }

    auto previous_vertex_index_map = make_shared<unordered_map<uint32_t, shared_ptr<multiple_core_pair_map_index>>>();
    {
        G->remove_edge_collection(insertion_edge_vector);
        auto pool = make_shared<thread_pool>(thread_number);

        auto vertex_edge_index_map = make_shared<unordered_map<uint32_t, shared_ptr<map<uint32_t, shared_ptr<unordered_set<uint32_t>>>>>>();
        branch_multiple_core_decomposition::init(G, vertex_edge_index_map, previous_vertex_index_map, pool);
        branch_multiple_core_decomposition::decompose(G, vertex_edge_index_map, previous_vertex_index_map, pool);
    }

    auto base_map = make_shared<unordered_map<uint32_t, shared_ptr<multiple_core_pair_map_index>>>();
    {
        G->insert_edge_collection(insertion_edge_vector);

        simple_timer decomposition_timer;
        auto pool = make_shared<thread_pool>(thread_number);
        auto vertex_mutex_map = make_shared<unordered_map<uint32_t, shared_ptr<mutex>>>();
        branch_multiple_core_decomposition::init(G, vertex_mutex_map, base_map, pool);
        branch_multiple_core_decomposition::decompose(G, vertex_mutex_map, base_map, pool);
        auto decomposition_time = decomposition_timer.get_elapse_second();

        LOG(logger, LOG_RANK::INFO) << "Branch Decomposition2," << decomposition_time << "\n";

        G->remove_edge_collection(insertion_edge_vector);
    }


    {
        auto base_map2 = make_shared<unordered_map<uint32_t, shared_ptr<multiple_core_pair_map_index>>>();
        G->insert_edge_collection(insertion_edge_vector);

        simple_timer decomposition_timer;
        auto pool = make_shared<thread_pool>(thread_number);
        auto vertex_mutex_map = make_shared<unordered_map<uint32_t, shared_ptr<mutex>>>();
        auto vertex_edge_index_map = make_shared<unordered_map<uint32_t, shared_ptr<map<uint32_t, shared_ptr<unordered_set<uint32_t>>>>>>();
        branch_multiple_core_decomposition::init(G, vertex_mutex_map, vertex_edge_index_map, base_map2, pool);
        branch_multiple_core_decomposition::decompose(G, vertex_mutex_map, vertex_edge_index_map, base_map2, pool);
        auto decomposition_time = decomposition_timer.get_elapse_second();

        LOG(logger, LOG_RANK::INFO) << "Branch Decomposition2*," << decomposition_time << "\n";

        G->remove_edge_collection(insertion_edge_vector);
    }

    {
        auto vertex_index_map = make_shared<unordered_map<uint32_t, shared_ptr<multiple_core_pair_map_index>>>();
        for (const auto &[u, u_index]: *previous_vertex_index_map) {
            vertex_index_map->insert({u, make_shared<multiple_core_pair_map_index>(u_index)});
        }

        auto insertion_edge_set = container_copy::to_unordered_set<shared_ptr<temporal_edge>>(
                insertion_edge_vector);

        simple_timer maintenance_timer;
        auto pool = make_shared<thread_pool>(thread_number);
        auto vertex_mutex_map = make_shared<unordered_map<uint32_t, shared_ptr<mutex>>>();
        auto vertex_degree_map = make_shared<unordered_map<uint32_t, uint32_t>>();

        hierarchy_multiple_core_maintenance::init(G, vertex_mutex_map, vertex_degree_map,
                                                  pool);
        hierarchy_multiple_core_maintenance::insert(G, vertex_mutex_map, insertion_edge_set, vertex_index_map, vertex_degree_map, pool);

        auto maintenance_time = maintenance_timer.get_elapse_second();

        if (multiple_core_compare::same(base_map, vertex_index_map)) {
            LOG(logger, LOG_RANK::INFO) << "Hierarchy Insertion2," << maintenance_time << "\n";
        }

        G->remove_edge_collection(insertion_edge_vector);
    }

    {
        auto vertex_index_map = make_shared<unordered_map<uint32_t, shared_ptr<multiple_core_pair_map_index>>>();
        for (const auto &[u, u_index]: *previous_vertex_index_map) {
            vertex_index_map->insert({u, make_shared<multiple_core_pair_map_index>(u_index)});
        }

        auto insertion_edge_set = container_convert::to_unordered_set<shared_ptr<temporal_edge>>(
                insertion_edge_vector);

        simple_timer maintenance_timer;
        auto pool = make_shared<thread_pool>(thread_number);
        auto vertex_mutex_map = make_shared<unordered_map<uint32_t, shared_ptr<mutex>>>();
        auto vertex_degree_map = make_shared<unordered_map<uint32_t, uint32_t>>();

        branch_multiple_core_maintenance::init(G, vertex_mutex_map, vertex_degree_map, pool);
        branch_multiple_core_maintenance::insert(G, vertex_mutex_map, insertion_edge_set, vertex_index_map, vertex_degree_map, pool);

        auto maintenance_time = maintenance_timer.get_elapse_second();

        if (multiple_core_compare::same(base_map, vertex_index_map)) {
            LOG(logger, LOG_RANK::INFO) << "Branch Insertion2," << maintenance_time << "\n";
        }

        G->remove_edge_collection(insertion_edge_vector);
    }

    {
        auto vertex_index_map = make_shared<unordered_map<uint32_t, shared_ptr<multiple_core_pair_map_index>>>();
        for (const auto &[u, u_index]: *previous_vertex_index_map) {
            vertex_index_map->insert({u, make_shared<multiple_core_pair_map_index>(u_index)});
        }

        auto insertion_edge_set = container_convert::to_unordered_set<shared_ptr<temporal_edge>>(
                insertion_edge_vector);

        simple_timer maintenance_timer;
        auto pool = make_shared<thread_pool>(thread_number);
        auto vertex_edge_index_map = make_shared<unordered_map<uint32_t, shared_ptr<map<uint32_t, shared_ptr<unordered_set<uint32_t>>>>>>();
        auto vertex_mutex_map = make_shared<unordered_map<uint32_t, shared_ptr<mutex>>>();
        auto vertex_degree_map = make_shared<unordered_map<uint32_t, uint32_t>>();

        branch_multiple_core_maintenance::init(G, vertex_mutex_map, vertex_edge_index_map,
                                               vertex_degree_map, pool);
        branch_multiple_core_maintenance::insert(G, vertex_mutex_map, insertion_edge_set, vertex_edge_index_map,
                                                 vertex_index_map, vertex_degree_map, pool);

        auto maintenance_time = maintenance_timer.get_elapse_second();

        if (multiple_core_compare::same(base_map, vertex_index_map)) {
            LOG(logger, LOG_RANK::INFO) << "Branch Insertion2*," << maintenance_time << "\n";
        }

        G->remove_edge_collection(insertion_edge_vector);
    }

    /**
     * @brief restore the graph
     */
    G->insert_edge_collection(insertion_edge_vector);
}

void removal_test(const shared_ptr<temporal_graph>& G,
                  const shared_ptr<vector<shared_ptr<temporal_edge>>>& edge_vector,
                  uint32_t thread_number,
                  const string& log_path,
                  const string& input_file_name){
    const string log_file_name = log_path + "removal_test.log";
    scnu::simple_logger logger(log_file_name);


    LOG(logger, LOG_RANK::INFO) << input_file_name << "\n";

    const uint64_t m = 5000;
    auto removal_edge_vector = make_shared<vector<shared_ptr<temporal_edge>>>();
    for (uint32_t i = 0; i < min(m, edge_vector->size()); ++i) {
        removal_edge_vector->push_back(edge_vector->at(i));
    }

    auto previous_vertex_index_map = make_shared<unordered_map<uint32_t, shared_ptr<multiple_core_pair_map_index>>>();
    {
        auto pool = make_shared<thread_pool>(thread_number);

        auto vertex_edge_index_map = make_shared<unordered_map<uint32_t, shared_ptr<map<uint32_t, shared_ptr<unordered_set<uint32_t>>>>>>();
        branch_multiple_core_decomposition::init(G, vertex_edge_index_map, previous_vertex_index_map, pool);
        branch_multiple_core_decomposition::decompose(G, vertex_edge_index_map, previous_vertex_index_map, pool);
    }

    auto base_map = make_shared<unordered_map<uint32_t,shared_ptr<multiple_core_pair_map_index>>>();
    {
        G->remove_edge_collection(removal_edge_vector);

        simple_timer decomposition_timer;
        auto pool = make_shared<thread_pool>(thread_number);
        auto vertex_mutex_map = make_shared<unordered_map<uint32_t, shared_ptr<mutex>>>();
        branch_multiple_core_decomposition::init(G, vertex_mutex_map, base_map, pool);
        branch_multiple_core_decomposition::decompose(G, vertex_mutex_map, base_map, pool);
        auto decomposition_time = decomposition_timer.get_elapse_second();

        LOG(logger, LOG_RANK::INFO) << "Branch Decomposition2," << decomposition_time << "\n";

        G->insert_edge_collection(removal_edge_vector);
    }


    {
        auto base_map2 = make_shared<unordered_map<uint32_t,shared_ptr<multiple_core_pair_map_index>>>();
        G->remove_edge_collection(removal_edge_vector);

        simple_timer decomposition_timer;
        auto pool = make_shared<thread_pool>(thread_number);
        auto vertex_mutex_map = make_shared<unordered_map<uint32_t, shared_ptr<mutex>>>();
        auto vertex_edge_index_map = make_shared<unordered_map<uint32_t, shared_ptr<map<uint32_t, shared_ptr<unordered_set<uint32_t>>>>>>();
        branch_multiple_core_decomposition::init(G, vertex_mutex_map, vertex_edge_index_map, base_map2, pool);
        branch_multiple_core_decomposition::decompose(G, vertex_mutex_map, vertex_edge_index_map, base_map2, pool);
        auto decomposition_time = decomposition_timer.get_elapse_second();

        LOG(logger, LOG_RANK::INFO) << "Branch Decomposition2*," << decomposition_time << "\n";

        G->insert_edge_collection(removal_edge_vector);
    }


    {
        auto vertex_index_map = make_shared<unordered_map<uint32_t,shared_ptr<multiple_core_pair_map_index>>>();
        for(const auto&[u, u_index]:*previous_vertex_index_map){
            vertex_index_map->insert({u, make_shared<multiple_core_pair_map_index>(u_index)});
        }

        auto removal_edge_set = container_convert::to_unordered_set<shared_ptr<temporal_edge>>(removal_edge_vector);

        simple_timer maintenance_timer;
        auto pool = make_shared<thread_pool>(thread_number);
        auto vertex_mutex_map = make_shared<unordered_map<uint32_t, shared_ptr<mutex>>>();
        auto vertex_degree_map = make_shared<unordered_map<uint32_t, uint32_t>>();

        hierarchy_multiple_core_maintenance::init(G, vertex_mutex_map, vertex_degree_map, pool);
        hierarchy_multiple_core_maintenance::remove(G, vertex_mutex_map, removal_edge_set, vertex_index_map, vertex_degree_map, pool);

        auto maintenance_time = maintenance_timer.get_elapse_second();

        if(multiple_core_compare::same(base_map, vertex_index_map)){
            LOG(logger, LOG_RANK::INFO) << "Hierarchy Removal2," << maintenance_time << "\n";
        }

        G->insert_edge_collection(removal_edge_vector);
    }

    {
        auto vertex_index_map = make_shared<unordered_map<uint32_t, shared_ptr<multiple_core_pair_map_index>>>();
        for (const auto &[u, u_index]: *previous_vertex_index_map) {
            vertex_index_map->insert({u, make_shared<multiple_core_pair_map_index>(u_index)});
        }

        auto removal_edge_set = container_copy::to_unordered_set<shared_ptr<temporal_edge>>(removal_edge_vector);

        simple_timer maintenance_timer;
        auto pool = make_shared<thread_pool>(thread_number);
        auto vertex_mutex_map = make_shared<unordered_map<uint32_t, shared_ptr<mutex>>>();
        auto vertex_degree_map = make_shared<unordered_map<uint32_t, uint32_t>>();

        branch_multiple_core_maintenance::init(G, vertex_mutex_map, vertex_degree_map, pool);
        branch_multiple_core_maintenance::remove(G, vertex_mutex_map, removal_edge_set, vertex_index_map, vertex_degree_map, pool);

        auto maintenance_time = maintenance_timer.get_elapse_second();

        if (multiple_core_compare::same(base_map, vertex_index_map)) {
            LOG(logger, LOG_RANK::INFO) << "Branch Removal2," << maintenance_time << "\n";
        }

        G->insert_edge_collection(removal_edge_vector);
    }

    {
        auto vertex_index_map = make_shared<unordered_map<uint32_t, shared_ptr<multiple_core_pair_map_index>>>();
        for (const auto &[u, u_index]: *previous_vertex_index_map) {
            vertex_index_map->insert({u, make_shared<multiple_core_pair_map_index>(u_index)});
        }

        auto removal_edge_set = container_copy::to_unordered_set<shared_ptr<temporal_edge>>(removal_edge_vector);

        simple_timer maintenance_timer;
        auto pool = make_shared<thread_pool>(thread_number);
        auto vertex_mutex_map = make_shared<unordered_map<uint32_t, shared_ptr<mutex>>>();
        auto vertex_edge_map = make_shared<unordered_map<uint32_t, shared_ptr<map<uint32_t, shared_ptr<unordered_set<uint32_t>>>>>>();
        auto vertex_degree_map = make_shared<unordered_map<uint32_t, uint32_t>>();

        branch_multiple_core_maintenance::init(G, vertex_mutex_map, vertex_edge_map, vertex_degree_map, pool);
        branch_multiple_core_maintenance::remove(G,vertex_mutex_map, removal_edge_set, vertex_edge_map, vertex_index_map, vertex_degree_map, pool);

        auto maintenance_time = maintenance_timer.get_elapse_second();

        if (multiple_core_compare::same(base_map, vertex_index_map)) {
            LOG(logger, LOG_RANK::INFO) << "Branch Removal2*," << maintenance_time << "\n";
        }

        G->insert_edge_collection(removal_edge_vector);
    }
}

void insertion_removal_test(const shared_ptr<temporal_graph>& G,
                            const shared_ptr<vector<shared_ptr<temporal_edge>>>& edge_vector,
                            uint32_t thread_number,
                            const string& log_path,
                            const string& input_file_name){
    const string log_file_name = log_path + "insertion_removal_test.log";
    scnu::simple_logger logger(log_file_name);


    LOG(logger, LOG_RANK::INFO) << input_file_name << "\n";

    const uint32_t m = 2500;
    auto previous_edge_vector = make_shared<vector<shared_ptr<temporal_edge>>>();
    auto insertion_edge_vector = make_shared<vector<shared_ptr<temporal_edge>>>();
    auto removal_edge_vector = make_shared<vector<shared_ptr<temporal_edge>>>();
    for (uint32_t i = 0; i < edge_vector->size(); ++i) {
        if (i < m) {
            removal_edge_vector->push_back(edge_vector->at(i));
        }
        if (i < edge_vector->size() - m) {
            previous_edge_vector->push_back(edge_vector->at(i));
        } else {
            insertion_edge_vector->push_back(edge_vector->at(i));
        }
    }

    auto previous_vertex_index_map = make_shared<unordered_map<uint32_t, shared_ptr<multiple_core_pair_map_index>>>();
    {
        G->remove_edge_collection(insertion_edge_vector);
        auto pool = make_shared<thread_pool>(thread_number);

        auto vertex_edge_index_map = make_shared<unordered_map<uint32_t, shared_ptr<map<uint32_t, shared_ptr<unordered_set<uint32_t>>>>>>();
        branch_multiple_core_decomposition::init(G, vertex_edge_index_map, previous_vertex_index_map, pool);
        branch_multiple_core_decomposition::decompose(G, vertex_edge_index_map, previous_vertex_index_map, pool);
    }

    auto base_map = make_shared<unordered_map<uint32_t, shared_ptr<multiple_core_pair_map_index>>>();
    {
        G->insert_edge_collection(insertion_edge_vector);
        G->remove_edge_collection(removal_edge_vector);

        simple_timer decomposition_timer;
        auto pool = make_shared<thread_pool>(thread_number);
        auto vertex_mutex_map = make_shared<unordered_map<uint32_t, shared_ptr<mutex>>>();
        branch_multiple_core_decomposition::init(G, vertex_mutex_map, base_map, pool);
        branch_multiple_core_decomposition::decompose(G, vertex_mutex_map, base_map, pool);
        auto decomposition_time = decomposition_timer.get_elapse_second();

        LOG(logger, LOG_RANK::INFO) << "Branch Decomposition2," << decomposition_time << "\n";

        G->remove_edge_collection(insertion_edge_vector);
        G->insert_edge_collection(removal_edge_vector);
    }


    {
        auto base_map2 = make_shared<unordered_map<uint32_t, shared_ptr<multiple_core_pair_map_index>>>();
        G->insert_edge_collection(insertion_edge_vector);
        G->remove_edge_collection(removal_edge_vector);

        simple_timer decomposition_timer;
        auto pool = make_shared<thread_pool>(thread_number);
        auto vertex_mutex_map = make_shared<unordered_map<uint32_t, shared_ptr<mutex>>>();
        auto vertex_edge_index_map = make_shared<unordered_map<uint32_t, shared_ptr<map<uint32_t, shared_ptr<unordered_set<uint32_t>>>>>>();
        branch_multiple_core_decomposition::init(G, vertex_mutex_map, vertex_edge_index_map, base_map2, pool);
        branch_multiple_core_decomposition::decompose(G, vertex_mutex_map, vertex_edge_index_map, base_map2, pool);
        auto decomposition_time = decomposition_timer.get_elapse_second();

        LOG(logger, LOG_RANK::INFO) << "Branch Decomposition2*," << decomposition_time << "\n";

        G->remove_edge_collection(insertion_edge_vector);
        G->insert_edge_collection(removal_edge_vector);
    }


    {
        auto vertex_index_map = make_shared<unordered_map<uint32_t, shared_ptr<multiple_core_pair_map_index>>>();
        for (const auto &[u, u_index]: *previous_vertex_index_map) {
            vertex_index_map->insert({u, make_shared<multiple_core_pair_map_index>(u_index)});
        }

        auto insertion_edge_set = container_convert::to_unordered_set<shared_ptr<temporal_edge>>(insertion_edge_vector);
        auto removal_edge_set = container_convert::to_unordered_set<shared_ptr<temporal_edge>>(removal_edge_vector);

        simple_timer maintenance_timer;
        auto pool = make_shared<thread_pool>(thread_number);
        auto vertex_mutex_map = make_shared<unordered_map<uint32_t, shared_ptr<mutex>>>();
        auto vertex_degree_map = make_shared<unordered_map<uint32_t, uint32_t>>();

        hierarchy_multiple_core_maintenance::init(G, vertex_mutex_map, vertex_degree_map, pool);
        hierarchy_multiple_core_maintenance::insert(G, vertex_mutex_map, insertion_edge_set, vertex_index_map, vertex_degree_map, pool);
        hierarchy_multiple_core_maintenance::remove(G, vertex_mutex_map, removal_edge_set, vertex_index_map, vertex_degree_map, pool);


        auto maintenance_time = maintenance_timer.get_elapse_second();

        if (multiple_core_compare::same(base_map, vertex_index_map)) {
            LOG(logger, LOG_RANK::INFO) << "Hierarchy Maintenance2," << maintenance_time << "\n";
        }

        G->remove_edge_collection(insertion_edge_vector);
        G->insert_edge_collection(removal_edge_vector);
    }

    {
        auto vertex_index_map = make_shared<unordered_map<uint32_t, shared_ptr<multiple_core_pair_map_index>>>();
        for (const auto &[u, u_index]: *previous_vertex_index_map) {
            vertex_index_map->insert({u, make_shared<multiple_core_pair_map_index>(u_index)});
        }

        auto insertion_edge_set = container_convert::to_unordered_set<shared_ptr<temporal_edge>>(insertion_edge_vector);
        auto removal_edge_set = container_convert::to_unordered_set<shared_ptr<temporal_edge>>(removal_edge_vector);

        simple_timer maintenance_timer;
        auto pool = make_shared<thread_pool>(thread_number);

        auto vertex_mutex_map = make_shared<unordered_map<uint32_t, shared_ptr<mutex>>>();
        auto vertex_degree_map = make_shared<unordered_map<uint32_t, uint32_t>>();

        branch_multiple_core_maintenance::init(G, vertex_mutex_map, vertex_degree_map, pool);
        branch_multiple_core_maintenance::insert(G, vertex_mutex_map, insertion_edge_set, vertex_index_map, vertex_degree_map,pool);
        branch_multiple_core_maintenance::remove(G, vertex_mutex_map, removal_edge_set, vertex_index_map, vertex_degree_map, pool);


        auto maintenance_time = maintenance_timer.get_elapse_second();

        if (multiple_core_compare::same(base_map, vertex_index_map)) {
            LOG(logger, LOG_RANK::INFO) << "Branch Maintenance2," << maintenance_time << "\n";
        }

        G->remove_edge_collection(insertion_edge_vector);
        G->insert_edge_collection(removal_edge_vector);
    }

    {
        auto vertex_index_map = make_shared<unordered_map<uint32_t, shared_ptr<multiple_core_pair_map_index>>>();
        for (const auto &[u, u_index]: *previous_vertex_index_map) {
            vertex_index_map->insert({u, make_shared<multiple_core_pair_map_index>(u_index)});
        }

        auto insertion_edge_set = container_convert::to_unordered_set<shared_ptr<temporal_edge>>(insertion_edge_vector);
        auto removal_edge_set = container_convert::to_unordered_set<shared_ptr<temporal_edge>>(removal_edge_vector);

        simple_timer maintenance_timer;
        auto pool = make_shared<thread_pool>(thread_number);
        auto vertex_edge_index_map = make_shared<unordered_map<uint32_t, shared_ptr<map<uint32_t, shared_ptr<unordered_set<uint32_t>>>>>>();
        auto vertex_mutex_map = make_shared<unordered_map<uint32_t, shared_ptr<mutex>>>();
        auto vertex_degree_map = make_shared<unordered_map<uint32_t, uint32_t>>();

        branch_multiple_core_maintenance::init(G, vertex_mutex_map, vertex_edge_index_map, vertex_degree_map, pool);
        branch_multiple_core_maintenance::insert(G, vertex_mutex_map, insertion_edge_set, vertex_edge_index_map,
                                                 vertex_index_map, vertex_degree_map, pool);
        branch_multiple_core_maintenance::remove(G, vertex_mutex_map, removal_edge_set, vertex_edge_index_map,
                                                 vertex_index_map, vertex_degree_map, pool);

        auto maintenance_time = maintenance_timer.get_elapse_second();

        if (multiple_core_compare::same(base_map, vertex_index_map)) {
            LOG(logger, LOG_RANK::INFO) << "Branch Maintenance2*," << maintenance_time << "\n";
        }

        G->remove_edge_collection(insertion_edge_vector);
        G->insert_edge_collection(removal_edge_vector);
    }

    /**
     * @brief restore the graph
     */
    G->insert_edge_collection(insertion_edge_vector);
}


int main(int argc, char **argv){
    if(argc < 4){
        cout<<"Usage: path, file_name, and thread_number!\n";
        exit(0);
    }
    string const path = argv[1];
    string const input_file_name = argv[2];
    uint32_t const thread_number = std::stoul(argv[3]);

    /**
     *@brief prepare log file
     **/
    string const log_path = path + "log/";
    auto directory = std::filesystem::path(log_path);
    if (!exists(directory)) {
        create_directory(directory);
    }

    /**
    * @brief load global edge vector
    */
    auto edge_vector = temporal_graph_io::get_edge_vector(path, input_file_name);
    auto G = shared_ptr<temporal_graph>();
    {
        auto pool = make_shared<thread_pool>(thread_number);
        G = temporal_graph_io::load_graph(edge_vector, pool);
    }

    decomposition_test(G, thread_number, log_path, input_file_name);

    insertion_test(G, edge_vector, thread_number, log_path, input_file_name);

    removal_test(G, edge_vector, thread_number, log_path, input_file_name);

    insertion_removal_test(G, edge_vector, thread_number, log_path, input_file_name);

    return 0;
}

