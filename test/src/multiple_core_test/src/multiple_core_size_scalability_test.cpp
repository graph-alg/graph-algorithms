
#include "multiple_core_test/multiple_core_test.h"

void insertion_size_test(const shared_ptr<vector<shared_ptr<temporal_edge>>>& total_edge_vector,
                         uint32_t thread_number,
                         const string& log_path,
                         const string &input_file_name){
    const string log_file_name = log_path + "insertion_size_test.log";
    scnu::simple_logger logger(log_file_name);

    const uint32_t m = 25000;

    const vector<double> ratio_list{0.2, 0.4, 0.6, 0.8, 1.0};
    for(const auto &r:ratio_list){

        LOG(logger,LOG_RANK::INFO) << input_file_name  << "," << r << "\n";

        auto size = uint32_t (std::ceil(total_edge_vector->size() * r));
        auto edge_vector = make_shared<vector<shared_ptr<temporal_edge>>>();
        for(uint32_t i = 0; i < size;++i){
            edge_vector->push_back(total_edge_vector->at(i));
        }

        auto previous_edge_vector = make_shared<vector<shared_ptr<temporal_edge>>>();
        auto insertion_edge_vector = make_shared<vector<shared_ptr<temporal_edge>>>();
        for(uint32_t i = 0; i < edge_vector->size(); ++i){
            if(i < edge_vector->size() - m){
                previous_edge_vector->push_back(edge_vector->at(i));
            }else{
                insertion_edge_vector->push_back(edge_vector->at(i));
            }
        }

        auto previous_vertex_index_map = make_shared<unordered_map<uint32_t,shared_ptr<multiple_core_pair_map_index>>>();
        shared_ptr<temporal_graph> G;
        {
            auto pool = make_shared<thread_pool>(thread_number);
            G = temporal_graph_io::load_graph(previous_edge_vector, pool);

            auto vertex_edge_index_map = make_shared<unordered_map<uint32_t, shared_ptr<map<uint32_t, shared_ptr<unordered_set<uint32_t>>>>>>();
            branch_multiple_core_decomposition::init(G, vertex_edge_index_map, previous_vertex_index_map, pool);
            branch_multiple_core_decomposition::decompose(G, vertex_edge_index_map, previous_vertex_index_map, pool);
        }

        auto base_map = make_shared<unordered_map<uint32_t,shared_ptr<multiple_core_pair_map_index>>>();
        {

            G->insert_edge_collection(insertion_edge_vector);

            simple_timer decomposition_timer;
            auto pool = make_shared<thread_pool>(thread_number);
            auto vertex_mutex_map = make_shared<unordered_map<uint32_t, shared_ptr<mutex>>>();
            auto vertex_edge_index_map = make_shared<unordered_map<uint32_t, shared_ptr<map<uint32_t, shared_ptr<unordered_set<uint32_t>>>>>>();
            branch_multiple_core_decomposition::init(G, vertex_mutex_map, vertex_edge_index_map, base_map, pool);
            branch_multiple_core_decomposition::decompose(G, vertex_mutex_map, vertex_edge_index_map, base_map, pool);
            auto decomposition_time = decomposition_timer.get_elapse_second();

            LOG(logger, LOG_RANK::INFO) << "Branch Decomposition2*," << decomposition_time << "\n";

            G->remove_edge_collection(insertion_edge_vector);
        }

        {
            auto vertex_index_map = make_shared<unordered_map<uint32_t,shared_ptr<multiple_core_pair_map_index>>>();
            for(const auto&[u, u_index]:*previous_vertex_index_map){
                vertex_index_map->insert({u, make_shared<multiple_core_pair_map_index>(u_index)});
            }

            auto insertion_edge_set = container_convert::to_unordered_set<shared_ptr<temporal_edge>>(insertion_edge_vector);

            simple_timer maintenance_timer;
            auto pool = make_shared<thread_pool>(thread_number);
            auto vertex_mutex_map = make_shared<unordered_map<uint32_t, shared_ptr<mutex>>>();
            auto vertex_degree_map = make_shared<unordered_map<uint32_t, uint32_t>>();

            branch_multiple_core_maintenance::init(G, vertex_mutex_map, vertex_degree_map, pool);
            branch_multiple_core_maintenance::insert(G, vertex_mutex_map, insertion_edge_set, vertex_index_map, vertex_degree_map, pool);

            auto maintenance_time = maintenance_timer.get_elapse_second();

            if(multiple_core_compare::same(base_map, vertex_index_map)){
                LOG(logger, LOG_RANK::INFO) << "Branch Insertion2," << maintenance_time << "\n";
            }

            G->remove_edge_collection(insertion_edge_vector);
        }

        {
            auto vertex_index_map = make_shared<unordered_map<uint32_t,shared_ptr<multiple_core_pair_map_index>>>();
            for(const auto&[u, u_index]:*previous_vertex_index_map){
                vertex_index_map->insert({u, make_shared<multiple_core_pair_map_index>(u_index)});
            }

            auto insertion_edge_set = container_convert::to_unordered_set<shared_ptr<temporal_edge>>(insertion_edge_vector);

            simple_timer maintenance_timer;
            auto pool = make_shared<thread_pool>(thread_number);
            auto vertex_edge_index_map = make_shared<unordered_map<uint32_t, shared_ptr<map<uint32_t, shared_ptr<unordered_set<uint32_t>>>>>>();
            auto vertex_mutex_map = make_shared<unordered_map<uint32_t, shared_ptr<mutex>>>();
            auto vertex_degree_map = make_shared<unordered_map<uint32_t, uint32_t>>();

            branch_multiple_core_maintenance::init(G, vertex_mutex_map, vertex_edge_index_map, vertex_degree_map, pool);
            branch_multiple_core_maintenance::insert(G, vertex_mutex_map, insertion_edge_set, vertex_edge_index_map, vertex_index_map, vertex_degree_map, pool);

            auto maintenance_time = maintenance_timer.get_elapse_second();

            if(multiple_core_compare::same(base_map, vertex_index_map)){
                LOG(logger, LOG_RANK::INFO) << "Branch Insertion2*," << maintenance_time << "\n";
            }

            G->remove_edge_collection(insertion_edge_vector);
        }
    }
}

void removal_size_test(const shared_ptr<vector<shared_ptr<temporal_edge>>>& total_edge_vector,
                       uint32_t thread_number,
                       const string& log_path,
                       const string &input_file_name){
    const string log_file_name = log_path + "removal_size_test.log";
    scnu::simple_logger logger(log_file_name);

    const uint64_t m = 25000;
    const vector<double> ratio_list{0.2, 0.4, 0.6, 0.8, 1.0};
    for(const auto &r:ratio_list) {

        LOG(logger,LOG_RANK::INFO) << input_file_name  << "," << r << "\n";

        auto size = uint32_t(std::ceil(total_edge_vector->size() * r));
        auto edge_vector = make_shared<vector<shared_ptr<temporal_edge>>>();
        for (uint32_t i = 0; i < size; ++i) {
            edge_vector->emplace_back(total_edge_vector->at(i));
        }

        auto removal_edge_vector = make_shared<vector<shared_ptr<temporal_edge>>>();
        for (uint32_t i = 0; i < min(m, edge_vector->size()); ++i) {
            removal_edge_vector->emplace_back(edge_vector->at(i));
        }

        auto previous_vertex_index_map = make_shared<unordered_map<uint32_t, shared_ptr<multiple_core_pair_map_index>>>();
        auto G = shared_ptr<temporal_graph>();
        {
            auto pool = make_shared<thread_pool>(thread_number);
            G = temporal_graph_io::load_graph(edge_vector, pool);

            auto vertex_edge_index_map = make_shared<unordered_map<uint32_t, shared_ptr<map<uint32_t, shared_ptr<unordered_set<uint32_t>>>>>>();
            branch_multiple_core_decomposition::init(G, vertex_edge_index_map, previous_vertex_index_map, pool);
            branch_multiple_core_decomposition::decompose(G, vertex_edge_index_map, previous_vertex_index_map, pool);
        }

        auto base_map = make_shared<unordered_map<uint32_t, shared_ptr<multiple_core_pair_map_index>>>();
        {

            G->remove_edge_collection(removal_edge_vector);

            simple_timer decomposition_timer;
            auto pool = make_shared<thread_pool>(thread_number);
            auto vertex_mutex_map = make_shared<unordered_map<uint32_t, shared_ptr<mutex>>>();
            auto vertex_edge_index_map = make_shared<unordered_map<uint32_t, shared_ptr<map<uint32_t, shared_ptr<unordered_set<uint32_t>>>>>>();
            branch_multiple_core_decomposition::init(G, vertex_mutex_map, vertex_edge_index_map, base_map, pool);
            branch_multiple_core_decomposition::decompose(G, vertex_mutex_map, vertex_edge_index_map, base_map, pool);
            auto decomposition_time = decomposition_timer.get_elapse_second();

            LOG(logger, LOG_RANK::INFO) << "Branch Decomposition2*," << decomposition_time << "\n";

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

            branch_multiple_core_maintenance::init(G, vertex_mutex_map,
                                                   vertex_degree_map, pool);
            branch_multiple_core_maintenance::remove(G, vertex_mutex_map, removal_edge_set,
                                                     vertex_index_map, vertex_degree_map, pool);

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
            auto vertex_edge_index_map = make_shared<unordered_map<uint32_t, shared_ptr<map<uint32_t, shared_ptr<unordered_set<uint32_t>>>>>>();
            auto vertex_degree_map = make_shared<unordered_map<uint32_t, uint32_t>>();

            branch_multiple_core_maintenance::init(G, vertex_mutex_map, vertex_edge_index_map,
                                                   vertex_degree_map, pool);
            branch_multiple_core_maintenance::remove(G, vertex_mutex_map, removal_edge_set, vertex_edge_index_map,
                                                     vertex_index_map, vertex_degree_map, pool);

            auto maintenance_time = maintenance_timer.get_elapse_second();

            if (multiple_core_compare::same(base_map, vertex_index_map)) {
                LOG(logger, LOG_RANK::INFO) << "Branch Removal2*," << maintenance_time << "\n";
            }

            G->insert_edge_collection(removal_edge_vector);
        }
    }
}

void insertion_removal_size_test(const shared_ptr<vector<shared_ptr<temporal_edge>>>& total_edge_vector,
                                 uint32_t thread_number,
                                 const string& log_path,
                                 const string &input_file_name){
    const string log_file_name = log_path + "insertion_removal_size_test.log";
    scnu::simple_logger logger(log_file_name);

    const uint32_t m = 12500;

    const vector<double> ratio_list{0.2, 0.4, 0.6, 0.8, 1.0};
    for(const auto &r:ratio_list) {

        LOG(logger,LOG_RANK::INFO) << input_file_name  << "," << r << "\n";

        auto size = uint32_t(std::ceil(total_edge_vector->size() * r));
        auto edge_vector = make_shared<vector<shared_ptr<temporal_edge>>>();
        for (uint32_t i = 0; i < size; ++i) {
            edge_vector->emplace_back(total_edge_vector->at(i));
        }

        auto previous_edge_vector = make_shared<vector<shared_ptr<temporal_edge>>>();
        auto insertion_edge_vector = make_shared<vector<shared_ptr<temporal_edge>>>();
        auto removal_edge_vector = make_shared<vector<shared_ptr<temporal_edge>>>();
        for (uint32_t i = 0; i < edge_vector->size(); ++i) {
            if (i < m) {
                removal_edge_vector->emplace_back(edge_vector->at(i));
            }
            if (i < edge_vector->size() - m) {
                previous_edge_vector->emplace_back(edge_vector->at(i));
            } else {
                insertion_edge_vector->emplace_back(edge_vector->at(i));
            }
        }

        auto previous_vertex_index_map = make_shared<unordered_map<uint32_t, shared_ptr<multiple_core_pair_map_index>>>();
        shared_ptr<temporal_graph> G;
        {
            auto pool = make_shared<thread_pool>(thread_number);
            G = temporal_graph_io::load_graph(previous_edge_vector, pool);

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
            auto vertex_edge_index_map = make_shared<unordered_map<uint32_t, shared_ptr<map<uint32_t, shared_ptr<unordered_set<uint32_t>>>>>>();
            branch_multiple_core_decomposition::init(G, vertex_mutex_map, vertex_edge_index_map, base_map, pool);
            branch_multiple_core_decomposition::decompose(G, vertex_mutex_map, vertex_edge_index_map, base_map, pool);
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

            auto insertion_edge_set = container_convert::to_unordered_set<shared_ptr<temporal_edge>>(
                    insertion_edge_vector);
            auto removal_edge_set = container_convert::to_unordered_set<shared_ptr<temporal_edge>>(removal_edge_vector);

            simple_timer maintenance_timer;
            auto pool = make_shared<thread_pool>(thread_number);
            auto vertex_mutex_map = make_shared<unordered_map<uint32_t, shared_ptr<mutex>>>();
            auto vertex_degree_map = make_shared<unordered_map<uint32_t, uint32_t>>();

            branch_multiple_core_maintenance::init(G, vertex_mutex_map,
                                                   vertex_degree_map, pool);
            branch_multiple_core_maintenance::insert(G, vertex_mutex_map, insertion_edge_set,
                                                     vertex_index_map, vertex_degree_map, pool);
            branch_multiple_core_maintenance::remove(G, vertex_mutex_map, removal_edge_set,
                                                     vertex_index_map, vertex_degree_map, pool);

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

            auto insertion_edge_set = container_convert::to_unordered_set<shared_ptr<temporal_edge>>(
                    insertion_edge_vector);
            auto removal_edge_set = container_convert::to_unordered_set<shared_ptr<temporal_edge>>(removal_edge_vector);

            simple_timer maintenance_timer;
            auto pool = make_shared<thread_pool>(thread_number);
            auto vertex_edge_index_map = make_shared<unordered_map<uint32_t, shared_ptr<map<uint32_t, shared_ptr<unordered_set<uint32_t>>>>>>();
            auto vertex_mutex_map = make_shared<unordered_map<uint32_t, shared_ptr<mutex>>>();
            auto vertex_degree_map = make_shared<unordered_map<uint32_t, uint32_t>>();

            branch_multiple_core_maintenance::init(G, vertex_mutex_map, vertex_edge_index_map,
                                                   vertex_degree_map, pool);
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
    }
}


int main(int argc, char **argv){
    if(argc < 4){
        cout<<"Usage: path, file_name, and thread_number!\n";
        exit(0);
    }

    const string path = argv[1];
    const string input_file_name = argv[2];
    const uint32_t thread_number = std::stoul(argv[3]);

    /**
     *@brief prepare log file
     **/
    const string log_path = path + "log/";
    auto directory = std::filesystem::path(log_path);
    if (!exists(directory)) {
        create_directory(directory);
    }

    /**
      * @brief load global edge vector
      */

    auto edge_vector = temporal_graph_io::get_edge_vector(path, input_file_name);
    {
//        auto rd = make_shared<random_device>();
//        std::shuffle(edge_vector->begin(), edge_vector->end(), *random_generator::get_default_engine(rd));
//        ofstream f(log_path + "insertion_compare_error", std::ios::out);
//        for (const auto &e: *edge_vector) {
//            f << e->get_source_vertex_id() << "," << e->get_destination_vertex_id() << "," << std::setprecision(0) <<
//              std::to_string(e->get_weight()) << "," << e->get_timestamp() << '\n';
//        }
    }

    {
        insertion_size_test(edge_vector, thread_number,log_path, input_file_name);
    }

    {
        removal_size_test(edge_vector, thread_number,log_path, input_file_name);
    }

    {
        insertion_removal_size_test(edge_vector, thread_number,log_path, input_file_name);
    }

    return 0;
}

