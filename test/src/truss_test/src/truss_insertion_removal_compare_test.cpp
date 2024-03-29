
#include "truss_test/truss_test.h"

shared_ptr<abstract_graph> load_graph(const string &path,
                                      const string &input_file_name,
                                      uint32_t thread_number,
                                      const shared_ptr<vector<shared_ptr<abstract_edge>>> &total_insertion_edge_vector,
                                      const shared_ptr<vector<shared_ptr<abstract_edge>>> &total_removal_edge_vector){
    auto edge_vector = abstract_graph_io::get_edge_vector(path, input_file_name);

    for (uint32_t i = 0; i < 100000; ++i) {
            total_removal_edge_vector->push_back(edge_vector->at(i));
    }

    for (uint32_t i = edge_vector->size() - 100000; i < edge_vector->size(); ++i) {
        total_insertion_edge_vector->push_back(edge_vector->at(i));
    }

    auto pool = make_shared<thread_pool>(thread_number);
    auto G = abstract_graph_io::load_graph(edge_vector, edge_vector->size() - 100000, pool);

    return  G;
}

void prepare(const shared_ptr<abstract_graph>&G,
             const shared_ptr<unordered_map<shared_ptr<abstract_edge>, uint32_t>> &previous_edge_truss_map,
             const shared_ptr<unordered_map<shared_ptr<abstract_edge>, uint32_t>> &previous_edge_truss_support_map,
             const shared_ptr<unordered_map<uint32_t, shared_ptr<extend_list<int, shared_ptr<abstract_edge>>>>> &previous_truss_order_map,
             const shared_ptr<unordered_map<shared_ptr<abstract_edge>, uint32_t>> &previous_rem,
             uint32_t thread_number){
    auto edge_mutex_map = make_shared<unordered_map<shared_ptr<abstract_edge>, shared_ptr<mutex>>>();
    auto edge_rank_map = make_shared<unordered_map<shared_ptr<abstract_edge>, uint32_t>>();
    auto edge_support_map = make_shared<unordered_map<shared_ptr<abstract_edge>, uint32_t>>();

    auto pool = make_shared<thread_pool>(thread_number);

    basic_truss_decomposition::init(G, edge_mutex_map, edge_rank_map, edge_support_map, previous_edge_truss_map,
                                    previous_edge_truss_support_map, previous_rem, pool);
    basic_truss_decomposition::decompose(G, edge_mutex_map, edge_rank_map, edge_support_map,
                                         previous_edge_truss_map, previous_edge_truss_support_map,
                                         previous_truss_order_map, previous_rem, pool);
}

double decompose(const shared_ptr<abstract_graph> &G,
                 const shared_ptr<vector<shared_ptr<abstract_edge>>> &insertion_edge_vector,
                 const shared_ptr<vector<shared_ptr<abstract_edge>>> &removal_edge_vector,
                 const shared_ptr<unordered_map<shared_ptr<abstract_edge>, uint32_t>> &contrast_edge_truss_map,
                 uint32_t thread_number) {
    G->insert_edge_collection(insertion_edge_vector);
    G->remove_edge_collection(removal_edge_vector);

    simple_timer decomposition_timer;
    auto edge_mutex_map = make_shared<unordered_map<shared_ptr<abstract_edge>, shared_ptr<mutex>>>();
    auto edge_rank_map = make_shared<unordered_map<shared_ptr<abstract_edge>, uint32_t>>();
    auto edge_support_map = make_shared<unordered_map<shared_ptr<abstract_edge>, uint32_t>>();

    auto pool = make_shared<thread_pool>(thread_number);
    basic_truss_decomposition::init(G, edge_mutex_map, edge_rank_map, edge_support_map, contrast_edge_truss_map,
                                    pool);
    basic_truss_decomposition::decompose(G, edge_mutex_map, edge_rank_map, edge_support_map,
                                         contrast_edge_truss_map, pool);
    auto decomposition_time = decomposition_timer.get_elapse_second();


    G->remove_edge_collection(insertion_edge_vector);
    G->insert_edge_collection(removal_edge_vector);

    return decomposition_time;
}

double order_maintenance(const shared_ptr<abstract_graph> &G,
                         const shared_ptr<vector<shared_ptr<abstract_edge>>> &insertion_edge_vector,
                         const shared_ptr<vector<shared_ptr<abstract_edge>>> &removal_edge_vector,
                         const shared_ptr<unordered_map<shared_ptr<abstract_edge>, uint32_t>> &order_edge_truss_map,
                         const shared_ptr<unordered_map<shared_ptr<abstract_edge>, uint32_t>> &previous_edge_truss_support_map,
                         const shared_ptr<unordered_map<uint32_t, shared_ptr<extend_list<int, shared_ptr<abstract_edge>>>>> &previous_truss_order_map,
                         const shared_ptr<unordered_map<shared_ptr<abstract_edge>, uint32_t>> &previous_rem){

    auto inserted_edge_set = container_copy::to_unordered_set<shared_ptr<abstract_edge>>(insertion_edge_vector);
    auto removal_edge_set = container_copy::to_unordered_set<shared_ptr<abstract_edge>>(removal_edge_vector);
    auto ext = make_shared<unordered_map<shared_ptr<abstract_edge>, uint32_t>>();
    auto rem = container_copy::to_unordered_map<shared_ptr<abstract_edge>, uint32_t>(previous_rem);
    auto truss_order_map = make_shared<unordered_map<uint32_t, shared_ptr<extend_list<int, shared_ptr<abstract_edge>>>>>();
    for (const auto &[k, k_list]: *previous_truss_order_map) {
        auto copy_list = make_shared<extend_list<int, shared_ptr<abstract_edge>>>();
        for (auto p = k_list->get_head(); p; p = p->get_next()) {
            copy_list->push_back(p->get_value());
        }
        truss_order_map->insert({k, copy_list});
    }
    auto ts = container_copy::to_unordered_map<shared_ptr<abstract_edge>, uint32_t>(
            previous_edge_truss_support_map);

    simple_timer maintenance_timer;
    order_truss_maintenance::init(order_edge_truss_map, ext);
    order_truss_maintenance::remove(G, removal_edge_set, order_edge_truss_map, ts, truss_order_map, rem);
    order_truss_maintenance::insert(G, inserted_edge_set, order_edge_truss_map, ts, truss_order_map, rem, ext);

    auto maintenance_time = maintenance_timer.get_elapse_second();

    G->remove_edge_collection(insertion_edge_vector);
    G->insert_edge_collection(removal_edge_vector);

    return maintenance_time;
}

double jes_order_maintenance(const shared_ptr<abstract_graph> &G,
                             const shared_ptr<vector<shared_ptr<abstract_edge>>> &insertion_edge_vector,
                             const shared_ptr<vector<shared_ptr<abstract_edge>>> &removal_edge_vector,
                             const shared_ptr<unordered_map<shared_ptr<abstract_edge>, uint32_t>> &parallel_edge_truss_map,
                             const shared_ptr<unordered_map<shared_ptr<abstract_edge>, uint32_t>> &previous_edge_truss_support_map,
                             const shared_ptr<unordered_map<uint32_t, shared_ptr<extend_list<int, shared_ptr<abstract_edge>>>>> &previous_truss_order_map,
                             const shared_ptr<unordered_map<shared_ptr<abstract_edge>, uint32_t>> &previous_rem,
                             uint32_t thread_number){
    auto inserted_edge_set = container_copy::to_unordered_set<shared_ptr<abstract_edge>>(insertion_edge_vector);
    auto removed_edge_set = container_copy::to_unordered_set<shared_ptr<abstract_edge>>(removal_edge_vector);
    auto truss_order_map = make_shared<unordered_map<uint32_t, shared_ptr<extend_list<int, shared_ptr<abstract_edge>>>>>();
    for (const auto &[k, k_list]: *previous_truss_order_map) {
        auto copy_list = make_shared<extend_list<int, shared_ptr<abstract_edge>>>();
        for (auto p = k_list->get_head(); p; p = p->get_next()) {
            copy_list->push_back(p->get_value());
        }
        truss_order_map->insert({k, copy_list});
    }
    auto TS = container_copy::to_unordered_map<shared_ptr<abstract_edge>, uint32_t>(
            previous_edge_truss_support_map);
    auto rem = container_copy::to_unordered_map<shared_ptr<abstract_edge>, uint32_t>(
            previous_rem);

    simple_timer maintenance_timer;
    auto pool = make_shared<thread_pool>(thread_number);
    jes_order_truss_maintenance::insert(G, inserted_edge_set, parallel_edge_truss_map, truss_order_map, rem, TS,
                                        pool);
    jes_order_truss_maintenance::remove(G, removed_edge_set, parallel_edge_truss_map, truss_order_map, rem, TS,
                                        pool);
    auto maintenance_time = maintenance_timer.get_elapse_second();

    G->remove_edge_collection(insertion_edge_vector);
    G->insert_edge_collection(removal_edge_vector);

    return maintenance_time;
}

int main(int argc, char **argv) {
    if (argc < 4) {
        std::cout << "Usage: input path, file name, and thread number!";
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

    string const log_file_name = log_path + "insertion_removal_compare_test.log";
    simple_logger logger(log_file_name);

    {
//        auto rd = make_shared<random_device>();
//        shuffle(edge_vector->begin(), edge_vector->end(), *random_generator::get_default_engine(rd));
//
//        for (const auto &e: *edge_vector) {
//            std::cout << e->get_source_vertex_id() << "," << e->get_destination_vertex_id() << '\n';
//        }
    }

    auto total_insertion_edge_vector = make_shared<vector<shared_ptr<abstract_edge>>>();
    auto total_removal_edge_vector = make_shared<vector<shared_ptr<abstract_edge>>>();

    auto G = load_graph(path, input_file_name, thread_number, total_insertion_edge_vector, total_removal_edge_vector);

    auto previous_edge_truss_map = make_shared<unordered_map<shared_ptr<abstract_edge>, uint32_t>>();
    auto previous_edge_truss_support_map = make_shared<unordered_map<shared_ptr<abstract_edge>, uint32_t>>();
    auto previous_truss_order_map = make_shared<unordered_map<uint32_t, shared_ptr<extend_list<int, shared_ptr<abstract_edge>>>>>();
    auto previous_rem = make_shared<unordered_map<shared_ptr<abstract_edge>, uint32_t>>();
    {
        prepare(G, previous_edge_truss_map, previous_edge_truss_support_map, previous_truss_order_map, previous_rem, thread_number);
    }

    vector<uint32_t> m_array{20000, 40000, 60000, 80000, 100000};
    for(const auto &m:m_array){
        LOG(logger, LOG_RANK::INFO) << input_file_name << "," << m << "\n";

        auto insertion_edge_vector = make_shared<vector<shared_ptr<abstract_edge>>>(m);
        auto removal_edge_vector = make_shared<vector<shared_ptr<abstract_edge>>>(m);
        for (uint32_t i = 0; i < m; ++i) {
            insertion_edge_vector->at(i) = total_insertion_edge_vector->at(i);
            removal_edge_vector->at(i) = total_removal_edge_vector->at(i);
        }

//        auto contrast_edge_truss_map = make_shared<unordered_map<shared_ptr<abstract_edge>, uint32_t>>();
//        {
//            auto decomposition_time = decompose(G, insertion_edge_vector, removal_edge_vector, contrast_edge_truss_map,
//                                                thread_number);
//            LOG(logger, LOG_RANK::INFO) << "Decomposition," << decomposition_time << "\n";
//        }

        auto order_edge_truss_map = container_copy::to_unordered_map<shared_ptr<abstract_edge>, uint32_t>(
                previous_edge_truss_map);
        {
            auto maintenance_time = order_maintenance(G, insertion_edge_vector, removal_edge_vector,
                                                      order_edge_truss_map,
                                                      previous_edge_truss_support_map, previous_truss_order_map,
                                                      previous_rem);

            //if (truss_compare::same_associative_map(order_edge_truss_map, contrast_edge_truss_map)) {
            LOG(logger, LOG_RANK::INFO) << "Order Maintenance," << maintenance_time << "\n";
            //}
        }

        {
            auto parallel_edge_truss_map = container_copy::to_unordered_map<shared_ptr<abstract_edge>, uint32_t>(
                    previous_edge_truss_map);

            auto maintenance_time = jes_order_maintenance(G, insertion_edge_vector, removal_edge_vector, parallel_edge_truss_map,
                                                                previous_edge_truss_support_map, previous_truss_order_map, previous_rem, thread_number);

            if (truss_compare::same_associative_map(parallel_edge_truss_map, order_edge_truss_map)) {
                LOG(logger, LOG_RANK::INFO) << "Parallel Maintenance," << maintenance_time << "\n";
            }
        }
    }
}

