
#include "wing_test/wing_test.h"
#include <iostream>
shared_ptr<abstract_bipartite_graph> load_graph(const string &path,
                                                const string &input_file_name,
                                                uint32_t m,
                                                uint32_t thread_number,
                                                const shared_ptr<vector<shared_ptr<abstract_bipartite_edge>>> &total_insertion_edge_vector){
    auto edge_vector = abstract_bipartite_graph_io::get_edge_vector(path, input_file_name);
    auto previous_edge_vector = make_shared<vector<shared_ptr<abstract_bipartite_edge>>>();
    for (uint32_t i = 0; i < edge_vector->size(); i++) {
        if (i < edge_vector->size() - m) {
            previous_edge_vector->push_back(edge_vector->at(i));
        } else {
            total_insertion_edge_vector->push_back(edge_vector->at(i));
        }
    }

    auto pool = make_shared<thread_pool>(thread_number);
    auto G = abstract_bipartite_graph_io::load_graph(previous_edge_vector, pool);
    return  G;
}

void prepare(const shared_ptr<abstract_bipartite_graph>&G,
             const shared_ptr<unordered_map<shared_ptr<abstract_bipartite_edge>, uint32_t>> &previous_edge_wing_map,
             const shared_ptr<unordered_map<shared_ptr<abstract_bipartite_edge>, uint32_t>> &previous_WS,
             uint32_t &previous_k_max,
             uint32_t thread_number){
    auto edge_mutex_map = make_shared<unordered_map<shared_ptr<abstract_bipartite_edge>, shared_ptr<mutex>>>();
    auto edge_rank_map = make_shared<unordered_map<shared_ptr<abstract_bipartite_edge>, uint32_t>>();
    auto edge_support_map = make_shared<unordered_map<shared_ptr<abstract_bipartite_edge>, uint32_t>>();
    auto pool = make_shared<thread_pool>(thread_number);
    index_wing_decomposition::init(G, edge_mutex_map, edge_rank_map, edge_support_map,
                                   previous_edge_wing_map, pool);
    index_wing_decomposition::decompose(G, edge_mutex_map, edge_rank_map, edge_support_map,
                                        previous_edge_wing_map, pool);

    previous_k_max = index_wing_decomposition::WS_init(G, edge_mutex_map, edge_rank_map, previous_edge_wing_map, previous_WS, pool);
    std::cout<<"pre k max: "<<previous_k_max<<std::endl;
}

double decompose(const shared_ptr<abstract_bipartite_graph> &G,
                 const shared_ptr<vector<shared_ptr<abstract_bipartite_edge>>> &insertion_edge_vector,
                 const shared_ptr<unordered_map<shared_ptr<abstract_bipartite_edge>, uint32_t>> &contrastive_edge_wing_map,
                 uint32_t thread_number) {
    G->insert_edge_collection(insertion_edge_vector);



    auto edge_mutex_map = make_shared<unordered_map<shared_ptr<abstract_bipartite_edge>, shared_ptr<mutex>>>();
    auto edge_rank_map = make_shared< unordered_map<shared_ptr<abstract_bipartite_edge>, uint32_t>>();
    auto edge_support_map = make_shared<unordered_map<shared_ptr<abstract_bipartite_edge>, uint32_t>>();
    auto pool = make_shared<thread_pool>(thread_number);
    simple_timer bottom_up_timer;
    index_wing_decomposition::init(G, edge_mutex_map, edge_rank_map, edge_support_map,
                                   contrastive_edge_wing_map, pool);
    index_wing_decomposition::decompose(G, edge_mutex_map, edge_rank_map, edge_support_map,
                                        contrastive_edge_wing_map, pool);
    auto bottom_up_time = bottom_up_timer.get_elapse_second();

    G->remove_edge_collection(insertion_edge_vector);

    return bottom_up_time;
}

double BS_maintenance(const shared_ptr<abstract_bipartite_graph> &G,
                      const shared_ptr<vector<shared_ptr<abstract_bipartite_edge>>> &insertion_edge_vector,
                      const shared_ptr<unordered_map<shared_ptr<abstract_bipartite_edge>, uint32_t>> &quasi_edge_wing_map,
                      const shared_ptr<unordered_map<shared_ptr<abstract_bipartite_edge>, uint32_t>> &previous_WS,
                      uint32_t previous_k_max,
                      uint32_t thread_number){
    auto WS = make_shared<unordered_map<shared_ptr<abstract_bipartite_edge>, uint32_t>> ();
    auto insertion_edge_set = container_copy::to_unordered_set<shared_ptr<abstract_bipartite_edge>>(
            insertion_edge_vector);
    auto k_max = make_shared<uint32_t>(previous_k_max);
    auto edge_mutex_map = make_shared<unordered_map<shared_ptr<abstract_bipartite_edge>, shared_ptr<mutex>>>();
    auto edge_rank_map = make_shared< unordered_map<shared_ptr<abstract_bipartite_edge>, uint32_t>>();
    auto edge_support_map = make_shared<unordered_map<shared_ptr<abstract_bipartite_edge>, uint32_t>>();
    auto pool = make_shared<thread_pool>(thread_number);

    simple_timer ws_timer;
    quasi_wing_maintenance::init(quasi_edge_wing_map, edge_mutex_map, edge_rank_map, edge_support_map, WS, pool);
    quasi_wing_maintenance::insert(G, edge_mutex_map, insertion_edge_set, edge_rank_map, edge_support_map, quasi_edge_wing_map,
                                   WS,
                                   k_max,
                                   pool);
    auto bs_time = ws_timer.get_elapse_second();

    G->remove_edge_collection(insertion_edge_vector);

    return bs_time;
}

double WS_maintenance(const shared_ptr<abstract_bipartite_graph> &G,
                      const shared_ptr<vector<shared_ptr<abstract_bipartite_edge>>> &insertion_edge_vector,
                      const shared_ptr<unordered_map<shared_ptr<abstract_bipartite_edge>, uint32_t>> &ws_edge_wing_map,
                      const shared_ptr<unordered_map<shared_ptr<abstract_bipartite_edge>, uint32_t>> &previous_WS,
                      uint32_t previous_k_max,
                      uint32_t thread_number){
    auto WS = container_copy::to_unordered_map<shared_ptr<abstract_bipartite_edge>, uint32_t>(previous_WS);
    auto insertion_edge_set = container_copy::to_unordered_set<shared_ptr<abstract_bipartite_edge>>(
            insertion_edge_vector);
    auto k_max = make_shared<uint32_t>(previous_k_max);
    auto edge_mutex_map = make_shared<unordered_map<shared_ptr<abstract_bipartite_edge>, shared_ptr<mutex>>>();
    auto edge_rank_map = make_shared< unordered_map<shared_ptr<abstract_bipartite_edge>, uint32_t>>();
    auto edge_support_map = make_shared<unordered_map<shared_ptr<abstract_bipartite_edge>, uint32_t>>();
    auto pool = make_shared<thread_pool>(thread_number);

    simple_timer ws_timer;
    quasi_wing_maintenance::init(ws_edge_wing_map, edge_mutex_map, edge_rank_map, edge_support_map, pool);
    quasi_wing_maintenance::insert2(G, edge_mutex_map, insertion_edge_set, edge_rank_map, edge_support_map, ws_edge_wing_map,
                                    WS,
                                    k_max,
                                    pool);
    auto ws_time = ws_timer.get_elapse_second();

    G->remove_edge_collection(insertion_edge_vector);

    return ws_time;
}

void wing_insertion_thread_func(const shared_ptr<abstract_bipartite_graph> &G,
                                const shared_ptr<vector<shared_ptr<abstract_bipartite_edge>>> &insertion_edge_vector,
                                const shared_ptr<unordered_map<shared_ptr<abstract_bipartite_edge>, uint32_t>> &previous_edge_wing_map,
                                const shared_ptr<unordered_map<shared_ptr<abstract_bipartite_edge>, uint32_t>> &contrastive_edge_wing_map,
                                const shared_ptr<unordered_map<shared_ptr<abstract_bipartite_edge>, uint32_t>> &previous_WS,
                                uint32_t previous_k_max,
                                simple_logger &logger,
                                uint32_t thread_number){
    LOG(logger, LOG_RANK::INFO)  << thread_number <<'\n';

    {
        auto bu_edge_wing_map = make_shared<unordered_map<shared_ptr<abstract_bipartite_edge>, uint32_t>>();
        auto bu_time = decompose(G,insertion_edge_vector,bu_edge_wing_map,thread_number);

        if (wing_compare::same_associative_map(bu_edge_wing_map, contrastive_edge_wing_map)) {
            LOG(logger, LOG_RANK::INFO) << "BU Decomposition," << bu_time << "\n";
        }

    }

    {
        auto quasi_edge_wing_map = container_copy::to_unordered_map<shared_ptr<abstract_bipartite_edge>, uint32_t>(
                previous_edge_wing_map);
        auto bs_time = BS_maintenance(G,insertion_edge_vector,quasi_edge_wing_map,previous_WS,previous_k_max,thread_number);

        if (wing_compare::same_associative_map(quasi_edge_wing_map, contrastive_edge_wing_map)) {
            LOG(logger, LOG_RANK::INFO) << "Basic Insertion," << bs_time << "\n";
        }

    }

    {
        auto ws_edge_wing_map = container_copy::to_unordered_map<shared_ptr<abstract_bipartite_edge>, uint32_t>(
                previous_edge_wing_map);
        auto ws_time = WS_maintenance(G,insertion_edge_vector,ws_edge_wing_map,previous_WS,previous_k_max,thread_number);

        if (wing_compare::same_associative_map(ws_edge_wing_map, contrastive_edge_wing_map)) {
            LOG(logger, LOG_RANK::INFO) << "WS Insertion," << ws_time << "\n";
        }

    }
}

int main(int argc, char **argv) {
    string path = argv[1];
    string input_file_name = argv[2];
    uint32_t total_thread_number = std::stoul(argv[3]);

    if(argc < 2){
        std::cout<<"Usage: require  path and filename!"<<std::endl;
    }

    string log_path = path + "log/";
    string log_file_name = log_path + input_file_name + "_insertion_thread_test.log";
    auto directory = std::filesystem::path(log_path.c_str());
    if (!exists(directory)) {
        create_directory(directory);
    }

    /**
     * @brief prepare log file
     */
    simple_logger logger(log_file_name);

    uint32_t m = 20000;

    auto insertion_edge_vector = make_shared<vector<shared_ptr<abstract_bipartite_edge>>>();
    auto G = load_graph(path, input_file_name, m, total_thread_number, insertion_edge_vector);


    auto previous_edge_wing_map = make_shared<unordered_map<shared_ptr<abstract_bipartite_edge>, uint32_t>>();
    auto previous_WS = make_shared<unordered_map<shared_ptr<abstract_bipartite_edge>, uint32_t>>();
    uint32_t previous_k_max = 0;
    {
        prepare(G,previous_edge_wing_map,previous_WS,previous_k_max,total_thread_number);
    }

    auto contrastive_edge_wing_map = make_shared<unordered_map<shared_ptr<abstract_bipartite_edge>, uint32_t>>();
    {
        decompose(G,insertion_edge_vector,contrastive_edge_wing_map,total_thread_number);
    }

    /**
     * @brief common parts are pre-calculated
     */
    vector<uint32_t> thread_array{1};
    for (const auto &thread_number :thread_array) {
        wing_insertion_thread_func(G,insertion_edge_vector,previous_edge_wing_map,contrastive_edge_wing_map,previous_WS,previous_k_max,logger,thread_number);
    }
    return 0;
}

