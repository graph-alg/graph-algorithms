
#include "wing_test/wing_test.h"
#include <iostream>

int main(int argc, char **argv) {
    std::string path = argv[1];
    std::string input_file_name = argv[2];
    uint32_t thread_number = std::stoul(argv[3]);

    if(argc < 3){
        std::cout<<"Usage: require  path, filename and thread number!"<<std::endl;
    }

    std::string log_path = path + "log/";
    std::string log_file_name = log_path + input_file_name + "_insertion_test.log";
    auto directory = std::filesystem::path(log_path.c_str());
    if (!exists(directory)) {
        create_directory(directory);
    }

    simple_logger logger(log_file_name);

    auto edge_vector = abstract_bipartite_graph_io::get_edge_vector(path, input_file_name);
    {
//        auto rd = make_shared<random_device>();
//        std::shuffle(edge_vector->begin(), edge_vector->end(), *random_generator::get_default_engine(rd));
    }

    /**
     * @brief the number of inserted edges
     */
    auto m = 10000;

    auto previous_edge_vector = make_shared<vector<shared_ptr<abstract_bipartite_edge>>>();
    auto insertion_edge_vector = make_shared<vector<shared_ptr<abstract_bipartite_edge>>>();
    for (uint32_t i = 0; i < edge_vector->size(); i++) {
        if (i < edge_vector->size() - m) {
            previous_edge_vector->push_back(edge_vector->at(i));
        } else {
            insertion_edge_vector->push_back(edge_vector->at(i));
        }
    }

    auto G = shared_ptr<abstract_bipartite_graph>();
    {
        auto pool = make_shared<thread_pool>(thread_number);
        G = abstract_bipartite_graph_io::load_graph(previous_edge_vector, pool);
    }


    auto previous_edge_wing_map = make_shared<unordered_map<shared_ptr<abstract_bipartite_edge>, uint32_t>>();
    auto previous_wing_order_map = make_shared<unordered_map<uint32_t, shared_ptr<extend_list<double, shared_ptr<abstract_bipartite_edge>>>>>();
    auto previous_WS = make_shared<unordered_map<shared_ptr<abstract_bipartite_edge>, uint32_t>>();
    auto previous_rem = make_shared<unordered_map<shared_ptr<abstract_bipartite_edge>, uint32_t>>();
    auto previous_WL = make_shared<unordered_map<uint32_t,shared_ptr<map<uint32_t,shared_ptr<unordered_set<uint32_t>>>>>>();


    uint32_t previous_k_max = 0;
    {
        auto pool = make_shared<thread_pool>(thread_number);
        auto edge_mutex_map =
                make_shared < unordered_map < shared_ptr < abstract_bipartite_edge >, shared_ptr<mutex>>>();
        auto edge_rank_map = make_shared < unordered_map < shared_ptr < abstract_bipartite_edge >, uint32_t>>();
        auto edge_support_map = make_shared < unordered_map < shared_ptr < abstract_bipartite_edge >, uint32_t>>();
        basic_wing_decomposition::init(G, edge_mutex_map, edge_rank_map, edge_support_map, previous_edge_wing_map, previous_WS, previous_rem, pool);
        previous_k_max = basic_wing_decomposition::decompose(G, edge_mutex_map, edge_rank_map, edge_support_map, previous_edge_wing_map,previous_WS, previous_wing_order_map,
                                                             previous_rem,  pool);
    }

    {
        auto WS = make_shared<unordered_map<shared_ptr<abstract_bipartite_edge>, uint32_t>>();
        auto edge_mutex_map = make_shared<unordered_map<shared_ptr<abstract_bipartite_edge>, shared_ptr<mutex>>>();
        auto edge_rank_map = make_shared<unordered_map<shared_ptr<abstract_bipartite_edge>, uint32_t>>();
        auto edge_support_map = make_shared<unordered_map<shared_ptr<abstract_bipartite_edge>, uint32_t>>();
        auto pool = make_shared<thread_pool>(thread_number);
        index_wing_decomposition::init(G, edge_mutex_map, edge_rank_map, edge_support_map,
                                       previous_edge_wing_map, pool);
        /*index_wing_decomposition::decompose(G, edge_mutex_map, edge_rank_map, edge_support_map,
                                            previous_edge_wing_map, pool);*/

        index_wing_decomposition::WS_WL_init(G, edge_mutex_map, edge_rank_map, previous_edge_wing_map, WS, previous_WL, pool);
    }

    auto contrastive_edge_wing_map =
            make_shared < unordered_map < shared_ptr < abstract_bipartite_edge >, uint32_t>>();
    {
        G->insert_edge_collection(insertion_edge_vector);

        simple_timer bottom_up_timer;
        auto edge_mutex_map =
                make_shared < unordered_map < shared_ptr < abstract_bipartite_edge >, shared_ptr<mutex>>>();
        auto edge_rank_map = make_shared < unordered_map < shared_ptr < abstract_bipartite_edge >, uint32_t>>();
        auto edge_support_map = make_shared < unordered_map < shared_ptr < abstract_bipartite_edge >, uint32_t>>();

        auto pool = make_shared<thread_pool>(thread_number);
        index_wing_decomposition::init(G, edge_mutex_map, edge_rank_map, edge_support_map,
                                       contrastive_edge_wing_map, pool);
        index_wing_decomposition::decompose(G, edge_mutex_map, edge_rank_map, edge_support_map,
                                            contrastive_edge_wing_map, pool);
        auto bottom_up_time = bottom_up_timer.get_elapse_second();

        LOG(logger, LOG_RANK::INFO) << "Bottom-Up Decomposition," << bottom_up_time << '\n';

        G->insert_edge_collection(insertion_edge_vector);
    }

    /**
     * @brief baseline
     */
    {
        auto insertion_edge_set = container_copy::to_unordered_set<shared_ptr<abstract_bipartite_edge>>(
                insertion_edge_vector);
        auto order_edge_wing_map = container_copy::to_unordered_map<shared_ptr<abstract_bipartite_edge>, uint32_t>(
                previous_edge_wing_map);

        auto wing_order_map = make_shared<unordered_map<uint32_t,shared_ptr<extend_list<double,shared_ptr<abstract_bipartite_edge>>>>>();

        for(const auto&[k,k_list]:*previous_wing_order_map){
            wing_order_map->insert({k, make_shared<extend_list<double,shared_ptr<abstract_bipartite_edge>>>(k_list)});
        }

        auto ext = make_shared<unordered_map<shared_ptr<abstract_bipartite_edge>, uint32_t>>();
        auto rem = container_copy::to_unordered_map<shared_ptr<abstract_bipartite_edge>, uint32_t>(
                previous_rem);
        auto ts = container_copy::to_unordered_map<shared_ptr<abstract_bipartite_edge>, uint32_t>(
                previous_WS);
        simple_timer order_timer;
        order_wing_maintenance::init(order_edge_wing_map, ext);
        order_wing_maintenance::insert(G, insertion_edge_set, order_edge_wing_map, ts,
                                       wing_order_map, rem, ext);
        auto order_time = order_timer.get_elapse_second();

        if (wing_compare::same_associative_map(order_edge_wing_map, contrastive_edge_wing_map)) {
            LOG(logger, LOG_RANK::INFO) << "Order Insertion," << order_time << "\n";
        }

        G->remove_edge_collection(insertion_edge_vector);
    }

    {
        auto insertion_edge_set = container_copy::to_unordered_set<shared_ptr<abstract_bipartite_edge>>(
                insertion_edge_vector);
        auto quasi_edge_wing_map = container_copy::to_unordered_map<shared_ptr<abstract_bipartite_edge>, uint32_t>(
                previous_edge_wing_map);
        auto previous_k_max3 = make_shared<uint32_t>(previous_k_max);
        auto WS = make_shared<unordered_map<shared_ptr<abstract_bipartite_edge>, uint32_t>>();

        simple_timer improved_timer;
        auto edge_mutex_map = make_shared<unordered_map<shared_ptr<abstract_bipartite_edge>, shared_ptr<mutex>>>();
        auto edge_rank_map = make_shared< unordered_map<shared_ptr<abstract_bipartite_edge>, uint32_t>>();
        auto edge_support_map = make_shared<unordered_map<shared_ptr<abstract_bipartite_edge>, uint32_t>>();
        auto pool = make_shared<thread_pool>(thread_number);
        quasi_wing_maintenance::init(quasi_edge_wing_map, edge_mutex_map, edge_rank_map, edge_support_map, WS, pool);
        quasi_wing_maintenance::insert(G, edge_mutex_map, insertion_edge_set, edge_rank_map, edge_support_map, quasi_edge_wing_map,
                                       WS,
                                       previous_k_max3,
                                       pool);
        auto improved_time = improved_timer.get_elapse_second();

        if (wing_compare::same_associative_map(quasi_edge_wing_map, contrastive_edge_wing_map)) {
            LOG(logger, LOG_RANK::INFO) << "Basic Insertion," << improved_time << "\n";
        }

        G->remove_edge_collection(insertion_edge_vector);
    }


    {
        auto WS = container_copy::to_unordered_map<shared_ptr<abstract_bipartite_edge>, uint32_t>(previous_WS);
        auto insertion_edge_set = container_copy::to_unordered_set<shared_ptr<abstract_bipartite_edge>>(
                insertion_edge_vector);
        auto quasi_edge_wing_map = container_copy::to_unordered_map<shared_ptr<abstract_bipartite_edge>, uint32_t>(
                previous_edge_wing_map);
        auto previous_k_max3 = make_shared<uint32_t>(previous_k_max);

        simple_timer improved_timer;
        auto edge_mutex_map = make_shared<unordered_map<shared_ptr<abstract_bipartite_edge>, shared_ptr<mutex>>>();
        auto edge_rank_map = make_shared< unordered_map<shared_ptr<abstract_bipartite_edge>, uint32_t>>();
        auto edge_support_map = make_shared<unordered_map<shared_ptr<abstract_bipartite_edge>, uint32_t>>();
        auto pool = make_shared<thread_pool>(thread_number);
        quasi_wing_maintenance::init(quasi_edge_wing_map, edge_mutex_map, edge_rank_map, edge_support_map, pool);
        quasi_wing_maintenance::insert2(G, edge_mutex_map, insertion_edge_set, edge_rank_map, edge_support_map, quasi_edge_wing_map,
                                        WS,
                                        previous_k_max3,
                                        pool);
        auto improved_time = improved_timer.get_elapse_second();

        if (wing_compare::same_associative_map(quasi_edge_wing_map, contrastive_edge_wing_map)) {
            LOG(logger, LOG_RANK::INFO) << "WS Insertion," << improved_time << "\n";
        }

        G->remove_edge_collection(insertion_edge_vector);
    }

//    {
//        auto WS = container_copy::to_unordered_map<shared_ptr<abstract_bipartite_edge>, uint32_t>(previous_WS);
//        auto insertion_edge_set = container_copy::to_unordered_set<shared_ptr<abstract_bipartite_edge>>(
//                insertion_edge_vector);
//        auto bloom_edge_wing_map = container_copy::to_unordered_map<shared_ptr<abstract_bipartite_edge>, uint32_t>(
//                previous_edge_wing_map);
//        auto k_max = make_shared<uint32_t>(previous_k_max);
//
//
//        simple_timer improved_timer;
//        auto edge_mutex_map = make_shared<unordered_map<shared_ptr<abstract_bipartite_edge>, shared_ptr<mutex>>>();
//        auto edge_rank_map = make_shared< unordered_map<shared_ptr<abstract_bipartite_edge>, uint32_t>>();
//        auto edge_support_map = make_shared<unordered_map<shared_ptr<abstract_bipartite_edge>, uint32_t>>();
//        auto pool = make_shared<thread_pool>(thread_number);
//        quasi_wing_maintenance::init(bloom_edge_wing_map, edge_mutex_map, edge_rank_map, edge_support_map, pool);
//        quasi_wing_maintenance::insert3(G, edge_mutex_map, insertion_edge_set, edge_rank_map, edge_support_map, bloom_edge_wing_map,
//                                        WS,
//                                        k_max,
//                                        pool);
//        auto improved_time = improved_timer.get_elapse_second();
//
//        if (wing_compare::same_associative_map(bloom_edge_wing_map, contrastive_edge_wing_map)) {
//            LOG(logger, LOG_RANK::INFO) << "Bloom Insertion," << improved_time << "\n";
//        }
//
//        G->remove_edge_collection(insertion_edge_vector);
//    }
    {
        auto WS = container_copy::to_unordered_map<shared_ptr<abstract_bipartite_edge>, uint32_t>(previous_WS);
        auto WL = make_shared<unordered_map<uint32_t,shared_ptr<map<uint32_t,shared_ptr<unordered_set<uint32_t>>>>>>();
        {
            for(const auto&[source,k_map]:*previous_WL) {
                WL->insert({source, shared_ptr<map<uint32_t, shared_ptr<unordered_set<uint32_t>>>>()});
            }
            auto pool = make_shared<thread_pool>(thread_number);
            auto location_vector = pool->split_task(previous_WL);
            for(uint32_t i = 0; i < thread_number; ++i){
                pool->submit_task([=]{
                    for(auto iter = *location_vector->at(i); iter!=*location_vector->at(i+1); ++iter){
                        auto &[source,k_map] = *iter;
                        WL->at(source) = make_shared<map<uint32_t, shared_ptr<unordered_set<uint32_t>>>>();
                        for(const auto&[wing_number, v_set]:*k_map){
                            WL->at(source)->insert({wing_number, make_shared<unordered_set<uint32_t>>(*v_set)});
                        }
                    }
                });
            }
            pool->barrier();
        }
        auto wl_edge_wing_map = container_copy::to_unordered_map<shared_ptr<abstract_bipartite_edge>, uint32_t>(
                previous_edge_wing_map);
        auto insertion_edge_set = container_copy::to_unordered_set<shared_ptr<abstract_bipartite_edge>>(
                insertion_edge_vector);
        auto k_max = make_shared<uint32_t>(previous_k_max);
        auto edge_mutex_map = make_shared<unordered_map<shared_ptr<abstract_bipartite_edge>, shared_ptr<mutex>>>();
        auto edge_rank_map = make_shared< unordered_map<shared_ptr<abstract_bipartite_edge>, uint32_t>>();
        auto edge_support_map = make_shared<unordered_map<shared_ptr<abstract_bipartite_edge>, uint32_t>>();
        auto pool = make_shared<thread_pool>(thread_number);

        simple_timer wl_timer;
        quasi_wing_maintenance::init(wl_edge_wing_map, edge_mutex_map, edge_rank_map, edge_support_map, pool);
        quasi_wing_maintenance::insert4(G, edge_mutex_map, insertion_edge_set, edge_rank_map, edge_support_map, wl_edge_wing_map,
                                        WS,WL,
                                        k_max,
                                        pool);
        auto wl_time = wl_timer.get_elapse_second();
        if (wing_compare::same_associative_map(wl_edge_wing_map, contrastive_edge_wing_map)) {
            LOG(logger, LOG_RANK::INFO) << "WL Insertion," << wl_time << "\n";
        }
        G->remove_edge_collection(insertion_edge_vector);

    }
}




