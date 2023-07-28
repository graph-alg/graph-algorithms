
#include "core_test/core_test.h"


int main(int argc, char **argv) {
    if(argc < 4){
        std::cout<<"Usage: please input path, file name, and thread number!";
    }

    string path = argv[1];
    string input_file_name = argv[2];
    uint32_t thread_number = std::stoul(argv[3]);
    /**
     * @brief prepare the log file
     */
    string log_path = path + "log/";
    string log_file_name = log_path  + "insertion_memory_test.log";
    simple_logger logger(log_file_name);


    auto edge_vector = abstract_graph_io::get_edge_vector(path, input_file_name);
    /**
     * @brief shuffle the edges
     */
//    {
//        auto rd = make_shared<random_device>();
//        shuffle(edge_vector->begin(), edge_vector->end(), *random_generator::get_default_engine(rd));
//    }


    auto previous_edge_vector = make_shared<vector<shared_ptr<abstract_edge>>>();
    auto insertion_edge_vector = make_shared<vector<shared_ptr<abstract_edge>>>();
    /**
     * @brief the number of inserted edges
     */
    uint32_t m = 10000;
    LOG(logger, LOG_RANK::INFO)  << input_file_name << ","<< m << "\n";

    for (auto i = 0; i < edge_vector->size(); ++i) {
        if (i < edge_vector->size() - m) {
            previous_edge_vector->emplace_back(edge_vector->at(i));
        } else {
            insertion_edge_vector->emplace_back(edge_vector->at(i));
        }
    }

    /**
     * @details prepare core number of vertices in previous graph
     */
    auto previous_graph = shared_ptr<abstract_graph>();
    auto previous_vertex_core_map = make_shared<unordered_map<uint32_t, uint32_t>>();
    auto k_order = make_shared<unordered_map<uint32_t, shared_ptr<double_list<uint32_t>>>>();
    auto node_map = make_shared<unordered_map<uint32_t, shared_ptr<double_node<uint32_t>>>>();
    shared_ptr<gadget::Treap> tree;
    shared_ptr<vector<long long>> root;
    {
        auto pool = make_shared<thread_pool>(thread_number);
        previous_graph = abstract_graph_io::load_graph(previous_edge_vector,pool);

        tree = make_shared<gadget::Treap>(2 * previous_graph->get_vertex_number() + 1);
        root = make_shared<vector<long long>>(2 * previous_graph->get_vertex_number() + 1,
                                              2 * previous_graph->get_vertex_number() + 1);

        basic_core_decomposition::decompose(previous_graph, previous_vertex_core_map,
                                            k_order, node_map, tree, root, thread_number);
    }

    /**
     * @brief decomposition for checking results
     */
    auto contrastive_vertex_core_map = make_shared<unordered_map<uint32_t, uint32_t>>();
    {
        previous_graph->insert_edge_collection(insertion_edge_vector);

        simple_timer decomposition_timer;
        basic_core_decomposition::decompose(previous_graph, contrastive_vertex_core_map, thread_number);
        double decomposition_time = decomposition_timer.get_elapse_second();

        LOG(logger, LOG_RANK::INFO) << "Decomposition," << decomposition_time <<"\n";

        /**
         * @brief restore the graph
         */
         previous_graph->remove_edge_collection(insertion_edge_vector);
    }




    /**
     * @brief ordered-based insert maintenance
     */
    {
        auto previous_vertex_core_map2 = container_copy::to_unordered_map<uint32_t,uint32_t>(previous_vertex_core_map);

        auto &k_order2 = k_order;
        auto &node_map2 = node_map;
        auto &tree2 = tree;
        auto &root2 = root;
        auto rem2 = make_shared<unordered_map<uint32_t, uint32_t>>();
        auto ext2 = make_shared<unordered_map<uint32_t, uint32_t>>();
        auto rcd2 = make_shared<unordered_map<uint32_t, shared_ptr<unordered_map<uint32_t, uint32_t>>>>();
        uint32_t n2 = 1;

        simple_timer maintenance_timer;
        order_core_maintenance::initialize(previous_graph, previous_vertex_core_map2, tree2, root2, rem2, ext2,
                                           rcd2, n2);
        for (const auto &e:*insertion_edge_vector) {
            order_core_maintenance::insert(previous_graph, e, previous_vertex_core_map2, k_order2,
                                           node_map2, tree2, root2, rem2, ext2, rcd2, n2);
        }
        auto maintenance_time = maintenance_timer.get_elapse_second();

        if (core_compare::same_associative_map(contrastive_vertex_core_map, previous_vertex_core_map2)) {
            auto memory_size = process_information::get_memory(getpid());
            LOG(logger, LOG_RANK::INFO) << "Order Insertion," << maintenance_time<<","<< double (memory_size)/1024/1024 << "\n";
        }

        /**
         * @brief restore the graph
         */
        previous_graph->remove_edge_collection(insertion_edge_vector);
    }
}
