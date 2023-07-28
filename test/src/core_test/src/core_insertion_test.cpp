
#include "core_test/core_test.h"


int main(int argc, char **argv) {
    if(argc < 4){
        std::cout<<"Usage: please input path, file name, and thread number!";
    }
    string path = argv[1];
    string input_file_name = argv[2];
    auto thread_number = std::stoul(argv[3]);

    /**
     * @brief prepare the log file
     */
    string log_path = path + "log/";
    string log_file_name = log_path + "insertion_test.log";
    simple_logger logger(log_file_name);


    auto edge_vector = abstract_graph_io::get_edge_vector(path, input_file_name);
    /**
     * @brief shuffle the edges, test!
     */
//    {
//        auto rd = make_shared<random_device>();
//        shuffle(edge_vector->begin(),edge_vector->end(),*random_generator::get_default_engine(rd));
//        for(const auto &e:*edge_vector){
//            printf("%lu,%lu\n",e->get_source_vertex_id(), e->get_destination_vertex_id());
//        }
//    }


    auto previous_edge_vector = make_shared<vector<shared_ptr<abstract_edge>>>();
    auto insertion_edge_vector = make_shared<vector<shared_ptr<abstract_edge>>>();
    /**
     * @brief the number of inserted edges
     */
    uint32_t m = 10000;
    LOG(logger, LOG_RANK::INFO)  << input_file_name << "," << m << "\n";

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
    auto tree = make_shared<gadget::Treap>(2 * previous_graph->get_vertex_number() + 1);
    auto root = make_shared<vector<long long>>(2 * previous_graph->get_vertex_number() + 1,
                                          2 * previous_graph->get_vertex_number() + 1);
    uint32_t  previous_k_max_value = 0;
    {
        auto pool = make_shared<thread_pool>(thread_number);
        previous_graph = abstract_graph_io::load_graph(previous_edge_vector, pool);
        previous_k_max_value = basic_core_decomposition::decompose(previous_graph, previous_vertex_core_map,
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
        auto decomposition_time = decomposition_timer.get_elapse_second();

        LOG(logger, LOG_RANK::INFO) << "Decomposition," << decomposition_time << "\n";

        /**
         * @brief restore the graph
         */
        previous_graph->remove_edge_collection(insertion_edge_vector);
    }

    /**
      * @brief traversal insert
      */
    {
        auto previous_vertex_core_map1 = container_copy::to_unordered_map<uint32_t,uint32_t>(previous_vertex_core_map);
        auto rcd1 = make_shared<unordered_map<uint32_t, shared_ptr<unordered_map<uint32_t, uint32_t>>>>();
        uint32_t n2 = 2;

        simple_timer maintenance_timer;
        traversal_core_maintenance::init_rcd(previous_graph, previous_vertex_core_map1, rcd1, n2);
        for (const auto &edge:*insertion_edge_vector) {
            traversal_core_maintenance::insert(previous_graph, previous_vertex_core_map1, rcd1, n2,
                                               edge);
        }
        auto maintenance_time = maintenance_timer.get_elapse_second();

        if (core_compare::same_associative_map(contrastive_vertex_core_map, previous_vertex_core_map1)) {
            LOG(logger, LOG_RANK::INFO) << "Traversal Insertion," << maintenance_time << "\n";
        }

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
            LOG(logger, LOG_RANK::INFO) << "Order Insertion," << maintenance_time << "\n";
        }

        /**
         * @brief restore the graph
         */
        previous_graph->remove_edge_collection(insertion_edge_vector);
    }

    /**
     * @brief joint insert
     */
    {
        auto previous_vertex_core_map3 = container_copy::to_unordered_map<uint32_t,uint32_t>(previous_vertex_core_map);
        auto EI3 = container_convert::to_unordered_set<shared_ptr<abstract_edge>>(insertion_edge_vector);

        simple_timer maintenance_timer;
        jes_core_maintenance::joint_insert(previous_graph, EI3, previous_vertex_core_map3, thread_number);
        auto maintenance_time = maintenance_timer.get_elapse_second();
        if (core_compare::same_associative_map(contrastive_vertex_core_map, previous_vertex_core_map3)) {
            LOG(logger, LOG_RANK::INFO) << "Joint Insertion," << maintenance_time << "\n";
        }
        /**
         * @brief restore the graph
         */
        previous_graph->remove_edge_collection(insertion_edge_vector);
    }


    /**
     * @brief quasi insert
     */
    {
        auto previous_vertex_core_map4 = container_copy::to_unordered_map<uint32_t,uint32_t>(previous_vertex_core_map);
        auto insertion_edge_set4 = container_convert::to_unordered_set<shared_ptr<abstract_edge>>(
                insertion_edge_vector);
        auto previous_k_max4 = make_shared<uint32_t>(previous_k_max_value);

        simple_timer maintenance_timer;
        quasi_core_maintenance::insert(previous_graph, insertion_edge_set4,
                                       previous_vertex_core_map4, previous_k_max4);
        auto maintenance_time = maintenance_timer.get_elapse_second();
        if (core_compare::same_associative_map(previous_vertex_core_map4, contrastive_vertex_core_map)) {
            LOG(logger, LOG_RANK::INFO) << "Quasi Insertion," << maintenance_time << "\n";
        }

        /**
         * @brief restore the graph
         */
        previous_graph->remove_edge_collection(insertion_edge_vector);
    }

    /**
     * @brief parallel insert
     */
    {
        auto insertion_edge_set5 = container_convert::to_unordered_set<shared_ptr<abstract_edge>>(
                insertion_edge_vector);
        auto previous_vertex_core_map5 = container_copy::to_unordered_map<uint32_t,uint32_t>(previous_vertex_core_map);
        auto previous_k_max5 = make_shared<uint32_t>(previous_k_max_value);

        simple_timer maintenance_timer;
        quasi_core_maintenance::insert(previous_graph, insertion_edge_set5,
                                       previous_vertex_core_map5, previous_k_max5,
                                       thread_number);
        auto maintenance_time = maintenance_timer.get_elapse_second();
        if (core_compare::same_associative_map(previous_vertex_core_map5, contrastive_vertex_core_map)) {
            LOG(logger, LOG_RANK::INFO) << "Parallel Insertion," << maintenance_time << "\n";
        }

        /**
         * @brief restore the graph
         */
        previous_graph->remove_edge_collection(insertion_edge_vector);
    }
}
