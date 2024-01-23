#include "multiple_core_test/multiple_core_test.h"

void information_test(const shared_ptr<temporal_graph>& G,
                      const shared_ptr<multiple_core>& mc,
                      const string& log_path,
                      const string& input_file_name){
    const string log_file_name = log_path + "graph_information_test.log";
    scnu::simple_logger logger(log_file_name);

    LOG(logger,LOG_RANK::INFO) << input_file_name << "\n";

    LOG(logger,LOG_RANK::INFO) << "Vertices," << G->get_vertex_size() << "\n";
    LOG(logger,LOG_RANK::INFO) << "Edges," << G->get_edge_size() << "\n";
    LOG(logger,LOG_RANK::INFO) << "Maximal Neighbor Size," << G->get_maximal_neighbor_vertex_size() << "\n";
    LOG(logger,LOG_RANK::INFO) << "Maximal Parallel Edge Size," << G->get_maximal_parallel_edge_size() << "\n";
    LOG(logger,LOG_RANK::INFO) << "Average Neighbor Size," << G->get_average_neighbor_vertex_size() << "\n";
    LOG(logger,LOG_RANK::INFO) << "Average Edge Size," << G->get_average_edge_size() << "\n";
    LOG(logger,LOG_RANK::INFO) << "KMax," << mc->get_maximal_k() << "\n";
    LOG(logger,LOG_RANK::INFO) << "Delta," << mc->get_maximal_delta() << "\n";
    LOG(logger,LOG_RANK::INFO) << "Core Size," << mc->get_core_pair_set(1)->size() << "\n";
}

void index_memory_test(const shared_ptr<multiple_core>& mc,
                       const shared_ptr<unordered_map<uint32_t, shared_ptr<multiple_core_number_set_index>>> &vertex_core_number_set_map,
                       const shared_ptr<unordered_map<uint32_t, shared_ptr<multiple_core_pair_set_index>>> &vertex_core_pair_set_map,
                       const shared_ptr<unordered_map<uint32_t, shared_ptr<multiple_core_pair_map_index>>> &vertex_core_pair_map_map,
                       const string& log_path,
                       const string& input_file_name){
    const string log_file_name = log_path  + "index_memory_test.log";
    scnu::simple_logger logger(log_file_name);

    LOG(logger,LOG_RANK::INFO) << input_file_name << "\n";

    {
        simple_timer creation_timer;
        mc->convert(vertex_core_number_set_map);
        auto creation_time = creation_timer.get_elapse_second();

        LOG(logger,LOG_RANK::INFO) << "Core Number Set Creation," << creation_time<< "\n";

        uint32_t memory_cost = sizeof(vertex_core_number_set_map) +
                               sizeof (*vertex_core_number_set_map->begin()) * vertex_core_number_set_map->size();
        for(const auto&[v,v_index]:*vertex_core_number_set_map)
        {
            memory_cost += v_index->get_memory_cost();
        }
        LOG(logger,LOG_RANK::INFO) << "Core Number Set Memory Cost (MB)," << double (memory_cost)/1024/1024<< "\n";
    }

    {
        simple_timer creation_timer;
        mc->convert(vertex_core_pair_set_map);
        auto creation_time = creation_timer.get_elapse_second();

        LOG(logger,LOG_RANK::INFO) << "Core Pair Set Creation," << creation_time<< "\n";


        uint32_t memory_cost = sizeof(vertex_core_pair_set_map) +
                               sizeof (*vertex_core_pair_set_map->begin()) * vertex_core_pair_set_map->size();
        for(const auto&[v,v_index]:*vertex_core_pair_set_map)
        {
            memory_cost += v_index->get_memory_cost();
        }
        LOG(logger,LOG_RANK::INFO) << "Core Pair Set Memory Cost (MB)," << double (memory_cost)/1024/1024<< "\n";
    }


    {

        simple_timer creation_timer;
        mc->convert(vertex_core_pair_map_map);
        auto creation_time = creation_timer.get_elapse_second();

        LOG(logger,LOG_RANK::INFO) << "Core Pair Map Creation," << creation_time<< "\n";

        uint32_t memory_cost = sizeof(vertex_core_pair_map_map) +
                               sizeof (*vertex_core_pair_map_map->begin()) * vertex_core_pair_map_map->size();
        for(const auto&[v,v_index]:*vertex_core_pair_map_map)
        {
            memory_cost += v_index->get_memory_cost();
        }
        LOG(logger,LOG_RANK::INFO) << "Core Pair Map Memory Cost (MB)," << double (memory_cost)/1024/1024<< "\n";
    }

}

void query_test(const shared_ptr<multiple_core>& mc,
                const shared_ptr<unordered_map<uint32_t, shared_ptr<multiple_core_number_set_index>>> &vertex_core_number_set_map,
                const shared_ptr<unordered_map<uint32_t, shared_ptr<multiple_core_pair_set_index>>> &vertex_core_pair_set_map,
                const shared_ptr<unordered_map<uint32_t, shared_ptr<multiple_core_pair_map_index>>> &vertex_core_pair_map_map,
                const string& log_path,
                const string& input_file_name){
    const string log_file_name = log_path + "query_test.log";
    scnu::simple_logger logger(log_file_name);

    LOG(logger,LOG_RANK::INFO) << input_file_name << "\n";

    auto core_pair_vector = container_copy::to_vector<pair<uint32_t, uint32_t>>(mc->get_core_pair_set(1));
    auto rd = make_shared<random_device>();
    std::shuffle(core_pair_vector->begin(),core_pair_vector->end(),*random_generator::get_default_engine(rd));

    {
        auto core_pair_set = make_shared<unordered_set<pair<uint32_t, uint32_t>, hash_pair>>();
        for(uint32_t i = 0; i < 60; i++) {
            core_pair_set->insert(core_pair_vector->at(i));
        }

        auto base_core_pair_map =make_shared<unordered_map<pair<uint32_t, uint32_t>, shared_ptr<unordered_set<uint32_t>>, hash_pair>>();
        {
            simple_timer query_timer;
            multiple_core_query::query_multiple_core(vertex_core_number_set_map, core_pair_set, base_core_pair_map);
            auto query_time = query_timer.get_elapse_second();

            LOG(logger, LOG_RANK::INFO) << "Core Number Set Query," << query_time << "\n";
        }


        {
            auto sub_core_pair_map =make_shared<unordered_map<pair<uint32_t, uint32_t>, shared_ptr<unordered_set<uint32_t>>, hash_pair>>();

            simple_timer query_timer;
            multiple_core_query::query_multiple_core(vertex_core_pair_set_map, core_pair_set, sub_core_pair_map);
            auto query_time = query_timer.get_elapse_second();

            if(multiple_core_compare::same(base_core_pair_map, sub_core_pair_map)){
                LOG(logger,LOG_RANK::INFO) << "Core Pair Set Query," << query_time << "\n";
            }
        }


        {
            auto sub_core_pair_map =make_shared<unordered_map<pair<uint32_t, uint32_t>, shared_ptr<unordered_set<uint32_t>>, hash_pair>>();

            simple_timer query_timer;
            multiple_core_query::query_multiple_core(vertex_core_pair_map_map, core_pair_set, sub_core_pair_map);
            auto query_time = query_timer.get_elapse_second();

            if(multiple_core_compare::same(base_core_pair_map, sub_core_pair_map)){
                LOG(logger,LOG_RANK::INFO) << "Core Pair Map Query," << query_time << "\n";
            }
        }
    }
}

int main(int argc, char ** argv){
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
    * @brief load global edge vector and shuffle it
    */
    auto edge_vector = temporal_graph_io::get_edge_vector(path, input_file_name);



    auto pool = make_shared<thread_pool>(thread_number);
    auto metric_map = make_shared<unordered_map<uint32_t,shared_ptr<multiple_core_pair_map_index>>>();
    auto G = temporal_graph_io::load_graph(edge_vector, pool);
    {
        auto vertex_edge_index_map = make_shared<unordered_map<uint32_t, shared_ptr<map<uint32_t, shared_ptr<unordered_set<uint32_t>>>>>>();
        auto vertex_mutex_map = make_shared<unordered_map<uint32_t, shared_ptr<mutex>>>();
        branch_multiple_core_decomposition::init(G, vertex_mutex_map, vertex_edge_index_map, metric_map, pool);
        branch_multiple_core_decomposition::decompose(G, vertex_mutex_map, vertex_edge_index_map, metric_map, pool);
    }
    auto mc = make_shared<multiple_core>(metric_map);

    auto vertex_core_number_set_map = make_shared<unordered_map<uint32_t, shared_ptr<multiple_core_number_set_index>>>();
    auto vertex_core_pair_set_map = make_shared<unordered_map<uint32_t, shared_ptr<multiple_core_pair_set_index>>>();
    auto vertex_core_pair_map_map = make_shared<unordered_map<uint32_t, shared_ptr<multiple_core_pair_map_index>>>();


    information_test(G, mc, log_path, input_file_name);

    index_memory_test(mc, vertex_core_number_set_map, vertex_core_pair_set_map, vertex_core_pair_map_map, log_path, input_file_name);

    query_test(mc, vertex_core_number_set_map, vertex_core_pair_set_map, vertex_core_pair_map_map, log_path, input_file_name);
    return 0;
}

