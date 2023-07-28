
#include "core_test/core_test.h"

uint32_t get_exponent(uint32_t  value, uint32_t base){
    uint32_t exponent = 0;
    while(value  / base > 0){
        value = value / base;
        ++exponent;
    }
    return exponent;
}

int main(int argc, char **argv) {
    if(argc < 4){
        std::cout<<"Usage: input path, file name, and thread number!";
    }

    string path = argv[1];
    string input_file_name = argv[2];
    uint32_t  thread_number = std::stoul(argv[3]);

    /**
     *@brief prepare log file
     **/
    string log_path = path + "log/";
    auto directory = std::filesystem::path(log_path);
    if (!exists(directory)) {
        create_directory(directory);
    }

    string log_file_name = log_path  + "decomposition_count_test.log";
    simple_logger logger(log_file_name);

    auto edge_vector = abstract_graph_io::get_edge_vector(path, input_file_name);
    LOG(logger, LOG_RANK::INFO)  << input_file_name << "\n";

    uint32_t  k_max = 0;
    auto graph = shared_ptr<abstract_graph>();
    auto vertex_core_map1 = make_shared<unordered_map<uint32_t,uint32_t>>();
    {
        auto pool = make_shared<thread_pool>(thread_number);
        graph = abstract_graph_io::load_graph(edge_vector,pool);

        simple_timer decomposition_timer;
        k_max = basic_core_decomposition::decompose(graph, vertex_core_map1,thread_number);
        auto decomposition_time = decomposition_timer.get_elapse_second();

        LOG(logger,LOG_RANK::INFO) << "Base decomposition," << decomposition_time << "\n";
    }

    uint32_t  m = 200000;
    auto previous_edge_vector =  make_shared<vector<shared_ptr<abstract_edge>>>();
    auto insertion_edge_vector = make_shared<vector<shared_ptr<abstract_edge>>>();
    auto removal_edge_vector = make_shared<vector<shared_ptr<abstract_edge>>>();
    for(uint32_t i  = 0;i < edge_vector->size();++i){
        if(i < m){
            removal_edge_vector->push_back(edge_vector->at(i));
        }
        if (i < edge_vector->size() - m) {
            previous_edge_vector->push_back(edge_vector->at(i));
        } else {
            insertion_edge_vector->push_back(edge_vector->at(i));
        }
    }

    {
        auto insertion_vertex_set = make_shared<unordered_set<uint32_t>>();
        {
            auto insertion_core_number_vector = make_shared<vector<uint32_t>>(get_exponent(k_max, 2) + 1, 0);
            for(const auto &e:*insertion_edge_vector){
                auto u = e->get_source_vertex_id();
                auto v = e->get_destination_vertex_id();
                insertion_vertex_set->insert(u);
                insertion_vertex_set->insert(v);

                auto min_k = std::min(vertex_core_map1->at(u), vertex_core_map1->at(v));
                ++insertion_core_number_vector->at(get_exponent(min_k, 2));
            }
            ofstream insertion_core_number_file(log_path + "insertion_core_number.csv",std::ios::app);
            insertion_core_number_file << input_file_name << ",";
            for(uint32_t i = 0; i < insertion_core_number_vector->size(); ++i) {
                insertion_core_number_file << i << ":" << insertion_core_number_vector->at(i) << ",";
            }
            insertion_core_number_file<<"\n";
        }

        {
            auto insertion_degree_vector = make_shared<vector<uint32_t>>(get_exponent(graph->get_maximal_degree(), 10)+1, 0);
            for(const auto &v:*insertion_vertex_set){
                auto degree = graph->get_vertex(v)->get_degree();
                ++insertion_degree_vector->at(get_exponent(degree, 10));
            }

            ofstream  insertion_degree_file(log_path+"insertion_degree.csv", std::ios::app);
            insertion_degree_file << input_file_name << ",";
            for(uint32_t i = 0; i< insertion_degree_vector->size();++i){
                insertion_degree_file << i << ":" << insertion_degree_vector->at(i) << ",";
            }
            insertion_degree_file <<"\n";
        }
    }

    {
        auto removal_vertex_set = make_shared<unordered_set<uint32_t>>();
        {
            auto removal_core_number_vector = make_shared<vector<uint32_t>>(get_exponent(k_max, 2) + 1, 0);
            for(const auto &e:*removal_edge_vector){
                auto u = e->get_source_vertex_id();
                auto v = e->get_destination_vertex_id();
                removal_vertex_set->insert(u);
                removal_vertex_set->insert(v);

                auto min_k = std::min(vertex_core_map1->at(u), vertex_core_map1->at(v));
                ++removal_core_number_vector->at(get_exponent(min_k, 2));
            }
            ofstream removal_core_number_file(log_path + "removal_core_number.csv",std::ios::app);
            removal_core_number_file << input_file_name << ",";
            for(uint32_t i = 0; i < removal_core_number_vector->size(); ++i) {
                removal_core_number_file << i << ":" << removal_core_number_vector->at(i) << ",";
            }
            removal_core_number_file << "\n";
        }

        auto removal_degree_vector = make_shared<vector<uint32_t>>(get_exponent(graph->get_maximal_degree(), 10) + 1, 0);
        {
            for(const auto &v:*removal_vertex_set){
                auto degree = graph->get_vertex(v)->get_degree();
                ++removal_degree_vector->at(get_exponent(degree, 10));
            }

            ofstream  removal_degree_file(log_path + "removal_degree.csv",std::ios::app);
            removal_degree_file << input_file_name <<",";
            for(uint32_t i = 0; i < removal_degree_vector->size(); ++i){
                removal_degree_file << i <<":" << removal_degree_vector->at(i) << ",";
            }
            removal_degree_file <<"\n";
        }
    }

    LOG(logger,LOG_RANK::INFO) << "NumberOfVertices," << graph->get_vertex_number()<<"\n";
    LOG(logger,LOG_RANK::INFO) << "NumberOfEdges," << graph->get_edge_number()<<"\n";
    LOG(logger,LOG_RANK::INFO) << "MaximalDegree," << graph->get_maximal_degree()<<"\n";
    LOG(logger,LOG_RANK::INFO) << "kMax," << k_max<<"\n";
}