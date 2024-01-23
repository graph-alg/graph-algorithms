
#include "wing_test/wing_test.h"

int main(int argc, char **argv)
{
    if(argc < 3)
    {
        cout<<"Usage: Please input the path and the file name!\n";
    }

    string path = argv[1];
    string input_file_name = argv[2];
    uint32_t thread_number = std::stoul(argv[3]);

    string log_path = path + "log/";
    string log_file_name = log_path + input_file_name + "_decomposition_test.log";
    auto directory = std::filesystem::path(log_path);
    if (!exists(directory))
    {
        create_directory(directory);
    }
    scnu::simple_logger logger(log_file_name);

    auto edge_vector = scnu::abstract_bipartite_graph_io::get_edge_vector(path, input_file_name);
    /**
     * @brief shuffle all edges, test!!!
     */
//    {
//        auto rd = make_shared<random_device>();
//        std::shuffle(edge_vector->begin(),edge_vector->end(),*random_generator::get_default_engine(rd));
//    }

    auto G = make_shared<abstract_bipartite_graph>();
    uint32_t k_max = 0;
    {
        auto pool = make_shared<thread_pool>(thread_number);
        G =  abstract_bipartite_graph_io::load_graph(edge_vector, pool);
    }


//    auto basic_edge_wing_map = make_shared<unordered_map<shared_ptr<scnu::abstract_bipartite_edge>,uint32_t>>();
//    {
//        simple_timer basic_timer;
//        auto edge_mutex_map = make_shared<unordered_map<shared_ptr<abstract_bipartite_edge>, shared_ptr<mutex>>>();
//        auto edge_rank_map = make_shared<unordered_map<shared_ptr<abstract_bipartite_edge>, uint32_t>>();
//        auto edge_support_map = make_shared<unordered_map<shared_ptr<abstract_bipartite_edge>, uint32_t>>();
//        auto pool = make_shared<thread_pool>(thread_number);
//
//        basic_wing_decomposition::init(G,edge_mutex_map, edge_rank_map, edge_support_map, basic_edge_wing_map, pool);
//        k_max = basic_wing_decomposition::decompose(G, edge_mutex_map, edge_rank_map, edge_support_map, basic_edge_wing_map, pool);
//
//        auto basic_time = basic_timer.get_elapse_second();
//        LOG(logger,LOG_RANK::INFO)<< "Basic Decomposition,"<<basic_time<<"\n";
//    }


    {
        auto index_edge_wing_map = make_shared<unordered_map<shared_ptr<scnu::abstract_bipartite_edge>,uint32_t>>();

        simple_timer bottom_up_timer;
        auto edge_mutex_map = make_shared<unordered_map<shared_ptr<abstract_bipartite_edge>, shared_ptr<mutex>>>();
        auto edge_rank_map = make_shared<unordered_map<shared_ptr<abstract_bipartite_edge>, uint32_t>>();
        auto edge_support_map = make_shared<unordered_map<shared_ptr<abstract_bipartite_edge>, uint32_t>>();
        auto edge_wing_map = make_shared<unordered_map<shared_ptr<abstract_bipartite_edge>, uint32_t>>();
        auto pool = make_shared<thread_pool>(thread_number);

        index_wing_decomposition::init(G, edge_mutex_map, edge_rank_map, edge_support_map, edge_wing_map, pool);
        index_wing_decomposition::decompose(G, edge_mutex_map, edge_rank_map, edge_support_map, index_edge_wing_map,
                                            pool);
        auto bottom_time = bottom_up_timer.get_elapse_second();

        for(const auto &[e, wing_number]:*index_edge_wing_map){
            if(wing_number > k_max){
                k_max = wing_number;
            }
        }

//        if (scnu::wing_compare::same_associative_map(index_edge_wing_map, basic_edge_wing_map))
//        {
            LOG(logger,LOG_RANK::INFO)<< "Bottom-Up Decomposition,"<<bottom_time << "\n";
 //       }
    }


    /**
     * @brief output the basic information of graphs
     */
    LOG(logger,LOG_RANK::INFO) << "LeftVertices:" << G->get_left_vertex_number() << '\n';
    LOG(logger,LOG_RANK::INFO) << "RightVertices:" << G->get_right_vertex_number() << '\n';
    LOG(logger,LOG_RANK::INFO) << "Edges:" << G->get_edge_number() << '\n';
    LOG(logger,LOG_RANK::INFO) << "MaxLeftDegree:" << G->get_maximal_left_vertex_degree() << '\n';
    LOG(logger,LOG_RANK::INFO) << "MaxRightDegree:" << G->get_maximal_right_vertex_degree() << '\n';
    LOG(logger,LOG_RANK::INFO) << "MaxWingNumber:" << k_max<<'\n';

    return 0;
}
