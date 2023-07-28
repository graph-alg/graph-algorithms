
#include "bipartite_core_test/bipartite_core_test.h"

int main(int argc, char **argv) {
    if(argc < 3){
        std::cout<<"Usage: input path, and thread number!";
    }

    const string input_path = argv[1];
    const uint32_t thread_number = std::stoul(argv[2]);

    auto rd = make_shared<random_device>();
    auto directory = path(input_path);
    auto pool = make_shared<thread_pool>(thread_number);
    for (auto &file_iter: std::filesystem::directory_iterator(input_path)) {
        if (!std::filesystem::is_regular_file(file_iter)) {
            continue;
        }
        pool->submit_task([=] {
            const auto &file = file_iter.path();
            auto input_file_name = file.filename().string();
            auto edge_vector = abstract_bipartite_graph_io::get_edge_vector(input_path, input_file_name);

            auto half_size = uint32_t (std::ceil(edge_vector->size() / 2));
            auto half_begin = edge_vector->begin();
            std::advance(half_begin, half_size);

            {
                auto sub_edge_vector1 = make_shared<vector<shared_ptr<abstract_bipartite_edge>>>();
                copy(edge_vector->begin(), half_begin, inserter(*sub_edge_vector1, sub_edge_vector1->end()));

                abstract_bipartite_graph_io::output_csv_file(sub_edge_vector1, input_path + input_file_name+"-half1");
            }

            {
                auto sub_edge_vector2 = make_shared<vector<shared_ptr<abstract_bipartite_edge>>>();
                copy(half_begin, edge_vector->end(), inserter(*sub_edge_vector2, sub_edge_vector2->end()));

                abstract_bipartite_graph_io::output_csv_file(sub_edge_vector2, input_path + input_file_name+"-half2");
            }
        });
    }
    pool->barrier();
}

