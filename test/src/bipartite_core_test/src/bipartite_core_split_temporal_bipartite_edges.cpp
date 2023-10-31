
#include "bipartite_core_test/bipartite_core_test.h"

int main(int argc, char **argv) {
    if (argc < 3) {
        std::cout << "Usage: input path, and thread number!";
    }

    const string input_path = argv[1];
    const double rate = std::stod(argv[2]);
    const uint32_t thread_number = std::stoul(argv[3]);


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

            auto size = uint32_t(std::ceil(edge_vector->size() * rate));

            {
                auto sub_edge_vector = make_shared<vector<shared_ptr<abstract_bipartite_edge>>>(size);

                for (uint32_t i = 0; i < size; ++i) {
                    sub_edge_vector->at(i) = edge_vector->at(i);
                }

                abstract_bipartite_graph_io::output_csv_file(sub_edge_vector, input_path + input_file_name + "-head");
            }

            {
                auto sub_edge_vector = make_shared<vector<shared_ptr<abstract_bipartite_edge>>>(size);
                for (uint32_t i = edge_vector->size() - size; i < edge_vector->size(); ++i) {
                    sub_edge_vector->at(i - edge_vector->size() + size) = edge_vector->at(i);
                }
                abstract_bipartite_graph_io::output_csv_file(sub_edge_vector, input_path + input_file_name + "-rear");
            }
        });
    }
    pool->barrier();
}

