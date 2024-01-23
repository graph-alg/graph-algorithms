
#include "bipartite_core_test/bipartite_core_test.h"

int main(int argc, char **argv) {
    if(argc < 2){
        std::cout<<"Usage: input path, and thread number!";
    }

    string input_path = argv[1];
    uint32_t  thread_number = std::stoul(argv[2]);

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

            shuffle(edge_vector->begin(), edge_vector->end(), *random_generator::get_default_engine(rd));
            abstract_bipartite_graph_io::output_csv_file(edge_vector, input_path + input_file_name);
        });
    }
    pool->barrier();
}

