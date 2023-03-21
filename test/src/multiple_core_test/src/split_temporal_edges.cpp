
#include "multiple_core_test/multiple_core_test.h"

int main(int argc, char **argv) {
    if(argc < 4){
        std::cout<<"Usage: input path, split size,  and thread number!";
    }

    const string input_path = argv[1];
    const double split_size = std::stod(argv[2]);
    const uint32_t thread_number = std::stoul(argv[3]);


    auto rd = make_shared<random_device>();
    auto directory = path(input_path);
    auto pool = make_shared<thread_pool>(thread_number);
    for (auto &file_iter: std::filesystem::directory_iterator(input_path)) {
        if (!std::filesystem::is_regular_file(file_iter)) {
            continue;
        }
        pool->submit_task([=, &split_size] {

            const auto &file = file_iter.path();
            auto input_file_name = file.filename().string();
            auto edge_vector = temporal_graph_io::get_edge_vector(input_path, input_file_name);

            auto length = uint32_t (std::ceil(edge_vector->size() * split_size));
            auto part_end = edge_vector->begin();
            std::advance(part_end, length);

            auto sub_edge_vector = make_shared<vector<shared_ptr<temporal_edge>>>();
            copy(edge_vector->begin(), part_end, inserter(*sub_edge_vector, sub_edge_vector->end()));

            temporal_graph_io::output_csv_file(sub_edge_vector, input_path + input_file_name + "-part");
        });
    }
    pool->barrier();
}

