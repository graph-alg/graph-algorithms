#include "temporal_graph_process_test/temporal_graph_process_test.h"

int main(int argc, char **argv){
    if(argc < 3){
        std::cout<<"Usage: input path, output path, and thread number!";
    }
    std::string input_path(argv[1]);
    std::string output_path(argv[2]);

    auto thread_number = std::stoul(argv[3]);
    auto pool = make_shared<scnu::thread_pool>(thread_number);
    temporal_graph_io::output_unique_graph(input_path, output_path, pool);
    return 0;
}