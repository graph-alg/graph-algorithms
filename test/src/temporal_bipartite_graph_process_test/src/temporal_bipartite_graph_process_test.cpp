
#include "temporal_bipartite_graph_process_test/temporal_bipartite_graph_process_test.h"

int main(int argc, char **argv){
    if(argc < 3){
        std::cout<<"Usage: input path, output path, and thread number!";
    }
    std::string input_path(argv[1]);
    std::string output_path(argv[2]);
    auto thread_number = std::stoul(argv[3]);
    temporal_bipartite_graph_io::store_graph(input_path, output_path, thread_number);
    return 0;
}
