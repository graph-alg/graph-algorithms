
#include "graph_process_test/graph_process_test.h"

int main(int argc, char **argv) {
    if(argc < 3){
        std::cout<<"Usage: please set input and output paths of the dataset.";
    }
    std::string input_path(argv[1]);
    std::string output_path(argv[2]);

    auto thread_number = std::stoul(argv[3]);
    abstract_graph_io::store_graph(input_path, output_path, thread_number);
    return 0;
}






