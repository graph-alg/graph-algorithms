#include "unique_temporal_graph_process_test/unique_temporal_graph_process_test.h"

int main(int argc, char **argv){
    if(argc < 3){
        std::cout<<"Usage: input path, output path, and thread number!";
    }
    std::string input_path(argv[1]);
    std::string output_path(argv[2]);

    auto thread_number = std::stoul(argv[3]);
    temporal_graph_io::unique_output_csv_file(input_path, output_path,thread_number);
    return 0;
}