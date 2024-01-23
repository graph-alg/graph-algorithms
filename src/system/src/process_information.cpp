
#include "system/process_information.h"

namespace scnu
{
    uint32_t process_information::get_memory(pid_t p) {
        string input_file_name="/proc/"+std::to_string(p)+"/status";
        ifstream input_file(input_file_name);
        string line;
        uint32_t line_count = 0;
        while (getline(input_file,line).good())
        {
            ++line_count;
            if(line.empty()){
                break;
            }
            if(line_count == VMRSS_LINE){
                auto line_vector = string_algorithm::regex_split(line,":");
                auto memory_size = std::stoul(line_vector->at(1));
                return memory_size;
            }
        }
        return 0;
    }
}