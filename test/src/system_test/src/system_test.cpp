
#include "system_test/system_test.h"

int main()
{
    auto memory_size = process_information::get_memory(getpid());
    cout<<double (memory_size)/1024<<"MB\n";

    auto num_vector = make_shared<vector<uint32_t>>(102400);
    for(uint32_t i = 0;i<102400;++i){
        num_vector->at(i) = i+1;
    }

    cout<< (sizeof(*num_vector)+sizeof(uint32_t)*num_vector->size())/1024<<"KB\n";

    memory_size = process_information::get_memory(getpid());
    cout<<double (memory_size)/1024<<"MB\n";
}

