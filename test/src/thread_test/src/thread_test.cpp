
#include "thread_test/thread_test.h"

void prepare()
{
   auto num_vector = make_shared<vector<uint32_t>>(10, 0);
   auto pool = make_shared<scnu::thread_pool>(4);
   for(uint32_t i = 0; i < num_vector->size(); ++i){
       pool->submit_task([=]{
           num_vector->at(i) = (i + 1) * (i + 1);
       });
   }
   pool->barrier();

   for(const auto &v:*num_vector){
       printf("%u ", v);
   }
    printf("\n");
}

int main() {
    prepare();
}

