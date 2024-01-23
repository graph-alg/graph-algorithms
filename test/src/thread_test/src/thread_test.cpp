
#include "thread_test/thread_test.h"

template<class container_type>
static auto split_task(const container_type& container,
                      const shared_ptr<scnu::thread_pool>& pool){
    uint32_t  thread_number = pool->get_thread_number();
    uint32_t task_count = container->size()/ thread_number;
    auto location_vector = make_shared<vector<decltype(container->begin())>>(thread_number + 1, container->begin());
    for(uint32_t i = thread_number - 1; i > 0;--i){
        pool->submit_task([=]{
            auto location = container->begin();
            std::advance(location, i * task_count);
            location_vector->at(i) = location;
        });
    }
    location_vector->at(thread_number) = container->end();
    pool->barrier();
    return location_vector;
}

static void prepare()
{
   auto num_set = make_shared<unordered_set<uint32_t>>();
   num_set->reserve(10);
   auto pool = make_shared<scnu::thread_pool>(5);
   for(uint32_t i = 0; i < 10; i++){
       pool->submit_task([=]{
          num_set->insert((i + 1)*(i + 1));
       });
   }
   pool->barrier();

   for(const auto&num:*num_set){
       printf("%d ", num);
   }
    printf("\n");
}

int main() {
    prepare();
}

