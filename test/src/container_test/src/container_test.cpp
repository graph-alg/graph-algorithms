
#include "container_test/container_test.h"
#include "thread/thread_pool.h"
#include "time/simple_timer.h"

void init(const shared_ptr<unordered_map<uint32_t, shared_ptr<mutex>>>& num_map,
          const shared_ptr<unordered_map<uint32_t, uint32_t>>& num_map2,
          const shared_ptr<unordered_map<uint32_t, uint32_t>>& num_map3,
          uint32_t size){
    for(uint32_t i = 0; i < size; ++i){
        num_map->insert({i, make_shared<mutex>()});
        num_map2->insert({i, i});
        num_map3->insert({i, i});
    }
}

void init(const shared_ptr<unordered_map<uint32_t, shared_ptr<mutex>>>& num_map,
          const shared_ptr<unordered_map<uint32_t, uint32_t>>& num_map2,
          const shared_ptr<unordered_map<uint32_t, uint32_t>>& num_map3,
          uint32_t size,
          const shared_ptr<scnu::thread_pool>& pool){
    for(uint32_t i = 0; i < size; ++i){
        num_map->insert({i, shared_ptr<mutex>()});
    }

    auto thread_number = pool->get_thread_number();

    auto task_count = size / thread_number + 1;
    vector<uint32_t> location_vector{0 * task_count, 1 * task_count, 2 * task_count, 3 * task_count, 4 * task_count, 5 * task_count, size};
    for(uint32_t i = 0; i < thread_number; i++){
        pool->submit_task([=]{

            auto sub_begin = location_vector.at(i);
            auto sub_end = location_vector.at(i + 1);

            for(auto iter = sub_begin; iter!=sub_end; ++iter){
                num_map->at(iter) = make_shared<mutex>();
            }
        });
    }
    for(uint32_t i = 0; i < size; ++i){
        num_map2->insert({i, i});
        num_map3->insert({i, i});
    }
    pool->barrier();
}

void init2(const shared_ptr<unordered_map<uint32_t, shared_ptr<mutex>>>& num_map,
           const shared_ptr<unordered_map<uint32_t, uint32_t>>& num_map2,
           const shared_ptr<unordered_map<uint32_t, uint32_t>>& num_map3,
           uint32_t size,
           const shared_ptr<scnu::thread_pool>& pool){

    pool->submit_task([=]{
        for(uint32_t i = 0; i < size; ++i){
            num_map->insert({i, shared_ptr<mutex>()});
        }
    });

    pool->submit_task([=]{
        for(uint32_t j = 0; j < size; ++j){
            num_map2->insert({j, j});
        }
    });

    pool->submit_task([=]{
        for(uint32_t j = 0; j < size; ++j){
            num_map3->insert({j, j});
        }
    });
    pool->barrier();
}

int main()
{
    uint32_t size = 100000000;
    uint32_t thread_number = 6;
    {
        scnu::simple_timer t;
        auto num_map = make_shared<unordered_map<uint32_t, shared_ptr<mutex>>>();
        auto num_map2 = make_shared<unordered_map<uint32_t, uint32_t>>();
        auto num_map3= make_shared<unordered_map<uint32_t, uint32_t>>();
        init(num_map, num_map2, num_map3, size);
        printf("1: %lf\n", t.get_elapse_second());
    }

    {
        scnu::simple_timer t;
        auto num_map = make_shared<unordered_map<uint32_t, shared_ptr<mutex>>>();
        auto num_map2 = make_shared<unordered_map<uint32_t, uint32_t>>();
        auto num_map3= make_shared<unordered_map<uint32_t, uint32_t>>();
        auto pool= make_shared<scnu::thread_pool>(thread_number);
        init(num_map, num_map2, num_map3, size, pool);
        printf("2: %lf\n", t.get_elapse_second());
    }

    {
        scnu::simple_timer t;
        auto num_map = make_shared<unordered_map<uint32_t, shared_ptr<mutex>>>();
        auto num_map2 = make_shared<unordered_map<uint32_t, uint32_t>>();
        auto num_map3= make_shared<unordered_map<uint32_t, uint32_t>>();
        auto pool= make_shared<scnu::thread_pool>(thread_number);
        init2(num_map, num_map2, num_map3, size, pool);
        printf("3: %lf\n", t.get_elapse_second());
    }
}