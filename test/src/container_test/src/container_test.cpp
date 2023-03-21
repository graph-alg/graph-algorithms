
#include "container_test/container_test.h"


int main()
{

    auto v = make_shared<vector<uint32_t>>(5, 1);

    auto copy_v = make_shared<vector<uint32_t>>(*v);

    v->at(0) = 2;

    for(const auto &n:*v){
        printf("%d ",n);
    }
    printf("\n");
}