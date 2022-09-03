
#include "logger/simple_logger.h"


void free_test_function()
{
    scnu::simple_logger default_logger("test");
    LOG(default_logger,scnu::LOG_RANK::INFO) << "Log Test Begin!";
    LOG(default_logger,scnu::LOG_RANK::INFO) << "Log Test End!";
}

int main()
{
    free_test_function();
    return 0;
}




