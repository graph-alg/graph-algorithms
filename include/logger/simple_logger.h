/******************************************************************************************************************
* @copyright  : Boost Software License
* @file       : simple_logger.h
* @brief      : A simple_logger for output information
* @version    : 1.1
* @date       : 2020/10/11
******************************************************************************************************************/
#pragma once
#include "logger/logger_utility.h"

namespace scnu {
    /**
     * @brief the status of simple_logger
     */
    enum LOG_RANK {
        INFO,
        WARNING,
        ERROR
    };

    /**
     * @class simple_logger
     * @details a simple_logger for output all information
     */
    class simple_logger {
    public:
        explicit simple_logger(const string &log_file_name);

        ~simple_logger();

        simple_logger& start_log(LOG_RANK log_rank);

        simple_logger& start_log(LOG_RANK log_rank,
                                      uint32_t line,
                                      const string &function);

        /**
         * @details overload << operator
         * @remarks if message is newline character, then output log
         * @tparam T
         * @param message
         * @return
         */
        template<typename T>
        simple_logger& operator<<(const T& message)
        {
            string_stream<<message;
            auto result_string = string_stream.str();
            if(*result_string.rbegin()=='\n')
            {
                log_file<<result_string;
                std::cerr<<result_string;
                string_stream.str("");
            }
            return *this;
        }

    private:
        vector<string> log_rank_string_vector;

        ofstream log_file;

        stringstream string_stream;
    };

    /**
     * @macro macros for simple_logger recording
     * @brief a public interface to write log
     * @param log
     * @param log_rank
     */
    #define  LOG(log,log_rank) log.start_log(log_rank)
}




