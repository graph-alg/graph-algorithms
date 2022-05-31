#include "logger/simple_logger.h"

namespace scnu
{
    /**
     * @brief a construct function to open log file
     * @param log_file_name
     */
   simple_logger::simple_logger(const string &log_file_name)
   {
       log_rank_string_vector = {"INFO", "WARNING", "ERROR"};
       log_file.open(log_file_name, ios_base::app);
   }

    /**
     * @brief close the log file
     */
    simple_logger::~simple_logger()
    {
        log_file.close();
    }

    /**
      * @details  prepare output stream
      * @param other_log_rank
      * @param line
      * @param function
      * @return
      */
    simple_logger& simple_logger::start_log(LOG_RANK log_rank)
    {

        /**
         * @brief start a new stream
         */
        auto now = system_clock::now();
        auto current_time = system_clock::to_time_t(now);
        auto second_count = std::chrono::duration_cast<std::chrono::seconds>(now.time_since_epoch()).count()*1e6;
        auto microsecond_count = std::chrono::duration_cast<std::chrono::microseconds>(now.time_since_epoch()).count();
        auto fraction = microsecond_count-second_count;

        string_stream << put_time(localtime(&current_time), "%Y-%m-%d %X")<<"."<<fraction << " " <<
                      log_rank_string_vector.at(int(log_rank)) << " ";
        return *this;
    }

    /**
      * @details  prepare output stream
      * @param other_log_rank
      * @param line
      * @param function
      * @return
      */
    simple_logger& simple_logger::start_log(LOG_RANK log_rank,
                                            uint32_t line,
                                            const string &function)
    {

        /**
         * @brief start a new stream
         */
        auto now = system_clock::now();
        auto current_time = system_clock::to_time_t(now);
        auto second_count = std::chrono::duration_cast<std::chrono::seconds>(now.time_since_epoch()).count()*1e6;
        auto microsecond_count = std::chrono::duration_cast<std::chrono::microseconds>(now.time_since_epoch()).count();
        auto fraction = microsecond_count-second_count;

        string_stream << put_time(localtime(&current_time), "%Y-%m-%d %X")<<"."<<fraction << " " <<
                      log_rank_string_vector.at(int(log_rank)) << " " << line << " " << function << " ";
        return *this;
    }
}



