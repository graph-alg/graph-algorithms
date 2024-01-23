/******************************************************************************************************************
* @copyright  : Boost Software License
* @file       : extend_string.cpp
* @brief      : A subclass of string for implementing extra functions
* @version    : 1.1
* @date       : 2019/10/09
******************************************************************************************************************/

#include "string/extend_string.h"

namespace scnu
{
    /**
     * @brief construct a extend_string with a given std::string
     * @param otherTextString
     */
    extend_string::extend_string(std::string  other_string):
        data_string(std::move(other_string))
    {

    }

    /**
     * @brief get the delimiter of this string
     * @return
     */
    std::string& extend_string::get_delimiter() {
        return delimiter;
    }

    /**
     * @brief get std::string of this extend_string
     * @return
     */
    std::string& extend_string::get_string()
    {
        return data_string;
    }

    /**
     * @brief split this string by the given other_delimiter
     * @param other_delimiter
     * @return
     */
    std::shared_ptr<std::vector<std::string>> extend_string::split(const std::string& other_delimiter)
    {
        std::istringstream input_stream(data_string);
        delimiter = other_delimiter;
        return std::make_shared<std::vector<std::string>>(std::istream_iterator<extend_string>(input_stream),
                                              std::istream_iterator<extend_string>());
    }

    /**
     * @brief overload >> for input
     * @param input_stream
     * @param other_string
     * @return
     */
    std::istream& operator>>(std::istream& input_stream, extend_string& other_string)
    {
        return std::getline(input_stream, other_string, *(extend_string::get_delimiter().c_str()));
    }

    /**
     * @brief overload << for output
     * @param output_stream
     * @param other_string
     * @return
     */
    std::ostream& operator<<(std::ostream & output_stream, extend_string& other_string)
    {
        return output_stream<<other_string.get_string();
    }
}
