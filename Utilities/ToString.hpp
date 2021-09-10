#ifndef TO_STRING_HPP
#define TO_STRING_HPP

#include <string>

template<class T>
std::string inline read_from_file(const T& object)
{
    return object.to_string();
}

#endif // TO_STRING_HPP