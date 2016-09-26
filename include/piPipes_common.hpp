#ifndef PIPIPES_PIPIPES_COMMON_HPP
#define PIPIPES_PIPIPES_COMMON_HPP

#include <boost/utility/string_ref.hpp>
#include <vector>
#include <sstream>

using StringView = boost::string_ref;

std::vector<StringView> Tokenize(StringView line, char delimiter);

template <class T>
bool StringViewTo(const StringView& s, T& t) {
    std::istringstream ss(s.to_string());
    ss >> t;
    return ss.fail() ? false : true;
}

#endif //PIPIPES_PIPIPES_COMMON_HPP
