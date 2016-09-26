#include "piPipes_common.hpp"

std::vector<StringView> Tokenize(StringView line, char delimiter) {
    std::vector<StringView> strrefs;
    auto iter1 = line.begin();
    auto iter2 = iter1;
    size_t l = 0;
    for (;;) {
        while (*iter1 != delimiter && iter1 != line.end()) {
            ++iter1;
            ++l;
        }
        if (iter1 == line.end()) {
            break;
        }
        if (l) {
            strrefs.emplace_back(iter2, l);
        }
        iter2 = ++iter1;
        l = 0;
    }
    strrefs.emplace_back(iter2, l);
    return strrefs;
}


template < >
bool StringViewTo(const StringView& s, int& ret) {
    if (s.empty()) return false;
    auto iter = s.cbegin();
    bool negative = false;
    if (*iter == '-') {
        negative = true;
        ++iter;
        if (iter == s.cend()) return false;
    }
    ret = 0;
    for (; iter != s.cend(); ++iter) {
        if (*iter < '0' || *iter > '9') return false;
        ret *= 10;
        ret += *iter - '0';
    }
    if (negative) ret = -ret;
    return true;
}
