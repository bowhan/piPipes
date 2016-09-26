#include <stdlib.h>
#include <string.h>
#include <iostream>
#include <fstream>
#include "piPipes_ping_pong.hpp"

boost::optional<int> multi_threading_queue::pop() {
    std::lock_guard<std::mutex> lk(mx_);
    boost::optional<int> ret{boost::none};
    if (!positions_.empty()) {
        ret = positions_.back();
        positions_.pop_back();
    }
    return ret;
}

void ParseBed(const std::string& file_name, pptype& ret) {
    using namespace std;
    ifstream in(file_name);
    if (!in) {
        cerr << "file " << file_name << " cannot be opened\n";
        exit(1);
    }
    char buffer[1024];
    char* saveptr;
    string chr;
    uint64_t start;
    uint64_t end;
    char strand;
    uint64_t nreads;
    uint64_t norigins;

    while (in.getline(buffer, 1024, '\n')) {
        chr = strtok_r(buffer, "\t", &saveptr);
        start = strtoul(strtok_r(NULL, "\t", &saveptr), NULL, 0);
        end = strtoul(strtok_r(NULL, "\t", &saveptr), NULL, 0);
        nreads = strtoul(strtok_r(NULL, "\t", &saveptr), NULL, 0);
        norigins = strtoul(strtok_r(NULL, "\t", &saveptr), NULL, 0);
        strand = strtok_r(NULL, "\t", &saveptr)[0];
        GenomicPosition g{chr
                          , strand == '+' ? (start) : (end - 1)
                          , strand
        };
        ret[g] += double(nreads) / norigins;
    }
}

void PingPong(const pptype& poola
              , const pptype& poolb
              , int distance
              , double* res
             ) {
    for (const auto& itera : poola) {
        const GenomicPosition& a = itera.first;
        if (a.pos() < distance) continue; // into the minus
        GenomicPosition b{a.chr()
                          , a.strand() == '+' ? (a.pos() + distance) : (a.pos() - distance)
                          , a.strand() == '+' ? '-' : '+'};
        auto hasb = poolb.find(b);
        if (hasb != poolb.end()) {
            *res += itera.second * hasb->second;
        }
    }
}

void Phasing(const pptype& poola
             , const pptype& poolb
             , int distance
             , double* res
            ) {
    for (const auto& itera : poola) {
        const GenomicPosition& a = itera.first;
        if (a.pos() < distance) continue; // into the minus
        GenomicPosition b{a.chr()
                          , a.strand() == '+' ? (a.pos() + distance) : (a.pos() - distance)
                          , a.strand() == '+' ? '+' : '-'};
        auto hasb = poolb.find(b);
        if (hasb != poolb.end()) {
            *res += std::min(itera.second, hasb->second);
        }
    }
}
