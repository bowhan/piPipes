#ifndef PIPIPES_PIPIPES_PING_PONG_HPP
#define PIPIPES_PIPIPES_PING_PONG_HPP

#include <deque>
#include <thread>
#include <mutex>
#include <unordered_map>
#include <boost/optional.hpp>

class multi_threading_queue {
    std::deque<int> positions_;
    std::mutex mx_;
public:
    explicit multi_threading_queue(int N)
        : positions_{} {
        for (int i{}; i < N; ++i) {
            positions_.push_back(i);
        }
    }

    boost::optional<int> pop();
};

template <typename C>
class PingPongPlayer {
    using queue_type = multi_threading_queue;
    using data_type = C;
    using function_type = void (*)(const data_type&, const data_type&, int, double* const);
private:
    const data_type& a_;
    const data_type& b_;
    queue_type& tasks_;
    double* answers_;
    function_type func_;

public:
    PingPongPlayer<C>(const data_type& a, const data_type& b, queue_type& tasks, double* answers, function_type f)
        :
        a_{a}, b_{b}, tasks_{tasks}, answers_{answers}, func_{f} {}

    void operator()() {
        while (1) {
            auto n = tasks_.pop();
            if (n) {
                (*func_)(a_, b_, *n, answers_ + *n);
            } else {
                break;
            }
        }
    }
};


class GenomicPosition {
private:
    std::string chr_;
    uint64_t pos_;
    char strand_;
public:

    GenomicPosition(const std::string& c, uint64_t p, char s)
        :
        chr_(c), pos_(p), strand_(s) {}

    const std::string& chr() const { return chr_; }

    std::string& chr() { return chr_; }

    const uint64_t& pos() const { return pos_; };

    uint64_t& pos() { return pos_; };

    const char& strand() const { return strand_; }

    char& strand() { return strand_; }


    bool operator==(const GenomicPosition& other) const {
        return chr_ == other.chr_
            && pos_ == other.pos_
            && strand_ == other.strand_;
    }

};

namespace std {
template < >
struct hash<GenomicPosition> {
    std::size_t operator()(const GenomicPosition& g) const {
        return
            hash<std::string>()(g.chr())
                ^ hash<uint64_t>()(g.pos())
                ^ hash<char>()(g.strand());
    }
}; // hash
} // namespace std

using pptype = std::unordered_map<GenomicPosition, double>;

using MyPingPongPlayer = PingPongPlayer<pptype>;

void ParseBed(const std::string& filename, pptype&);

void PingPong(const pptype&, const pptype&, int, double*);

void Phasing(const pptype&, const pptype&, int, double*);

#endif //PIPIPES_PIPIPES_PING_PONG_HPP
