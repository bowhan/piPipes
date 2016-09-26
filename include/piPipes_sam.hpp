#ifndef PIPIPES_SAM_HPP_
#define PIPIPES_SAM_HPP_

#include <iostream>
#include <string>
#include <vector>
#include <bitset>
#include <algorithm>
#include <iterator>
#include <unordered_map>
#include <boost/lexical_cast.hpp>
#include <boost/algorithm/string.hpp>


class Sam {
public:
    enum SAM_FLAG {
        PAIRED_END
        , EACH_END_ALIGNED
        , UNMAPPED
        , NEXT_UNMAPPED
        , REVERSE_COMPLEMENTED
        , NEXT_REVERSE_COMPLEMENTED
        , FIRST_SEG
        , SECOND_SEG
        , SECONDARY_ALIGNMENT
        , NOT_PASSING_QUALITY
        , PCR_DUP
        , FLAG_SIZE
    };

private:
    std::string qname_ = "";
    std::bitset<SAM_FLAG::FLAG_SIZE> flag_;
    std::string rname_ = "";
    uint64_t pos_ = 0;
    int mapq_ = 255;
    std::string cigar_ = "";
    std::string rnext_ = "";
    uint64_t pnext_ = 0;
    int64_t tlen_ = 0;
    std::string seq_ = "";
    std::string qual_ = "";
    std::unordered_map<std::string, std::string> optional_fields_;
public:
    /**! default constructor **/
    Sam() = default;

    /**! ctor from a line **/
    explicit Sam(const std::string& line);

    /**! copy ctor **/
    Sam(const Sam&) = default;

    /**! move ctor **/
    Sam(Sam&&) = default;

    /**! lvalue assignment operator **/
    Sam& operator=(const Sam& other) {
        if (this == &other) return *this;
        qname_ = other.qname_;
        flag_ = other.flag_;
        rname_ = other.rname_;
        pos_ = other.pos_;
        mapq_ = other.mapq_;
        cigar_ = other.cigar_;
        rnext_ = other.rnext_;
        pnext_ = other.pnext_;
        tlen_ = other.tlen_;
        seq_ = other.seq_;
        qual_ = other.qual_;
        optional_fields_ = other.optional_fields_;
        return *this;
    }

    /**! rvalue assignment operator **/
    Sam& operator=(Sam&& other) {
        assert (this != &other);
        qname_.swap(other.qname_);
        flag_ = other.flag_;
        rname_.swap(other.rname_);
        pos_ = other.pos_;
        mapq_ = other.mapq_;
        cigar_.swap(other.cigar_);
        rnext_.swap(other.rnext_);
        pnext_ = other.pnext_;
        tlen_ = other.tlen_;
        seq_.swap(other.seq_);
        qual_.swap(other.qual_);
        optional_fields_.swap(other.optional_fields_);
        other.~Sam();
        return *this;
    }

    /**! destructor **/
    ~Sam() = default;

    /** parse CIGAR **/
    std::vector<std::pair<char, uint32_t> > parseCIGAR() const;

    template <typename T>
    bool checkFlag(T f) const { return flag_.test(f); }

    template <typename T, typename... Args>
    bool checkFlag(T f, Args&& ... rest) const {
        if (sizeof ...(rest))
            return flag_.test(f) && checkFlag(std::forward<Args>(rest)...);
        return flag_.test(f);
    }

    std::string getTag(const std::string& key) const { return optional_fields_.at(key); }

    std::string mutation() const { return optional_fields_.at("MD"); }

    std::string getSeq() const { return seq_; }

    std::string getRevCompSeq() const;

    int getSize() const { return seq_.size(); }

    std::string getQueryName() const { return qname_; }

    friend std::ostream& operator<<(std::ostream& os, const Sam& sam);

};

#endif /* PIPIPES_SAM_HPP_ */
