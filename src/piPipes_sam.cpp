#include "piPipes_sam.hpp"

Sam::Sam(const std::string& line) {
    std::string::const_iterator iter1{line.cbegin()}, iter2{line.cbegin()};
    while (*++iter2 != '\t');
    qname_ = std::string {iter1, iter2};
    iter1 = ++iter2;

    while (*++iter2 != '\t');
    flag_ = std::bitset<SAM_FLAG::FLAG_SIZE> {std::stoul(std::string {iter1, iter2})};
    iter1 = ++iter2;

    while (*++iter2 != '\t');
    rname_ = std::string {iter1, iter2};
    iter1 = ++iter2;

    while (*++iter2 != '\t');
    pos_ = std::stoul(std::string {iter1, iter2});
    iter1 = ++iter2;

    while (*++iter2 != '\t');
    mapq_ = std::stoi(std::string {iter1, iter2});
    iter1 = ++iter2;

    while (*++iter2 != '\t');
    cigar_ = std::string {iter1, iter2};
    iter1 = ++iter2;

    while (*++iter2 != '\t');
    rnext_ = std::string {iter1, iter2};
    iter1 = ++iter2;

    while (*++iter2 != '\t');
    pnext_ = std::stoul(std::string {iter1, iter2});
    iter1 = ++iter2;

    while (*++iter2 != '\t');
    tlen_ = std::stol(std::string {iter1, iter2});
    iter1 = ++iter2;

    while (*++iter2 != '\t');
    seq_ = std::string {iter1, iter2};
    iter1 = ++iter2;

    while (*++iter2 != '\t');
    qual_ = std::string {iter1, iter2};
    std::string Tag, Value;
    while (iter2 != line.cend()) {
        iter1 = ++iter2;
        iter2 += 2;
        Tag = std::string {iter1, iter2};
        iter1 += 5;
        while ((*iter2 != '\t') && (iter2 != line.cend()))
            ++iter2;
        Value = std::string {iter1, iter2};
        optional_fields_.insert(std::make_pair(Tag, Value));
    }
}


std::vector<std::pair<char, uint32_t> > Sam::parseCIGAR() const {
    std::vector<std::pair<char, uint32_t> > CIGARs;
    uint32_t t = 0;
    for (char x : cigar_) {
        if (x >= '0' && x <= '9') {
            t = 10 * t + x - '0';
        } else {
            CIGARs.emplace_back(x, t);
            t = 0;
        }
    }
    return CIGARs;
}


std::string Sam::getRevCompSeq() const {
    std::string rcSeq{seq_.rbegin(), seq_.rend()};
    for (auto& c : rcSeq) {
        switch (c) {
            case 'A': c = 'T';
                break;
            case 'T': c = 'A';
                break;
            case 'C': c = 'G';
                break;
            case 'G': c = 'C';
                break;
            default : throw "illegal char";
        }
    } /* end of RC */
    return rcSeq;
}

std::ostream& operator<<(std::ostream& os, const Sam& sam) {
    os << sam.qname_ << '\t'
       << sam.flag_.to_string() << '\t'
       << sam.rname_ << '\t'
       << sam.pos_ << '\t'
       << sam.mapq_ << '\t'
       << sam.cigar_ << '\t'
       << sam.rnext_ << '\t'
       << sam.pnext_ << '\t'
       << sam.tlen_ << '\t'
       << sam.seq_ << '\t'
       << sam.qual_;
    for (const auto& p : sam.optional_fields_)
        os << '\t' << p.first << ':' << p.second;
    os << '\n';
    return os;
}


