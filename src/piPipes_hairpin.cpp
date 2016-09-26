#include "piPipes_hairpin.hpp"

miRBed::miRBed(const std::string& line) {
    std::string::const_iterator iter1{line.cbegin()}, iter2{line.cbegin()};
    while (*++iter2 != '\t');
    hairpin_name_ = std::string(iter1, iter2);
    iter1 = ++iter2;
    while (*++iter2 != '\t');
    start_ = stoi(std::string(iter1, iter2));
    iter1 = ++iter2;
    while (*++iter2 != '\t');
    end_ = stoi(std::string(iter1, iter2));
    iter1 = ++iter2;
    while (*++iter2 != '\t');
    reads_ = stoi(std::string(iter1, iter2));
    while (*++iter2 != '\t');
    while (*++iter2 != '\t');
    ++iter2;
    sequence = std::string(iter2, line.cend());

}


Hairpin::Hairpin(const std::string& line) {
    int start{0}, end{0};
    std::string matureName;
    std::string::const_iterator iter1{line.cbegin()}, iter2{line.cbegin()};
    while (*++iter2 != '\t');
    name_ = std::string(iter1, iter2);
    iter1 = ++iter2;
    while (*++iter2 != '\t');
    start = stoi(std::string(iter1, iter2));
    iter1 = ++iter2;
    while (*++iter2 != '\t');
    end = stoi(std::string(iter1, iter2));
    iter1 = ++iter2;
    while (*++iter2 != '\t');
    matureName = std::string(iter1, iter2);
    /**! this is important, in the ideal situation, every mature miRNA should have a 5p or 3p in their name
     * according to the newest version of miRBase, fly miRNA all have this done, but mouse doesn't
     * so the pipeline privide a small script fix_mature_missing_5p3p to fix this.
     * So ideally this problem should be solved
     */
    /// find 5p miRNA
    if (matureName.find("5p") != std::string::npos) {
        /// if fivePrimeStart has been settled, which means there are two line has 5p of this hairpin
        if (five_start_ != -1) {
            std::cerr << "Warning: fivePrime has already been set for " << name_ << " at line " << line
                      << std::endl;
        }
        /// annotate the 5' and 3' end of this hairpin
        five_start_ = start;
        five_end_ = end;
    } else if (matureName.find("3p") != std::string::npos) {
        if (three_start_ != -1) {
            std::cerr << "Warning: threePrime has already been set for " << name_ << " at line " << line
                      << std::endl;
        }
        three_start_ = start;
        three_end_ = end;
    }
        /// if there isn't any 5p or 3p
    else {
        std::cerr << "Warning: the mature miRNA has no arm information in its name " << line << std::endl;
    }
}

void Hairpin::update(const Hairpin& hp) {
    /// if 5' have not been annoated
    if (hp.five_start_ != -1 && five_start_ == -1) {
        five_start_ = hp.five_start_;
        five_end_ = hp.five_end_;
    }
        /// if 3' arm miRNA has not been annotated
    else if (hp.three_start_ != -1 && three_start_ == -1) {
        three_start_ = hp.three_start_;
        three_end_ = hp.three_end_;
    }
}

void Hairpin::insert(const std::string& line) {
    bed_data_.emplace_back(line);
    miRBed& b = bed_data_.back();
    /**! make the judgement whether this read is 5' arm or 3' arm or the loop **/
    auto overlapBaseWith5arm = std::min(b.end_, five_end_) - std::max(b.start_, five_start_);
    if (overlapBaseWith5arm > (five_end_ - five_start_) / 2) {
        b.hp = hairpinPosition::fivePrimeArm;
        return;
    }
    auto overlapBaseWith3arm = std::min(b.end_, three_end_) - std::max(b.start_, three_start_);
    if (overlapBaseWith3arm > (three_end_ - three_start_) / 2) {
        b.hp = hairpinPosition::threePrimeArm;
        return;
    }
    auto overlapBaseWithLoop = std::min(b.end_, three_start_) - std::max(b.start_, five_end_);
    if (overlapBaseWithLoop > (three_start_ - five_end_) / 2) {
        b.hp = hairpinPosition::loop;
        return;
    }
    b.hp = hairpinPosition::unknown;
    return;
}

std::ostream& Hairpin::outputAll(std::ostream& os) const {
    for (const auto& bed : bed_data_) {
        os << bed.hairpin_name_ << '\t';
        if (bed.hp == hairpinPosition::fivePrimeArm) {
            os << bed.start_ - five_start_ << '\t';
            os << bed.end_ - five_end_ << '\t';
        } else if (bed.hp == hairpinPosition::threePrimeArm) {
            os << bed.start_ - three_start_ << '\t';
            os << bed.end_ - three_end_ << '\t';
        } else if (bed.hp == hairpinPosition::loop) {
            os << bed.start_ - five_end_ << '\t';
            os << bed.end_ - three_start_ << '\t';
        } else {
            os << '*' << '\t';
            os << '*' << '\t';
        }
        os << bed.reads_ << '\t';
        os << bed.hp << "\t+\t";
        os << bed.sequence << '\n';
    }
    return os;;
}

void Hairpin::calculateHeterogeneity() {
    for (const auto& mir : bed_data_) {
        if (mir.hp == hairpinPosition::fivePrimeArm) {
            five_heter_five_arm_ += std::abs(mir.start_ - five_start_) * mir.reads_;
            three_heter_five_arm_ += std::abs(mir.end_ - five_end_) * mir.reads_;
            five_arm_counter_ += mir.reads_;
        } else if (mir.hp == hairpinPosition::threePrimeArm) {
            five_heter_three_arm_ += std::abs(mir.start_ - three_start_) * mir.reads_;
            three_heter_three_arm_ += std::abs(mir.end_ - three_end_) * mir.reads_;
            three_arm_counter_ += mir.reads_;
        }
    }
    if (five_arm_counter_) {
        five_heter_five_arm_ /= five_arm_counter_;
        three_heter_five_arm_ /= five_arm_counter_;
    }
    if (three_arm_counter_) {
        five_heter_three_arm_ /= three_arm_counter_;
        three_heter_three_arm_ /= three_arm_counter_;
    }
}
