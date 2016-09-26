#ifndef HAIRPIN_HPP_
#define HAIRPIN_HPP_

#include <string>
#include <iostream>
#include <vector>
#include <memory>
#include <algorithm>

enum hairpinPosition {
    unknown
    , fivePrimeArm
    , loop
    , threePrimeArm
};

/**
 * miRBed is a class that record each line of mature miRNA mapping to hairpin result in bed format
 */
class miRBed {
private:
    /**! index name being mapped to, which is the hairpin name **/
    std::string hairpin_name_ = "";
    /**! starting position on the hairpin been mapped to, zero based **/
    int start_ = -1;
    /**! ending position on the hairpin been mapped to, open **/
    int end_ = -1;
    /**! name of the mature miRNA **/
    //	std::string matureName = "";
    /**! 4th col of bed2 format is the number of reads **/
    int reads_ = 0;
    /// 5th col of this bed2 is ignored
    /// 6th col of this bed2 has to be '+'
    /**! 7th col of bed2 format is the sequence **/
    std::string sequence = "";
    /**! annotation of this mature miRNA being 5'arm/3'arm/loop **/
    hairpinPosition hp = hairpinPosition::unknown;

public:
    /**! explicit ctor for constructing this bed2 entry from one line (string) **/
    explicit miRBed(const std::string& line);

    /**! let class Hairpin to have access to miRBed **/
    friend class Hairpin;

    /**! outputing method, for debug purpose **/
    friend std::ostream& operator<<(std::ostream& os, miRBed& mb) {
        os << mb.hairpin_name_ << '\t'
           << mb.start_ << '\t'
           << mb.end_ << '\t'
           << mb.reads_ << '*' << '\t'
           << '+' << '\t'
           << mb.sequence << '\n';
        return os;
    }
};

/**
 * Hairpin is the class that store everything about one hairpin
 * it has an internal vector to store the raw bed2 format (miRBed)
 */
class Hairpin {
protected:
    /**! name of the hairpin **/
    std::string name_ = "";
    /**! the annotated starting point of the 5' mature miRNA **/
    int five_start_ = -1;
    /**! the annotated ending point of the 5' mature  miRNA **/
    int five_end_ = -1;
    /**! the annotated starting point of the 3' mature miRNA **/
    int three_start_ = -1;
    /**! the annoated ending point of the 3' mature miRNA **/
    int three_end_ = -1;
    /**! the internal vector to store every line of bed2 mapping to this hairpin **/
    std::vector<miRBed> bed_data_;
    /**! 5' heterogeneity of 5' arm miRNA **/
    double five_heter_five_arm_ = 0.0;
    /**! 3' heterogeneity of 5' arm miRNA **/
    double three_heter_five_arm_ = 0.0;
    /**! 5' heterogeneity of 3' arm miRNA **/
    double five_heter_three_arm_ = 0.0;
    /**! 3' heterogeneity of 3' arm miRNA **/
    double three_heter_three_arm_ = 0.0;
    /**! counter to count how many 5' miRNA mapped in this bed file **/
    int five_arm_counter_ = 0;
    /**! counter to count how many 3' miRNA mapped in this bed file **/
    int three_arm_counter_ = 0;
public:
    /**! constructing Hairpin from one line in bed2 format
    dme-mir-1	17	39	dme-miR-1-5p	255	+
    this is supposed to be used to annotate hairpin by reading mapping result of mature.fa to hairpin.fa
     **/
    explicit Hairpin(const std::string& line);

    /**! when reading this hairpin the second time (if the first time the mature was 5p and this time it is 3p
     *  cannot use ctor so use update
     */
    void update(const Hairpin& hp);

    /**! when reading real mapping data, we store that bed into the internal vector
     * using this function
     */
    void insert(const std::string& line);

    /**! get name of the hairpin **/
    std::string getName() const {
        return name_;
    }

    /**! get total size of bed2 in the internal vector **/
    unsigned long size() const {
        return bed_data_.size();
    }

    /**! output all the bed2s in the internal vector, adjusting their 5' end and 3' end by to the relative distance to the annoated mature **/
    std::ostream& outputAll(std::ostream& os) const;

    /**! calculate the 5' and 3' heterogeneity **/
    void calculateHeterogeneity();

    double getFivePrimeHeterogeneityOfFivePrimeArm() const {
        return five_heter_five_arm_;
    }

    double getThreePrimeHeterogeneityOfFivePrimArm() const {
        return three_heter_five_arm_;
    }

    double getFivePrimeHeterogeneityOfThreePrimeArm() const {
        return five_heter_three_arm_;
    }

    double getThreePrimeHeterogeneityOfThreePrimArm() const {
        return three_heter_three_arm_;
    }

    int getFivePrimeCounter() const {
        return five_arm_counter_;
    }

    int getThreePrimeCounter() const {
        return three_arm_counter_;
    }
};

#endif /* HAIRPIN_HPP_ */
