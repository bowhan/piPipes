/*
 * This script is part of the NGS pipeline developed in the Zamore lab and the Weng lab
 * Bo Han (bowhan@me.com)
 */

#ifndef HAIRPIN_HPP_
#define HAIRPIN_HPP_

#include <string>
#include <vector>
#include <ostream>
#include <algorithm>
#include <memory>
#include <boost/iterator/iterator_facade.hpp>

enum hairpinPosition {
	unknown,
	fivePrimeArm,
	loop,
	threePrimeArm
};
/**
 * miRBed is a class that record each line of mature miRNA mapping to hairpin result in bed format
 */
class miRBed
{
private:
	/**! index name being mapped to, which is the hairpin name **/
	std::string hairpinName = "";
	/**! starting position on the hairpin been mapped to, zero based **/
	int start = -1;
	/**! ending position on the hairpin been mapped to, open **/
	int end = -1;
	/**! name of the mature miRNA **/
//	std::string matureName = "";
	/**! 4th col of bed2 format is the number of reads **/
	int numOfSeq = 0;
	/// 5th col of this bed2 is ignored
	/// 6th col of this bed2 has to be '+'
	/**! 7th col of bed2 format is the sequence **/
	std::string sequence = "";
	/**! annotation of this mature miRNA being 5'arm/3'arm/loop **/
	hairpinPosition hp = hairpinPosition::unknown;

public:
	/**! explicit ctor for constructing this bed2 entry from one line (string) **/
	explicit miRBed (const std::string& line) {
		std::string::const_iterator iter1 {line.cbegin()}, iter2{line.cbegin()};
		while (*++iter2!='\t');
		hairpinName = std::string (iter1, iter2);
		iter1=++iter2;
		while (*++iter2!='\t');
		start = stoi (std::string (iter1, iter2));
		iter1=++iter2;
		while (*++iter2!='\t');
		end = stoi (std::string (iter1, iter2));
		iter1=++iter2;
		while (*++iter2!='\t');
		numOfSeq = stoi (std::string (iter1, iter2));
		while (*++iter2!='\t');
		while (*++iter2!='\t');
		++iter2;
		sequence = std::string (iter2, line.cend ());

	}
	/**! let class Hairpin to have access to miRBed **/
	friend class Hairpin;
	/**! outputing method, for debug purpose **/
	friend std::ostream & operator << (std::ostream& os, miRBed& mb) {
		os << mb.hairpinName << '\t' << mb.start << '\t' << mb.end << '\t' << mb.numOfSeq << '*' << '\t' << '+' << '\t' << mb.sequence << '\n';
		return os;
	}
};
/**
 * Hairpin is the class that store everything about one hairpin
 * it has an internal vector to store the raw bed2 format (miRBed)
 */
class Hairpin
{
private:
	/**! Hairpin has a vector inside, so this API also privide an iterator **/
	friend class HairpinIterator;
protected:
	/**! name of the hairpin **/
	std::string name = "";
	/**! the annotated starting point of the 5' mature miRNA **/
	int fivePrimeStart = -1;
	/**! the annotated ending point of the 5' mature  miRNA **/
	int fivePrimeEnd = -1;
	/**! the annotated starting point of the 3' mature miRNA **/
	int threePrimeStart = -1;
	/**! the annoated ending point of the 3' mature miRNA **/
	int threePrimeEnd = -1;
	/**! the internal vector to store every line of bed2 mapping to this hairpin **/
	std::vector<miRBed> beds;
	/**! 5' heterogeneity of 5' arm miRNA **/
	double fivePrimeHeterogeneityOfFivePrimeArm = 0.0;
	/**! 3' heterogeneity of 5' arm miRNA **/
	double threePrimeHeterogeneityOfFivePrimArm = 0.0;
	/**! 5' heterogeneity of 3' arm miRNA **/
	double fivePrimeHeterogeneityOfThreePrimeArm = 0.0;
	/**! 3' heterogeneity of 3' arm miRNA **/
	double threePrimeHeterogeneityOfThreePrimArm = 0.0;
	/**! counter to count how many 5' miRNA mapped in this bed file **/
	int fiveArmCounter = 0;
	/**! counter to count how many 3' miRNA mapped in this bed file **/
	int threeArmCounter = 0;
public:
	/**! constructing Hairpin from one line in bed2 format
	dme-mir-1	17	39	dme-miR-1-5p	255	+
	this is supposed to be used to annotate hairpin by reading mapping result of mature.fa to hairpin.fa
	 **/
	explicit Hairpin (const std::string& line) {
		int start {0}, end {0};
		std::string matureName;
		std::string::const_iterator iter1 {line.cbegin()}, iter2{line.cbegin()};
		while (*++iter2!='\t');
		name = std::string (iter1, iter2);
		iter1=++iter2;
		while (*++iter2!='\t');
		start = stoi (std::string (iter1, iter2));
		iter1=++iter2;
		while (*++iter2!='\t');
		end = stoi (std::string (iter1, iter2));
		iter1=++iter2;
		while (*++iter2!='\t');
		matureName = std::string (iter1, iter2);
		/**! this is important, in the ideal situation, every mature miRNA should have a 5p or 3p in their name
		 * according to the newest version of miRBase, fly miRNA all have this done, but mouse doesn't
		 * so the pipeline privide a small script fix_mature_missing_5p3p to fix this.
		 * So ideally this problem should be solved
		 */
		/// find 5p miRNA
		if (matureName.find ("5p")!=std::string::npos) {
		/// if fivePrimeStart has been settled, which means there are two line has 5p of this hairpin
			if (fivePrimeStart != -1) {
				std::cerr << "Warning: fivePrime has already been set for " << name << " at line " << line << std::endl;
			}
			/// annotate the 5' and 3' end of this hairpin
			fivePrimeStart = start;
			fivePrimeEnd = end;
		}
		else if (matureName.find ("3p")!=std::string::npos) {
			if (threePrimeStart != -1) {
				std::cerr << "Warning: threePrime has already been set for " << name << " at line " << line << std::endl;
			}
			threePrimeStart = start;
			threePrimeEnd = end;
		}
		/// if there isn't any 5p or 3p
		else {
			std::cerr << "Warning: the mature miRNA has no arm information in its name " << line << std::endl;
		}
	}
	/**! when reading this hairpin the second time (if the first time the mature was 5p and this time it is 3p
	 *  cannot use ctor so use update
	 */
	void update (const Hairpin& hp) {
		/// if 5' have not been annoated
		if (hp.fivePrimeStart!=-1 && fivePrimeStart==-1) {
			fivePrimeStart = hp.fivePrimeStart;
			fivePrimeEnd = hp.fivePrimeEnd;
		}
		/// if 3' arm miRNA has not been annotated
		else if (hp.threePrimeStart!=-1 && threePrimeStart==-1) {
			threePrimeStart = hp.threePrimeStart;
			threePrimeEnd = hp.threePrimeEnd;
		}
	}
	/**! when reading real mapping data, we store that bed into the internal vector
	 * using this function
	 */
	void insert (const std::string & line) {
		beds.emplace_back (line);
		miRBed& b = beds.back ();
		/**! make the judgement whether this read is 5' arm or 3' arm or the loop **/
		auto overlapBaseWith5arm = std::min (b.end, fivePrimeEnd) - std::max (b.start, fivePrimeStart);
		if ( overlapBaseWith5arm > (fivePrimeEnd - fivePrimeStart)/2 ) {
			b.hp = hairpinPosition::fivePrimeArm;
			return;
		}
		auto overlapBaseWith3arm = std::min (b.end, threePrimeEnd) - std::max (b.start, threePrimeStart);
		if ( overlapBaseWith3arm > (threePrimeEnd - threePrimeStart)/2) {
			b.hp = hairpinPosition::threePrimeArm;
			return;
		}
		auto overlapBaseWithLoop = std::min (b.end, threePrimeStart) - std::max (b.start, fivePrimeEnd);
		if ( overlapBaseWithLoop > (threePrimeStart - fivePrimeEnd)/2) {
			b.hp = hairpinPosition::loop;
			return;
		}
		b.hp = hairpinPosition::unknown;
		return;
	}
	/**! get name of the hairpin **/
	std::string getName () const {
		return name;
	}
	/**! get total size of bed2 in the internal vector **/
	unsigned int size () const {
		return beds.size ();
	}
	/**! output all the bed2s in the internal vector, adjusting their 5' end and 3' end by to the relative distance to the annoated mature **/
	std::ostream& outputAll (std::ostream& os) const {
		for (const auto bed : beds) {
			os << bed.hairpinName << '\t';
			if (bed.hp==hairpinPosition::fivePrimeArm) {
				os << bed.start - fivePrimeStart << '\t' ;
				os << bed.end - fivePrimeEnd << '\t' ;
			} else if (bed.hp==hairpinPosition::threePrimeArm) {
				os << bed.start - threePrimeStart << '\t' ;
				os << bed.end - threePrimeEnd << '\t' ;
			} else if (bed.hp==hairpinPosition::loop) {
				os << bed.start - fivePrimeEnd << '\t';
				os << bed.end - threePrimeStart << '\t';
			} else {
				os << '*' << '\t';
				os << '*' << '\t';
			}
			os << bed.numOfSeq << '\t';
			os << bed.hp << "\t+\t";
			os << bed.sequence << '\n';
		}
		return os;
	}
	/**! calculate the 5' and 3' heterogeneity **/
	void calculateHeterogeneity () {

		for (const auto& mir : beds ) {
			if (mir.hp == hairpinPosition::fivePrimeArm) {
				fivePrimeHeterogeneityOfFivePrimeArm += std::abs (mir.start - fivePrimeStart) * mir.numOfSeq;
				threePrimeHeterogeneityOfFivePrimArm += std::abs (mir.end - fivePrimeEnd) * mir.numOfSeq;
				fiveArmCounter += mir.numOfSeq;
			}
			else if (mir.hp == hairpinPosition::threePrimeArm) {
				fivePrimeHeterogeneityOfThreePrimeArm += std::abs (mir.start - threePrimeStart) * mir.numOfSeq;
				threePrimeHeterogeneityOfThreePrimArm += std::abs (mir.end - threePrimeEnd) * mir.numOfSeq;
				threeArmCounter += mir.numOfSeq;
			}
		}
		if (fiveArmCounter) {
			fivePrimeHeterogeneityOfFivePrimeArm /= fiveArmCounter;
			threePrimeHeterogeneityOfFivePrimArm /= fiveArmCounter;
		}
		if (threeArmCounter) {
			fivePrimeHeterogeneityOfThreePrimeArm /= threeArmCounter;
			threePrimeHeterogeneityOfThreePrimArm /= threeArmCounter;
		}
	}
	double getFivePrimeHeterogeneityOfFivePrimeArm () const {
		return fivePrimeHeterogeneityOfFivePrimeArm;
	}
	double getThreePrimeHeterogeneityOfFivePrimArm () const {
		return threePrimeHeterogeneityOfFivePrimArm;
	}
	double getFivePrimeHeterogeneityOfThreePrimeArm () const {
		return fivePrimeHeterogeneityOfThreePrimeArm;
	}
	double getThreePrimeHeterogeneityOfThreePrimArm () const {
		return threePrimeHeterogeneityOfThreePrimArm;
	}
	int getFivePrimeCounter () const {
		return fiveArmCounter;
	}
	int getThreePrimeCounter () const {
		return threeArmCounter;
	}
};

#endif /* HAIRPIN_HPP_ */
