/* 
# piper, a pipeline collection for PIWI-interacting RNA (piRNA) and transposon analysis
# Copyright (C) <2014>  <Bo Han, Wei Wang, Phillip Zamore, Zhiping Weng>
#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2 of the License, or
# (at your option) any later version.

# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.

# You should have received a copy of the GNU General Public License along
# with this program; if not, write to the Free Software Foundation, Inc.,
# 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
*/

#ifndef SAM_HPP_
#define SAM_HPP_

#include <iostream>
#include <string>
#include <vector>
#include <bitset>
#include <algorithm>
#include <iterator>
#include <unordered_map>
#include <boost/lexical_cast.hpp>
#include <boost/algorithm/string.hpp>


class Sam
{
public:
	enum SAM_FLAG {
		PAIRED_END = 0,
		EACH_END_ALIGNED = 1,
		UNMAPPED = 2,
		NEXT_UNMAPPED = 3,
		REVERSE_COMPLEMENTED = 4,
		NEXT_REVERSE_COMPLEMENTED = 5,
		FIRST_SEG = 6,
		SECOND_SEG = 7,
		SECONDARY_ALIGNMENT = 8,
		NOT_PASSING_QUALITY = 9,
		PCR_DUP = 10,
		FLAG_SIZE = 11
	};
private:
	std::string QNAME = "";
	std::bitset<SAM_FLAG::FLAG_SIZE> FLAG;
	std::string RNAME = "";
	uint64_t POS = 0;
	int MAPQ = 255;
	std::string CIGAR = "";
	std::string RNEXT = "";
	uint64_t PNEXT = 0;
	int64_t TLEN = 0;
	std::string SEQ = "";
	std::string QUAL = "";
	std::unordered_map<std::string, std::string> OPTIONAL_FIELDS ;
public:
	/**! default constructor **/
	Sam () = default ;
	/**! ctor from a line **/
	Sam (const std::string& line)
	{
		std::string::const_iterator iter1 {line.cbegin ()}, iter2 {line.cbegin ()};
		while (*++iter2 != '\t');
		QNAME = std::string {iter1, iter2};
		iter1 = ++iter2;

		while (*++iter2 != '\t');
		FLAG = std::bitset<SAM_FLAG::FLAG_SIZE> {std::stoul (std::string {iter1, iter2})};
		iter1 = ++iter2;

		while (*++iter2 != '\t');
		RNAME = std::string {iter1, iter2};
		iter1 = ++iter2;

		while (*++iter2 != '\t');
		POS = std::stoul (std::string {iter1, iter2});
		iter1 = ++iter2;

		while (*++iter2 != '\t');
		MAPQ = std::stoi (std::string {iter1, iter2});
		iter1 = ++iter2;

		while (*++iter2 != '\t');
		CIGAR = std::string {iter1, iter2};
		iter1 = ++iter2;

		while (*++iter2 != '\t');
		RNEXT = std::string {iter1, iter2};
		iter1 = ++iter2;

		while (*++iter2 != '\t');
		PNEXT = std::stoul (std::string {iter1, iter2});
		iter1 = ++iter2;

		while (*++iter2 != '\t');
		TLEN = std::stol (std::string {iter1, iter2});
		iter1 = ++iter2;

		while (*++iter2 != '\t');
		SEQ = std::string {iter1, iter2};
		iter1 = ++iter2;

		while (*++iter2 != '\t');
		QUAL = std::string {iter1, iter2};
		std::string Tag, Value;
		while (iter2 != line.cend ())
		{
			iter1 = ++iter2;
			iter2 += 2;
			Tag = std::string {iter1, iter2};
			iter1 += 5;
			while ( (*iter2 != '\t') && (iter2 != line.cend()))
				++iter2;
			Value = std::string {iter1, iter2};
			OPTIONAL_FIELDS.insert (std::make_pair (Tag, Value));
		}
	}
	/**! copy ctor **/
	Sam (const Sam&) = default;
	/**! move ctor **/
	Sam (Sam&&) =default;
	/**! lvalue assignment operator **/
	Sam& operator=(const Sam& other)
	{
		if (this == &other) return *this;
		QNAME = other.QNAME;
		FLAG  = other.FLAG;
		RNAME = other.RNAME;
		POS   = other.POS;
		MAPQ  = other.MAPQ;
		CIGAR = other.CIGAR;
		RNEXT = other.RNEXT;
		PNEXT = other.PNEXT;
		TLEN  = other.TLEN;
		SEQ   = other.SEQ;
		QUAL  = other.QUAL;
		OPTIONAL_FIELDS = other.OPTIONAL_FIELDS;
		return *this;
	}
	/**! rvalue assignment operator **/
	Sam& operator=(Sam&& other)
	{
		assert (this != &other);
		QNAME.swap (other.QNAME);
		FLAG = other.FLAG;
		RNAME.swap (other.RNAME);
		POS   = other.POS;
		MAPQ  = other.MAPQ;
		CIGAR.swap (other.CIGAR);
		RNEXT.swap (other.RNEXT);
		PNEXT = other.PNEXT;
		TLEN  = other.TLEN;
		SEQ.swap (other.SEQ);
		QUAL.swap (other.QUAL);
		OPTIONAL_FIELDS.swap (other.OPTIONAL_FIELDS);
		other.~Sam ();
		return *this;
	}
	/**! destructor **/
	~Sam () = default;

	/** parse CIGAR **/
	std::vector<std::pair <char, uint32_t> > parseCIGAR () const {
		std::vector<std::pair <char, uint32_t> > CIGARs;
		uint32_t t = 0;
		for (char x : CIGAR) {
			if (x >= '0'  && x <= '9') {
				t = 10*t + x - '0';
			} else {
				CIGARs.emplace_back (x, t);
				t = 0;
			}
		}
		return CIGARs;
	}

	template <typename T>
	bool checkFlag (T f) const
	{
		return FLAG.test (f);
	}

	template <typename T, typename ... Args>
	bool checkFlag(T f, Args ... rest) const
	{
		if (sizeof ...(rest))
			return FLAG.test(f) && checkFlag (rest...);
		return FLAG.test (f);
	}
	std::string getTag (const std::string& key) const
	{
		return OPTIONAL_FIELDS.at (key);
	}
	std::string mutation () const
	{
		return OPTIONAL_FIELDS.at ("MD");
	}
	std::string getSeq () const {
		return SEQ;
	}
	std::string getRevCompSeq () const {
		std::string rcSeq {SEQ.rbegin(), SEQ.rend()};
		for (auto & c : rcSeq) {
			switch (c) {
			case 'A': c = 'T'; break;
			case 'T': c = 'A'; break;
			case 'C': c = 'G'; break;
			case 'G': c = 'C'; break;
			default : throw "illegal char";
			}
		} /* end of RC */
		return rcSeq;
	}
	int getSize () const
	{
		return SEQ.size ();
	}

	std::string getQueryName () const {
		return QNAME;
	}
	friend std::ostream& operator<< (std::ostream& os, const Sam& sam)
	{
		os << sam.QNAME << '\t'
				<< sam.FLAG.to_string() << '\t'
				<< sam.RNAME << '\t'
				<< sam.POS << '\t'
				<< sam.MAPQ << '\t'
				<< sam.CIGAR << '\t'
				<< sam.RNEXT << '\t'
				<< sam.PNEXT << '\t'
				<< sam.TLEN << '\t'
				<< sam.SEQ << '\t'
				<< sam.QUAL;
		for (const auto& p : sam.OPTIONAL_FIELDS)
			os << '\t' << p.first << ':' << p.second;
		os << '\n';
		return os;
	}
};

#endif /* SAM_HPP_ */
