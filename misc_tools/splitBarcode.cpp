/*
# Copyright (C) 2014  Bo Han, Wei Wang, Zhiping Weng, Phillip Zamore
#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 3 of the License, or
# (at your option) any later version.

# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.

# You should have received a copy of the GNU General Public License along
# with this program; if not, write to the Free Software Foundation, Inc.,
# 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
*/

#include <iostream>
#include <fstream>
#include <string>
#include <memory>
#include <iterator>
#include <unordered_map>
#include <boost/program_options.hpp>
#include <boost/thread.hpp>
#include "gzstream.h"
using namespace std;
template<int OFFSET>
void parseOneAdaptor (string inputFile, const string bar, shared_ptr<ofstream> out, int maxMM) {
	istream* in1 {new ifstream {inputFile}};
	char firstChar = in1->peek ();
	if (firstChar == 0x1f) {
		delete in1;
		in1 = new igzstream (inputFile.c_str(), ios::in); /// reopen file using igzstream
	}
	int mm; // running mismatch recorder
	string header,  /// string record the header of each fastq
		seq,  /// string record sequence of each fastq
		qual; /// string record quality of each fastq
	while (getline (*in1, header) && getline (*in1, seq) && in1->ignore (10000, '\n') && getline (*in1, qual)) {
		mm = 0; /// initiate mm counter
		auto barIter = bar.crbegin(); /// reverse const iterator pointing to the end of the barcode
		auto strIter = header.crbegin (); /// reverse const iterator pointing to the end of the header, where barcode is stored
		advance (strIter, OFFSET); /// bypassing \1 or \2
		while (barIter != bar.crend ()) { /// while barcode sequence has not been consumed
			if (*barIter != *strIter) /// record mismatches
				++mm;
			++barIter; ++strIter; /// increment two iterators
		}
		if (mm <= maxMM) /// once pass the threshold, record it in the right place
			*(out) << header << '\n' << seq  << "\n+\n" << qual << '\n';

	}
	delete in1;
	out->close ();
}

int main (int argc, char** argv)
{
	
	int offset;
	string inputFile;
	string prefix;
	int maxMM;
	boost::program_options::options_description opts (R"(
This program split barcode from TruSeq barcoded library in fastq format, in paired end format
Currently, it support eight barcodes currently used in the Zamore Lab (2013/05):
	id1=ATCACG
	id2=CGATGT
	id3=TTAGGC
	id4=TGACCA
	id5=ACAGTG
	id6=GCCAAT
	id7=CAGATC
	id8=ACTTGA
updated on 09/2014.
@HWI-ST570:130:C4KR4ACXX:5:1101:1193:2329 1:N:0:CGATGTAT	  
The script will automatically determine whether the input file is in plain format or gzipped. 
Please contact bo.han@Umassmed.edu for any questions. 
)" 	);
	opts.add_options ()
		("help,h", "display this help message and exit")
		("input,i", boost::program_options::value<string>(&inputFile)->default_value("stdin"), "input fastq, can be plain format or gzipped format")
		("mismatch,m", boost::program_options::value<int>(&maxMM)->default_value(1), "number of mismatch allowed, [0,1,2], by default, 1 mismatch")
		("prefix,o", boost::program_options::value<string>(&prefix)->default_value("stdout"), "the prefix of the output file")
		("offset,s", boost::program_options::value<int>(&offset)->default_value(0), "distance of barcode to the right end of the line")
		;
	boost::program_options::variables_map vm;
	boost::program_options::store (boost::program_options::parse_command_line(argc, argv, opts), vm);
	boost::program_options::notify (vm);
	if (argc < 3 || vm.count("help"))	{
		cerr << opts << endl;
		exit (1);
	}
	/**! check input **/
	if (vm.find ("input")==vm.end()) {
		cerr << "Error: No input file specified" << endl;
		cerr << opts << endl;
		exit (1);
	}
	
	if  (prefix=="stdout") prefix = inputFile + ".barcode.";
	if (prefix.back () != '.') prefix+='.';
	/**! declare barcodes **/
	vector<string> barcodes {
		"ATCACG",
		"CGATGT",
		"TTAGGC",
		"TGACCA",
		"GCCAAT",
		"ACAGTG",
		"CAGATC",
		"ACTTGA"
	};
	unordered_map<string, shared_ptr<ofstream> > bar2out;
	for (const auto bar: barcodes) {
		bar2out.insert (make_pair (bar, unique_ptr<ofstream> {new ofstream {prefix+bar+".fq"}}));
	}
	vector<unique_ptr<boost::thread> > threads;
	clog << "offset: " << offset << endl;
	for (const auto & bar : barcodes) {
		switch(offset) {
			case 0: threads.emplace_back (new boost::thread {parseOneAdaptor<0>, inputFile, bar, bar2out[bar], maxMM}); break;
			case 1: threads.emplace_back (new boost::thread {parseOneAdaptor<1>, inputFile, bar, bar2out[bar], maxMM}); break;
			case 2: threads.emplace_back (new boost::thread {parseOneAdaptor<2>, inputFile, bar, bar2out[bar], maxMM}); break;
			case 3: threads.emplace_back (new boost::thread {parseOneAdaptor<3>, inputFile, bar, bar2out[bar], maxMM}); break;
			case 4: threads.emplace_back (new boost::thread {parseOneAdaptor<4>, inputFile, bar, bar2out[bar], maxMM}); break;
			case 5: threads.emplace_back (new boost::thread {parseOneAdaptor<5>, inputFile, bar, bar2out[bar], maxMM}); break;
			case 6: threads.emplace_back (new boost::thread {parseOneAdaptor<6>, inputFile, bar, bar2out[bar], maxMM}); break;
			case 7: threads.emplace_back (new boost::thread {parseOneAdaptor<7>, inputFile, bar, bar2out[bar], maxMM}); break;
			case 8: threads.emplace_back (new boost::thread {parseOneAdaptor<8>, inputFile, bar, bar2out[bar], maxMM}); break;
			case 9: threads.emplace_back (new boost::thread {parseOneAdaptor<9>, inputFile, bar, bar2out[bar], maxMM}); break;
			default: cerr << "offset not supported" << endl; break;		
		}
	}
	for (auto & t : threads)
		if (t->joinable ())
			t->join();
	/**! close all opened output **/
	for (const auto& bar: barcodes) {
		(bar2out[bar])->close ();
	}
}
