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

#include <iostream>
#include <fstream>
#include <string>
#include <unordered_map>
#include <boost/algorithm/string.hpp>
#include <boost/program_options.hpp>

using namespace std;

template <bool Five>
inline void countBed2 (istream*, unordered_map<string, int>&, unordered_map<string, double* >&, unordered_map<string, double* >&, int );

int main (int argc, char** argv) {
	std::string usage = R"(
======================================================
This program takes bed2 file as input and summarize
the hits on each chromosomes. The chromsomes were divided
into different number of bins. 
Different from similar tools like bigWigSummary (kent),
this script doesn't require the input to be sorted.
As a consequence, it uses much more memory.

Please contact bo.han@Umassmed.edu for any questions. 
======================================================
)";
	/** options **/
	std::string inputBed2FileName;
	std::string chromSizeFileName;
	std::string outputFileName;
	int numOfBins {};
	bool onlyCountFivePrimeEnd {};
	boost::program_options::options_description opts {usage};
	/** parsing cmdline **/
	try {
		opts.add_options ()
		("help,h", "Display This Help Message And Exit;")
		("input,i", boost::program_options::value<std::string>(&inputBed2FileName)->default_value(std::string("stdin")), "Input BED2 file, using stdin by default; Also can be specified by stdin or -")
		("chrom,c", boost::program_options::value<std::string>(&chromSizeFileName)->required(), "Chrom file")
		("bin,b", boost::program_options::value<int>(&numOfBins)->default_value(1000), "Number of bins to devide the chromsome to")
		("output,o", boost::program_options::value<std::string>(&outputFileName)->default_value(std::string("stdout")), "Output file, using stdout by default; Also can be specified by stdout or +")
		("five,5", boost::program_options::bool_switch(&onlyCountFivePrimeEnd)->default_value(false), "Only count the five prime end?")
		;
		boost::program_options::variables_map vm;
		boost::program_options::store (boost::program_options::parse_command_line(argc, argv, opts), vm);
		boost::program_options::notify (vm);
		if (vm.count("help"))	{ std::cerr << opts << std::endl; exit (1); }
	} catch (std::exception& e) {
		std::cerr << "Error: " << e.what() << std::endl;
		std::cerr << opts << std::endl;
		exit (1);
	} catch (...) {
		std::cerr << "Unknown error!" << std::endl;
		std::cerr << opts << std::endl;
		exit (1);
	} /** end of cmdline parsing **/
	/** input stream specify **/
	std::istream* in1 {&std::cin};
	if (inputBed2FileName!="stdin" && inputBed2FileName!="-") {
		in1 = new std::ifstream {inputBed2FileName};
		if (!*in1) {
			std::cerr << "Error: cannot open file " << inputBed2FileName << " for reading." << std::endl;
			exit (1);
		}
	}
	/** output stream specify **/
	std::ostream* out1 {&std::cout};
	if (outputFileName!="stdout" && outputFileName!="-") {
		out1 = new std::ofstream {outputFileName};
		if (!*out1) {
			std::cerr << "Error: cannot open file " << outputFileName << " for writing." << std::endl;
			exit (1);
		}
	}
	/* global variables */
	string line;
	vector<string> container;
	container.reserve (2);
	/* read chrom size file */
	unordered_map<string, int> chromSizes;
	unordered_map<string, double* > WatsonCounter;
	unordered_map<string, double* > CrickCounter;
	ifstream chromFh {chromSizeFileName};
	if (!chromFh) {
		cerr << "Error: cannot open " << chromSizeFileName << " for reading" << endl;
		exit (1);
	}
	while (getline (chromFh, line, '\n')) {
		boost::split (container, line, boost::is_any_of ("\t"));
		chromSizes[container[0]] = stoi (container[1]);
		WatsonCounter[container[0]] = new double [numOfBins];
		CrickCounter[container[0]] = new double [numOfBins];
		for (int i=0; i<numOfBins;++i) {
			WatsonCounter[container[0]][i] = 0.0;
			CrickCounter[container[0]][i] = 0.0;
		}
	}
	/* read input bed2 file */
	if (onlyCountFivePrimeEnd)
		countBed2<true> (in1, chromSizes, WatsonCounter, CrickCounter, numOfBins);
	else
		countBed2<false> (in1, chromSizes, WatsonCounter, CrickCounter, numOfBins);
	for (const auto & chr_counts : WatsonCounter ) {
		for ( int i=0; i<numOfBins;++i)
			*out1 << chr_counts.first << '\t' << i+1 << '\t' << chr_counts.second[i] << "\t-" << CrickCounter[chr_counts.first][i] << '\n';
	}
	/** close input stream if it is not connecting to the stdin **/
	if (in1 != &std::cin) {
		static_cast<std::ifstream*>(in1)->close ();
		delete in1;
	}
	/** close output stream if it is not connecting to the stdout **/
	if (out1 != &std::cout) {
		static_cast<std::ofstream*>(out1)->close ();
		delete out1;
	}
	return 0;
}


template <bool Five>
inline void countBed2 (istream* in, unordered_map<string, int>& chromSizes, unordered_map<string, double* >& Watson, unordered_map<string, double* >& Crick, int numOfBins) {
	vector<string> container; container.reserve (7);
	string line;
	while (getline (*in, line)) {
		boost::split (container, line, boost::is_any_of ("\t"));
		if (chromSizes.find (container[0]) == chromSizes.end ()) {
			cerr << "Error: " << container[0] << " does't exist in the chrom size file" << endl;
			exit (1); // memory has not been released.
		}
		if (container[5]=="+") { // Watson
			if (Five) {
				Watson[container[0]][stoi(container[1])%numOfBins] += stod (container[3])/stod(container[4]);
			}
			if (!Five) {
				for (int i = stoi (container[1]); i<stoi (container[2]);++i)
					Watson[container[0]][i%numOfBins] += stod (container[3])/stod(container[4]);
			}
		} else { // Crick
			if (Five) {
				Crick[container[0]][(stoi(container[2])-1)%numOfBins] += stod (container[3])/stod(container[4]);
			}
			if (!Five){
				for (int i = stoi (container[1]); i<stoi (container[2]);++i)
					Crick[container[0]][i%numOfBins] += stod (container[3])/stod(container[4]);
			}
		}
	}
}



