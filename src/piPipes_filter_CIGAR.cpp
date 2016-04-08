/* 
 # piPipes, a set of pipelines for PIWI-interacting RNA (piRNA) and transposon analysis
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
#include <set>
#include <boost/program_options.hpp>
#include "piPipes_sam.hpp"

using namespace std;
int main (int argc, char** argv) {
	std::string usage = R"(
This program filter the CIGAR field of a sam file; and is part of the piPipes package developped in the Zamore and Zlab in UMASS Med School. 
Contact Bo.Han@umassmed.edu or piPipesQ@gmail.com for questions, suggestions and bugs. 
Thank you!
)";
	/** options **/
	std::string inputFileName;
	std::string outputFileName;
	bool filter5, filter3;
	bool rev;
	uint32_t minLen, maxLen;
	char S;
	boost::program_options::options_description opts {usage};
	/** parsing cmdline **/
	try {
		opts.add_options ()
		("help,h", "Display This Help Message And Exit;")
		("input,i", boost::program_options::value<string>(&inputFileName)->default_value(std::string("stdin")), "Input file, using stdin if not specified")
		("output,o", boost::program_options::value<string>(&outputFileName)->default_value(std::string("stdout")), "Output file, using stdout if not specified")
		("filter5,5", boost::program_options::bool_switch(&filter5)->default_value(false), "Whether to filter out the sam entry with soft-clip at the 5' end; if --reverse flag is set, keep it")
		("filter3,3", boost::program_options::bool_switch(&filter3)->default_value(false), "Whether to filter out the sam entry with soft-clip at the 3' end; if --reverse flag is set, keep it")
		("target,t", boost::program_options::value<char>(&S)->default_value('S'), "The target CIGAR letter for filtering, use S by default")
		("reverse,r", boost::program_options::bool_switch(&rev)->default_value(false), "Whether to reverse the behavior: if this is on, KEEP the reads satisfying the condition instead of DUMPING them")
		("minlen,l", boost::program_options::value<uint32_t>(&minLen)->default_value(1), "The minimal length of the CIGAR target to keep, only valid when combined with r")
		("maxlen,L", boost::program_options::value<uint32_t>(&maxLen)->default_value(3), "The maximal length of the CIGAR target to keep, only valid when combined with r")
		;
		boost::program_options::variables_map vm;
		boost::program_options::store (boost::program_options::parse_command_line(argc, argv, opts), vm);
		boost::program_options::notify (vm);
		if (vm.count("help") || argc < 1)	{ std::cerr << opts << std::endl; exit (1); }
	} catch (std::exception& e) {
		std::cerr << "Error: " << e.what() << std::endl;
		std::cerr << opts << std::endl;
		exit (1);
	} catch (...) {
		std::cerr << "Unknown error!" << std::endl;
		std::cerr << opts << std::endl;
		exit (1);
	}
	if (!filter5 && !filter3) {
		std::cerr << "Error: Have to specify at least one of filter5 or filter3" << std::endl;
		std::cerr << opts << std::endl;
		exit (1);
	}
	if (minLen > maxLen && rev) {
		std::cerr << "Error: nothing will pass if --minlen is larger than --maxlen" << std::endl;
		std::cerr << opts << std::endl;
		exit (1);
	}
	if (!rev && (minLen != 1 || maxLen != 3)) {
		std::cerr << "Warning: --minlen and --maxlen are only valid when combined with --reverse" << std::endl;
	}
	std::set<char> validCIGAR {'M','I','D','N','S','H','P','=','X'};
	if (validCIGAR.find (S) == validCIGAR.end ()) {
		std::cerr << "Warning: the CIGAR letter you specified " << S << " is NOT a valid CIGAR letter. Valid CIGAR letters includes: " << std::endl;
		for (auto x : validCIGAR)
			std::cerr << x << '\n';
	}
	/** end of cmdline parsing **/

	/** output **/
	std::ostream* out1 {&std::cout};
	if (outputFileName!="stdout" && outputFileName!="-") {
		out1 = new std::ofstream {outputFileName};
	}
	string line;
	std::istream* in {&std::cin};
	if (inputFileName != "stdin" && inputFileName != "-" ) {
		in = new std::ifstream {inputFileName};
	}
	/* reading head */
	while (in->peek ()=='@' && getline (*in, line)) {
		*out1 << line << '\n';
	}
	/* reading reads */
	while (getline (*in, line)) {
		Sam sam {line};
		auto cigars = sam.parseCIGAR ();
		auto& l = cigars.front (); /// the left end of CIGAR
		auto& r = cigars.back (); /// the right end of CIGAR
		if (filter5) { /// trim 5' end
			if (sam.checkFlag(Sam::SAM_FLAG::REVERSE_COMPLEMENTED)) { /// 5' end is at the right of the CIGAR
				if (r.first == S) { // found the target
					if (rev) {
						if (r.second >= minLen && r.second <= maxLen)
							*out1 << line << '\n';
						continue;
					}
					else
						continue;
				}
			} else {
				if (l.first == S) { // found the target
					if (rev) {
						if (l.second >= minLen && l.second <= maxLen)
							*out1 << line << '\n';
						continue;
					}
					else
						continue;
				}
			}
		}
		if (filter3) {
			if (sam.checkFlag(Sam::SAM_FLAG::REVERSE_COMPLEMENTED)) { /// 3' end is at the left of the CIGAR
				if (l.first == S) { // found the target
					if (rev) {
						if (l.second >= minLen && l.second <= maxLen)
							*out1 << line << '\n';
						continue;
					}
					else
						continue;
				}
			} else {
				if (r.first == S) { // found the target
					if (rev) {
						if (r.second >= minLen && r.second <= maxLen)
							*out1 << line << '\n';
						continue;
					}
					else
						continue;
				}
			}
		}
		if (!rev)
			*out1 << line << '\n';
	}
	/** close instream if it is not connecting to the stdin **/
	if (in != &std::cin) {
		static_cast<std::ifstream*> (in)->close ();
		delete in;
	}
	/** close output stream if it is not connecting to the stdout **/
	if (out1 != &std::cout) {
		static_cast<std::ofstream*>(out1)->close ();
		delete out1;
	}
	return 0;
}


