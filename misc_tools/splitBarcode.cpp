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
#include <thread>
#include <iterator>
#include <unordered_map>
#include <boost/program_options.hpp>
#include <boost/iostreams/filtering_stream.hpp>
#include <boost/iostreams/device/file.hpp>
#include <boost/iostreams/filter/gzip.hpp>
#include <boost/iostreams/filter/bzip2.hpp>

using namespace std;

template<int OFFSET>
void
parseOneAdaptor (const string& inputFile, const string& bar, boost::iostreams::filtering_ostream* out, int maxMM) noexcept
{
	boost::iostreams::filtering_istream in;
	istream* pIn {new std::ifstream {inputFile}};
	char magic_number[4] = "\0\0\0";
	pIn->get(magic_number, 3);
	if (magic_number[0] == '\037' && magic_number[1] == (char)'\213') {in.push(boost::iostreams::gzip_decompressor());}
	else if (magic_number[0] == 'B' && magic_number[1] == 'Z') {in.push(boost::iostreams::bzip2_decompressor());}
	else if (magic_number[0] == '@') { }
	else {cerr << "unknown format" << endl; return;}
	pIn->seekg(0, pIn->beg);
	in.push(*pIn);
	int mm;
	string header, seq, qual;
	while (getline (in, header) && getline (in, seq) && in.ignore (10000, '\n') && getline (in, qual))
	{
		mm = 0; /// initiate mm counter
		auto barIter = bar.crbegin(); /// reverse const iterator pointing to the end of the barcode
		auto strIter = header.crbegin (); /// reverse const iterator pointing to the end of the header, where barcode is stored
		advance (strIter, OFFSET); /// bypassing \1 or \2
		while (barIter != bar.crend ())
		{ /// while barcode sequence has not been consumed
			if (*barIter != *strIter) ++mm;
			++barIter; ++strIter; /// increment two iterators
		}
		if (mm <= maxMM) *(out) << header << '\n' << seq  << "\n+\n" << qual << '\n';
	}
}

int main (int argc, char** argv)
{
	std::ios::sync_with_stdio(false);
	int offset;
	string inputFile;
	string prefix;
	int maxMM;
	boost::program_options::options_description opts (R"(
This program split barcode from TruSeq barcoded library in fastq format, in paired end format
Currently, it support 18 (R701-R718) barcodes currently used in the Zamore Lab (2015/04):
It also supports gzip/bzip2 as input and output gzip.
Please contact bo.han@Umassmed.edu for any questions.
)" 	);
	opts.add_options ()
	("help,h", "display this help message and exit")
	("input,i", boost::program_options::value<string>(&inputFile)->default_value("stdin"), "input fastq, can be plain/gzip/bzip2 format")
	("mismatch,m", boost::program_options::value<int>(&maxMM)->default_value(1), "number of mismatch allowed, [0,1,2], by default, 1 mismatch")
	("prefix,o", boost::program_options::value<string>(&prefix)->default_value("stdout"), "the prefix of the output files")
	("offset,s", boost::program_options::value<int>(&offset)->default_value(0), "distance of barcode to the right end of the line")
	;
	boost::program_options::variables_map vm;
	try
	{
		boost::program_options::store (boost::program_options::parse_command_line(argc, argv, opts), vm);
		boost::program_options::notify (vm);
	}
	catch(std::exception& e)
	{
		std::cerr << "Error: " << e.what() << std::endl;
		std::cerr << opts << std::endl;
		exit (1);
	}
	catch(...)
	{
		std::cerr << "Unknown error!" << std::endl;
		std::cerr << opts << std::endl;
		exit (1);
	}
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
		 "ATCACG"
		,"CGATGT"
		,"TTAGGC"
		,"TGACCA"
		,"GCCAAT"
		,"ACAGTG"
		,"CAGATC"
		,"ACTTGA"
		,"GATCAG"
		,"TAGCTT"
		,"GGCTAC"
		,"CTTGTA"
		,"AGTCAA"
		,"AGTTCC"
		,"ATGTCA"
		,"CCGTCC"
		,"GTAGAG"
		,"GTCCGC"
	};
	unordered_map<string, boost::iostreams::filtering_ostream*> bar2out;
	for (const auto bar: barcodes) {
		auto istret = bar2out.insert (make_pair (bar, new boost::iostreams::filtering_ostream));
		auto iter = istret.first;
		iter->second->push(boost::iostreams::gzip_compressor());
		iter->second->push(boost::iostreams::file_sink(prefix+bar+".fq.gz"));
	}
	vector<unique_ptr<std::thread> > threads;
	for (const auto & bar : barcodes) {
		switch(offset) {
			case 0: threads.emplace_back (new std::thread {parseOneAdaptor<0>, inputFile, bar, bar2out[bar], maxMM}); break;
			case 1: threads.emplace_back (new std::thread {parseOneAdaptor<1>, inputFile, bar, bar2out[bar], maxMM}); break;
			case 2: threads.emplace_back (new std::thread {parseOneAdaptor<2>, inputFile, bar, bar2out[bar], maxMM}); break;
			case 3: threads.emplace_back (new std::thread {parseOneAdaptor<3>, inputFile, bar, bar2out[bar], maxMM}); break;
			case 4: threads.emplace_back (new std::thread {parseOneAdaptor<4>, inputFile, bar, bar2out[bar], maxMM}); break;
			case 5: threads.emplace_back (new std::thread {parseOneAdaptor<5>, inputFile, bar, bar2out[bar], maxMM}); break;
			case 6: threads.emplace_back (new std::thread {parseOneAdaptor<6>, inputFile, bar, bar2out[bar], maxMM}); break;
			case 7: threads.emplace_back (new std::thread {parseOneAdaptor<7>, inputFile, bar, bar2out[bar], maxMM}); break;
			case 8: threads.emplace_back (new std::thread {parseOneAdaptor<8>, inputFile, bar, bar2out[bar], maxMM}); break;
			case 9: threads.emplace_back (new std::thread {parseOneAdaptor<9>, inputFile, bar, bar2out[bar], maxMM}); break;
			default: cerr << "offset not supported" << endl; break;
		}
	}
	for (auto & t : threads)
		if (t->joinable ())
			t->join();
	/**! close all opened output **/
	for (const auto& bar: barcodes) {
		delete bar2out[bar];
	}
}