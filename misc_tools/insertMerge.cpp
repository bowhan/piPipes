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
#include <vector>
#include <unordered_map>
#include <set>
#include <thread>
#include "boost/program_options.hpp"
#include "insertMerge.hpp"

using namespace std;

int main(int argc, char** argv)
{
	string usage=R"(
	This program merge several files based on a common column.
	Similar to merge function in R.
	)";
	vector<string> inputFileNames;
	string outputFileName;
	string defaultValue;
	int commonField = 0, mergeField = 0;
	boost::program_options::options_description opts {usage};
	opts.add_options()
	("help,h", "display help message")
	("inputs,i", boost::program_options::value<vector<string>>(&inputFileNames)->multitoken(), "input files")
	("output,o", boost::program_options::value<string>(&outputFileName)->default_value(string("stdout")), "output file")
	("common,c", boost::program_options::value<int>(&commonField)->default_value(1), "commond field to be used as key for merging, 1-based")
	("target,t", boost::program_options::value<int>(&mergeField)->default_value(0), "target field to be merged, 1-based")
	("default,d", boost::program_options::value<string>(&defaultValue), "the default target value to be used if this key does not exist")
	;
	boost::program_options::variables_map vm;
	try {
		boost::program_options::store(boost::program_options::parse_command_line(argc, argv, opts), vm);
		boost::program_options::notify(vm);
	} catch (std::exception& e) {
		cerr << "error: " << e.what() << endl;
		cerr << opts << endl;
		exit (1);
	} catch (...) {
		cerr << "unknown error" << endl;
		cerr << opts << endl;
		exit (1);
	}
	if (vm.count("help") || argc < 2) {
		cerr << opts << endl;
		exit (1);
	}
	if (commonField < 0 || mergeField < 0) {
		cerr << "-c and -t has to be positive value" << endl;
		exit (2);
	}
	std::ostream* out1 {& std::cout};
	if (outputFileName!="stdout" && outputFileName!="-") {
		out1 = new std::ofstream {outputFileName};
	}
	if (!*out1) {
		cerr << "cannot open output file " << outputFileName << endl;
		exit (1);
	}
	--commonField;
	--mergeField;
	size_t numberOfFiles = inputFileNames.size();
	unordered_map<string,string>* maps = new unordered_map<string,string>[numberOfFiles];
	vector<std::thread> threads (numberOfFiles);
	for (size_t i = 0; i < numberOfFiles; ++i) {
		unordered_map<string,string>* m = new unordered_map<string,string>;
		threads[i] = std::thread {readFile, boost::ref(inputFileNames[i]), m, commonField, mergeField};
	}
	for(auto& t : threads) {
        if(t.joinable()) {
            t.join();
        }
    }
	set<string> keys;
	for (size_t i = 0; i < numberOfFiles; ++i) {
		for (const auto& x : *(maps + i)) {
			keys.insert (x.first);
		}
	}
	for (const auto& key : keys) {
		*out1 << key;
		for (size_t i = 0; i < numberOfFiles; ++i) {
			auto found = (maps+i)->find(key);
			if (found == (maps+i)->end()) {
				*out1 << '\t' << defaultValue;
			} else {
				*out1 << '\t' << found->second;
			}
		}
		*out1 << '\n';
	}
	delete [] maps;
	if (out1 != &std::cout) {
		static_cast<std::ofstream*>(out1)->close();
		delete out1;
	}
}
