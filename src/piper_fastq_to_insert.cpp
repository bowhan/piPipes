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
#include <memory>
#include <unordered_map>
#include "gzstream.h"

using namespace std;

int main (int argc, char** argv) {
	string usage=R"(
  usage:
    fastqToInsert input.fq/input.fq.gz input.insert

)";

	if (argc != 3) {
		cerr << usage ;
		exit (1);
	}
	unique_ptr<istream> in {new ifstream {argv[1]}};

	char firstBit = in->peek ();

	if (firstBit == 0x1f) {
		in.reset (new igzstream (argv[1], ios::in));
	}
	if (!(*in)) {
		cerr << "error: cannot open file " << argv[1] << "for reading" << endl;
		exit (1);
	}
	ofstream  out {argv[2]};
	if (!out) {
		cerr << "error: cannot open file " << argv[2] << " for writting" << endl;
		exit (1);
	}
	string line, sequence;
	unordered_map<string, int> counter;
	while (in->ignore (10000,'\n') && getline (*in, sequence) && in->ignore (10000,'\n') && in->ignore (10000,'\n')) {
		++counter[sequence];
	}
	for (const auto & seq : counter) {
		out << seq.first << '\t' << seq.second << '\n';
	}
}



