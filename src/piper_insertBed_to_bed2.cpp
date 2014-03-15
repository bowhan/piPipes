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
#include <deque>
#include <unordered_map>

using namespace std;
int main (int argc, char** argv) {
	string usage=R"(
  This program takes two inputs:
  1. the insert file with two fields, first being the sequence and second being the number of times this seuqence been read;
  2. a standard BED file produced by bowtie reading the 1st file as input.
  and produces BED2 format, with 1-3 fields being name of chromosome, start and end, same as standard BED format.
  But the 4th column is the number of times this sequence been read, 5th column being the number of times this sequence being
  mapped. the 6th column is the sequence itself.
  usage:
  istBed2Bed2  input.insert  mapped.bed  >  output.bed2
  Please contact bowhan@me.com for any questions or bugs.
)";
	if (argc != 3) {
		cerr << usage;
		exit (1);
	}
	ifstream ist {argv[1]};
	ifstream bed {argv[2]};
	deque <int> istReads {};
	deque <string> istSeq {};
	string line {}, sequence {};
	while (getline (ist, line)) {
		auto iter = line.cbegin ();
		while (*++iter!= '\t') ;
		istSeq.emplace_back (line.cbegin(), iter);
		istReads.emplace_back (stoi (string(++iter, line.cend())));
	}
	auto istSize = istReads.size ();
	allocator<int> alloc;
	int* istNTM = alloc.allocate (istSize);
	int * p = istNTM;
	while (p!=istNTM+istSize) {
		*p++=0;
	}
	while (getline (bed, line)) {
		auto iter1 = line.cbegin ();
		while (*++iter1!='\t');
		while (*++iter1!='\t');
		while (*++iter1!='\t');
		auto iter2 = ++iter1;
		while (*++iter1!='\t');
		int n = stoi (string {iter2, iter1});
		if (istSize < n) {
			cerr << "Error: column 4 of " << line << " is larger than the lines of insert file..." << endl;
			exit (1);
		}
		++istNTM[n];
	}
	bed.clear ();
	bed.seekg (0, ios::beg);
	while (getline (bed, line)) {
		auto iter1 = line.cbegin ();
		while (*++iter1!='\t');
		while (*++iter1!='\t');
		while (*++iter1!='\t');
		auto iter2 = ++iter1;
		cout << string {line.cbegin(), iter1};
		while (*++iter1!='\t');
		int n = stoi (string {iter2, iter1});
		while (*++iter1!='\t');
		cout << istReads[n] << '\t' << istNTM[n] << '\t' << *(++iter1) << '\t' << istSeq[n] << '\n';
	}
	alloc.deallocate(istNTM, istSize);
}


