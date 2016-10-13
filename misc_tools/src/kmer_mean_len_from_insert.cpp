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
#include <unordered_map>
#include <boost/algorithm/string.hpp>

using namespace std;

void usage(char *p) {
    cerr << "usage: " << p
         << "\tinsert k " << endl;
}


int main(int argc, char **argv) {
    if (argc < 3) {
        usage(argv[0]);
        exit(1);
    }
    string fileName1{argv[1]};
    int k_mer = atoi(argv[2]);
    ifstream in1{fileName1};
    string line, prefix, seq, reads;
    string::iterator iter1, iter2;
    unordered_map<string, int> values_counts;
    unordered_map<string, int> species_counts;
    while (getline(in1, line)) {
        iter1 = iter2 = line.begin();
        for (int i = 0; i < k_mer; ++i) {
            ++iter2;
        }
        prefix = string {iter1, iter2};
        while (*++iter2 != '\t');
        seq = string {iter1, iter2};
        while (*iter2++ != '\t');
        reads = string {iter2, line.end()};
        values_counts[prefix] += stoi(reads) * seq.size();
        species_counts[prefix] += stoi(reads);
        // clog << seq << '\t' << prefix << '\t' << reads << endl;
    }
    for (const auto& p : values_counts) {
        cout << p.first << '\t' << ((double) p.second) / (double) species_counts[p.first] << endl;
    }
}
