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

#include "../include/chiSquared.hpp"

using namespace std;

const string usage = R"(

A program to calculate the chiSquared value from a
tab-delimited file. The first column of this file
should be the name of the values and the rest be
the values in different samples. It can have optional
number of fields.
usage:
 chiSquared file

contact bo.han@umassmed.edu for questions

)";

int main(int argc, char **argv) {
    if (argc != 2) {
        cerr << usage << endl;
        exit(1);
    }
    ifstream in{argv[1]};
    string firstLine;
    getline(in, firstLine);
    uint32_t N = 0;
    for (char c : firstLine) {
        if (c == '\t')
            ++N;
    }
    in.clear();
    in.seekg(ios_base::beg); /// return to the beginning
    chiSquared<double> cs{N};
    cs.readFile(&in);
}