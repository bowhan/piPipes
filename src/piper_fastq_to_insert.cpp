/*
 * This script is part of the piper, https://github.com/bowhan/Piper
 * Bo Han (bowhan@me.com)
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



