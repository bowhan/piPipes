/*
 * This script is part of the piper, https://github.com/bowhan/Piper
 * Bo Han (bowhan@me.com)
 */

#include <iostream>
#include <fstream>
#include <unordered_map>
#include <vector>
#include <iomanip>
#include "boost/algorithm/string.hpp"
#include "boost/program_options.hpp"
#include "boost/thread.hpp"
using namespace std;

const string watson="Watson.bedGraph";
const string crick="Crick.bedGraph";

template <typename T, char S>
class bed2_bedgraph_converter {
private:
	string _in_file;
	ostream* _out;
	unsigned long _max_size;
	unordered_map<string, unsigned long>* _size_table;
public:
	bed2_bedgraph_converter (string in_file, ostream* out, unsigned long size, unordered_map<string, unsigned long>* size_table):
		_in_file {in_file}, _out {out}, _max_size {size}, _size_table {size_table} {
		}
		void operator () () {
			*_out << std::setprecision(5);
			*_out << std::fixed;
			ifstream* _in = new ifstream {_in_file};
			T* _container = new T [_max_size];
			T* p = _container;
			T previous_value {};
			unsigned long start_point {};
			for (unsigned long i = 0; i < _max_size; ++ i)
				*p = T{};
			string line;
			vector<string> fields;
			string chr;
			unsigned long pos;
			T num_of_reads;
			T num_of_mapping;
			getline (*_in, line);
			boost::split (fields, line, boost::is_any_of ("\t"));
			chr = fields[0];
			_in->seekg (ios_base::beg); /// C++11 says it also clear eofbit flag
			while (getline (*_in, line)) {
				boost::split (fields, line, boost::is_any_of ("\t"));
				if (fields[5][0]==S) { /// if this is the strand to do
					if ( !chr.empty() && chr!=fields[0]) { /// new chromosome
						/* begin to output */
						p = _container;
						previous_value = *p++;
						start_point = 0;
						unsigned long i {1};
						for (; i < (*_size_table)[chr]; ++ i) {
							if (*p != previous_value) { /// new value, need to output previous value
								if (previous_value) {
									if (S=='+') /// compiling time decision
										*_out << chr << '\t' << start_point << '\t' << i << '\t' << previous_value << '\n';
									if (S=='-')
										*_out << chr << '\t' << start_point << '\t' << i << "\t-" << previous_value << '\n';
								}
								start_point = i;
								previous_value = *p;
							}
							*p++ = T{};
						}
						if (previous_value) {
							if (S=='+') /// compiling time decision
								*_out << chr << '\t' << start_point << '\t' << i << '\t' << previous_value << '\n';
							if (S=='-')
								*_out << chr << '\t' << start_point << '\t' << i << "\t-" << previous_value << '\n';
						}
						/* end of output */
						for (; i < (*_size_table)[fields[0]]; ++i) /// continue to blank container
							*p++ = T{};
						chr = fields[0];
					}
					num_of_reads = boost::lexical_cast<T> (fields[3]);
					num_of_mapping = boost::lexical_cast<T> (fields[4]);
					if (fields[5]=="+")
						pos = stoul (fields[1]);
					else if (fields[5]=="-")
						pos = stoul (fields[2])-1;
					else {
						cerr << "Error: unknown strand " << fields[5] << endl;
						exit;
					}
					if (pos < (*_size_table)[chr])
						_container[pos] +=num_of_reads/num_of_mapping;
					else
						cerr << "Warning: chr size overflowing \n" << line << endl;
				}
			}
			p = _container;
			previous_value = *p++;
			start_point = 0;
			unsigned long i {1};
			for (; i < (*_size_table)[chr]; ++ i, ++p) {
				if (*p != previous_value) { /// new value, need to output previous value
					if (previous_value) {
						if (S=='+') /// compiling time decision
							*_out << chr << '\t' << start_point << '\t' << i << '\t' << previous_value << '\n';
						if (S=='-')
							*_out << chr << '\t' << start_point << '\t' << i << "\t-" << previous_value << '\n';
					}
					start_point = i;
					previous_value = *p;
				}
				//*p++ = T{}; no need to re-initialize
			}
			if (previous_value) {
				if (S=='+') /// compiling time decision
					*_out << chr << '\t' << start_point << '\t' << i << '\t' << previous_value << '\n';
				if (S=='-')
					*_out << chr << '\t' << start_point << '\t' << i << "\t-" << previous_value << '\n';
			}
			delete [] _container;
			_in->close ();
			delete _in;
		}
};

int main (int argc, char** argv) {
	string usage=R"(
This program convert BED2 format to bedGraph format.
=========
BED2 format is defined in the piper pipeline developed 
the Zamore and Weng labs to describe small RNA coordinates
with number_of_time_sequenced and number_of_times_mapped
been the 4th and 5th field of the normal BED format.
This script convert this format to bedGraph format, by 
incrementing $4/$5 for each line; And it will ONLY count
the 5' end, since this is designed for small RNA.

IMPORTANT: this scripts requires the input to be sorted by chr.

=========
it will automatically create Prefix.Watson.bedGraph and Prefix.Crick.bedGraph.

)";
	string inputBED2  {};
	string outputPrefix {};
	string chromInfor {};
	size_t nthread    {};
	boost::program_options::options_description opts {usage};
	try {
		opts.add_options ()
		("help,h", "display this help message and exit")
		("input,i", boost::program_options::value<string>(&inputBED2)->required(), "Input BED2 file;")
		("chrom,c", boost::program_options::value<string>(&chromInfor)->required(), "ChromInfo file, tab-delimited, with two fields; first one is the name of the ref, second one is the size")
		("out_prefix,o", boost::program_options::value<string>(&outputPrefix)->default_value(string{"stdout"}), "Prefix for the output bedGraph; Watson.bedGraph and Crick.bedGraph will be appended;")
		("thread,p", boost::program_options::value<size_t>(&nthread)->default_value(1), "Number of thread to use; If available, it will use two threads to calculate Watson and Crick. More than 2 is unnecessary;")
		;
		boost::program_options::variables_map vm;
		boost::program_options::store (boost::program_options::parse_command_line(argc, argv, opts), vm);
		boost::program_options::notify(vm);
		if (vm.count("help") || argc < 3)	{ cerr << opts << endl; exit (1); }
	}
	catch (exception& e) {
		cerr << "Error: " << e.what() << endl;
		cerr << opts << endl;
		exit (1);
	} catch (...) {
		cerr << "Unknown error!" << endl;
		cerr << opts << endl;
		exit (1);
	}

	ifstream fh_input_chrom {chromInfor};
	if (!fh_input_chrom) {
		cerr << "Error: cannot open file " << fh_input_chrom << " for reading." << endl;
		exit (1);
	}
	if (outputPrefix=="stdout") {
		unsigned found = inputBED2.rfind (".bed2");
		if (found != string::npos)
			outputPrefix = inputBED2.substr (0, inputBED2.size () - found);
	}
	ofstream fh_watson {outputPrefix + watson};
	ofstream fh_crick {outputPrefix + crick};

	unordered_map <string, unsigned long> ref_sizes;
	string line;
	vector<string> fields;
	unsigned long tmp {0}, max_len {0};
	while (getline (fh_input_chrom, line)) {
		boost::split (fields, line, boost::is_any_of ("\t"));
		tmp = stoul (fields[1]);
		ref_sizes[fields[0]] = tmp;
		if (tmp > max_len) max_len = tmp;
	}
	if (nthread>1) {
		boost::thread t1 = boost::thread { bed2_bedgraph_converter<long double,'+'>  {inputBED2, &fh_watson, max_len, &ref_sizes} };
		boost::thread t2 = boost::thread { bed2_bedgraph_converter<long double,'-'>  {inputBED2, &fh_crick, max_len, &ref_sizes} };
		t1.join ();
		t2.join ();
	} else {
		bed2_bedgraph_converter<long double,'+'> wat_cvt {inputBED2, &fh_watson, max_len, &ref_sizes};
		wat_cvt ();
		bed2_bedgraph_converter<long double,'-'> crk_cvt {inputBED2, &fh_crick, max_len, &ref_sizes};
		crk_cvt ();
	}
	return 0;
}





