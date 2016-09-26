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
#include <unordered_map>
#include <vector>
#include <iomanip>
#include <thread>
#include <boost/program_options.hpp>
#include "piPipes_bed2_to_bedgraph.hpp"

using namespace std;

const string watson = "Watson.bedGraph";
const string crick = "Crick.bedGraph";

int main(int argc, char** argv) {
    string usage = R"(
This program convert BED2 format to bedGraph format.
=========
BED2 format is defined in the piPipes pipeline developed 
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
    string inputBED2{};
    string outputPrefix{};
    string chromInfor{};
    size_t nthread{};
    boost::program_options::options_description opts{usage};
    try {
        opts.add_options()
                ("help,h", "display this help message and exit")
                ("input,i", boost::program_options::value<string>(&inputBED2)->required(), "Input BED2 file;")
                ("chrom,c"
                 , boost::program_options::value<string>(&chromInfor)->required()
                 , "ChromInfo file, tab-delimited, with two fields; first one is the name of the ref, second one is the size")
                ("out_prefix,o"
                 , boost::program_options::value<string>(&outputPrefix)->default_value(string{"stdout"})
                 , "Prefix for the output bedGraph; Watson.bedGraph and Crick.bedGraph will be appended;")
                ("thread,p"
                 , boost::program_options::value<size_t>(&nthread)->default_value(1)
                 , "Number of thread to use; If available, it will use two threads to calculate Watson and Crick. More than 2 is unnecessary;");
        boost::program_options::variables_map vm;
        boost::program_options::store(boost::program_options::parse_command_line(argc, argv, opts), vm);
        boost::program_options::notify(vm);
        if (vm.count("help") || argc < 3) {
            cerr << opts << endl;
            exit(1);
        }
    }
    catch (exception& e) {
        cerr << "Error: " << e.what() << endl;
        cerr << opts << endl;
        exit(1);
    } catch (...) {
        cerr << "Unknown error!" << endl;
        cerr << opts << endl;
        exit(1);
    }

    ifstream fh_input_chrom{chromInfor};
    if (!fh_input_chrom) {
        cerr << "Error: cannot open file " << chromInfor << " for reading." << endl;
        exit(1);
    }
    if (outputPrefix == "stdout") {
        unsigned found = inputBED2.rfind(".bed2");
        if (found != string::npos)
            outputPrefix = inputBED2.substr(0, inputBED2.size() - found);
    }
    ofstream fh_watson{outputPrefix + watson};
    ofstream fh_crick{outputPrefix + crick};

    unordered_map<string, unsigned long> ref_sizes;
    string line;
    vector<string> fields;
    unsigned long tmp{0}, max_len{0};
    while (getline(fh_input_chrom, line)) {
        boost::split(fields, line, boost::is_any_of("\t"));
        tmp = stoul(fields[1]);
        ref_sizes[fields[0]] = tmp;
        if (tmp > max_len) max_len = tmp;
    }
    if (nthread > 1) {
        std::thread
            t1 = std::thread {bed2_bedgraph_converter<long double, '+'>  {inputBED2, &fh_watson, max_len, &ref_sizes}};
        std::thread
            t2 = std::thread {bed2_bedgraph_converter<long double, '-'>  {inputBED2, &fh_crick, max_len, &ref_sizes}};
        t1.join();
        t2.join();
    } else {
        bed2_bedgraph_converter<long double, '+'> wat_cvt{inputBED2, &fh_watson, max_len, &ref_sizes};
        wat_cvt();
        bed2_bedgraph_converter<long double, '-'> crk_cvt{inputBED2, &fh_crick, max_len, &ref_sizes};
        crk_cvt();
    }
    return 0;
}





