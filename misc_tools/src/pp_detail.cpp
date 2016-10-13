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
#include <iomanip>
#include <unordered_map>
#include <boost/program_options.hpp>
#include <boost/algorithm/string.hpp>

using namespace std;

void read_file_into_unordered_map(const string& file_name, unordered_map<string, unordered_map<char, unordered_map<
    uint64_t
    , double>>> *umap_five
                                 ) {
    ifstream in(file_name);
    if (!in) {
        cerr << "file " << file_name << " cannot be opened\n";
        exit(1);
    }
    char buffer[1024];
    char *saveptr;
    string chr;
    uint64_t num_of_reads, num_of_origin;
    uint64_t start, end;
    char strand;
    while (in.getline(buffer, 1024, '\n')) {
        char *saveptr;
        chr = strtok_r(buffer, "\t", &saveptr);
        start = strtoul(strtok_r(NULL, "\t", &saveptr), NULL, 0);
        end = strtoul(strtok_r(NULL, "\t", &saveptr), NULL, 0);
        num_of_reads = strtoul(strtok_r(NULL, "\t", &saveptr), NULL, 0);
        num_of_origin = strtoul(strtok_r(NULL, "\t", &saveptr), NULL, 0);
        strand = strtok_r(NULL, "\t", &saveptr)[0];
        if (strand == '+') {
            (*umap_five)[chr][strand][start] += (double(num_of_reads)) / double(num_of_origin);
            // (*umap_three) [chr][strand][end-1] += (double(num_of_reads))/double(num_of_origin);
        } else {
            (*umap_five)[chr][strand][end - 1] += (double(num_of_reads)) / double(num_of_origin);
            // (*umap_three) [chr][strand][start] += (double(num_of_reads))/double(num_of_origin);
        }
    }
}

int main(int argc, char *argv[]) {
    boost::program_options::options_description opts(R"(
this program calculate the ping-pong score for each genomic posisitons.
)");
    /** options **/
    string file_a;
    string file_b;
    string chrom_file;
    int distance;
    int num_of_threads;
    int upper_limit;
    try {
        opts.add_options()
                ("help,h", "display this help message and exit")
                ("file_a,a"
                 , boost::program_options::value<string>(&file_a)->default_value("stdin")
                 , "first file used for ping-pong in BED2 format")
                ("file_b,b"
                 , boost::program_options::value<string>(&file_b)->default_value("stdout")
                 , "second file used for ping-pong in BED2 format")
                ("chrom,c", boost::program_options::value<string>(&chrom_file)->required(), "chromosome size file")
                ("distance,d", boost::program_options::value<int>(&distance)->default_value(10), "Ping-Pong distance");
        boost::program_options::variables_map vm;
        boost::program_options::store(boost::program_options::parse_command_line(argc, argv, opts), vm);
        boost::program_options::notify(vm);
        if (vm.count("help")) {
            std::cerr << opts << std::endl;
            exit(1);
        }
    } catch (std::exception& e) {
        std::cerr << "Error: " << e.what() << std::endl;
        std::cerr << opts << std::endl;
        exit(1);
    } catch (...) {
        std::cerr << "Unknown error!" << std::endl;
        std::cerr << opts << std::endl;
        exit(1);
    } /** end of cmdline parsing **/
    unordered_map<string, uint64_t> chrom_size;
    clog << "reading chromosome file\n";
    ifstream ifs{chrom_file};
    if (!ifs) {
        cerr << "error, cannot open chromosome size file " << chrom_file << endl;
        exit(1);
    }
    string temp_str;
    while (getline(ifs, temp_str, '\n')) {
        string chr;
        uint64_t chr_size;
        auto str_iter = temp_str.begin();
        while (*++str_iter != '\t' && str_iter != temp_str.end());
        chr = string{temp_str.begin(), str_iter};
        auto str_iter2 = ++str_iter;
        while (*++str_iter != '\t' && str_iter != temp_str.end());
        chr_size = stoul(string{str_iter2, str_iter});
        if (chrom_size.find(chr) == chrom_size.end()) {
            chrom_size.insert(make_pair(chr, chr_size));
        } else {
            cerr << "error: duplicated chr: " << chr << endl;
            exit(1);
        }
    }
    unordered_map<string, unordered_map<char, unordered_map<uint64_t, double>>> *chrPos5endValue1;
    unordered_map<string, unordered_map<char, unordered_map<uint64_t, double>>> *chrPos5endValue2;
    chrPos5endValue1 = new unordered_map<string, unordered_map<char, unordered_map<uint64_t, double>>>;
    clog << "reading bed2 file a\n";
    read_file_into_unordered_map(file_a, chrPos5endValue1);
    if (file_a != file_b) {
        clog << "reading bed2 file b\n";
        chrPos5endValue2 = new unordered_map<string, unordered_map<char, unordered_map<uint64_t, double>>>;
        read_file_into_unordered_map(file_b, chrPos5endValue2);
    } else {
        chrPos5endValue2 = chrPos5endValue1;
        // chrPos3endValue2 = chrPos3endValue1;
    }
    ifstream ifa{file_a};
    vector<string> vec_;
    string chr;
    uint64_t _left;
    uint64_t _right;
    uint64_t _nor;
    int _ntm;
    double _score;
    double _pp_score;
    double _phasing_score;
    std::cout << std::fixed << std::setprecision(0);
    clog << "reading bed2 file a and output pp and phasing score\n";
    --distance; //
    while (getline(ifa, temp_str, '\n')) {
        _pp_score = 0.;
        _phasing_score = 0.;
        boost::split(vec_, temp_str, boost::is_any_of("\t"));
        chr = vec_[0];
        _left = stoul(vec_[1]);
        _right = stoul(vec_[2]) - 1;
        _nor = stoul(vec_[3]);
        _ntm = stoi(vec_[4]);
        _score = double(_nor) / _ntm;
        if (vec_[5] == "+") {
            if (_left + distance < chrom_size[chr])
                _pp_score = _score * (*chrPos5endValue2)[chr]['-'][_left + distance];
            if (_right + distance + 1 < chrom_size[chr])
                _phasing_score = min(_score, (double) (*chrPos5endValue2)[chr]['+'][_right + 1]);
        } else {
            if (_right - distance - 1 >= 0) _pp_score = _score * (*chrPos5endValue2)[chr]['+'][_right - distance];
            if (_left - 1 >= 0) _phasing_score = min(_score, (double) (*chrPos5endValue2)[chr]['-'][_left - 1]);
        }
        cout << temp_str << '\t' << _pp_score << '\t' << _phasing_score << '\n';
    }
    if (chrPos5endValue1 != chrPos5endValue2) { delete chrPos5endValue2; }
    delete chrPos5endValue1;
}
