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
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <unordered_map>
#include <thread>
#include <mutex>
#include <boost/program_options.hpp>
#include <boost/optional.hpp>
#include "piPipes_ping_pong.hpp"

#define DEFAULT_DISTANCE 100

using namespace std;

int main(int argc, char *argv[])
{
    boost::program_options::options_description opts(R"(
This program calculate local phasing scores from bed2 format.
BED2 format modifies BED format in two ways: 
  Firstly, it replace the name field with number of times this read is been sequenced in the library. 
  Secondly, the score field has been replaced by the number of times this read been mapped to the genome. 
  a.k.a, how many other loci, including itself, has the same sequence.
)");
    /** options **/
    string file_a;
    string file_b;
    int num_of_threads;
    int upper_limit;
    try {
        opts.add_options()
            ("help,h", "display this help message and exit")
            ("file_a,a",
             boost::program_options::value<string>(&file_a)->required(),
             "first file used for phasing calculation in BED2 format")
            ("file_b,b",
             boost::program_options::value<string>(&file_b)->required(),
             "second file used for phasing calculation in BED2 format")
            ("thread,p",
             boost::program_options::value<int>(&num_of_threads)->default_value(1),
             "the number of threads to use")
            ("upper_limit,u",
             boost::program_options::value<int>(&upper_limit)->default_value(DEFAULT_DISTANCE),
             "the maximal steps to go when calculating phasing scores");
        boost::program_options::variables_map vm;
        boost::program_options::store(boost::program_options::parse_command_line(argc, argv, opts), vm);
        boost::program_options::notify(vm);
        if (vm.count("help")) {
            std::cerr << opts << std::endl;
            exit(1);
        }
    } catch (std::exception &e) {
        std::cerr << "Error: " << e.what() << std::endl;
        std::cerr << opts << std::endl;
        exit(1);
    } catch (...) {
        std::cerr << "Unknown error!" << std::endl;
        std::cerr << opts << std::endl;
        exit(1);
    } /** end of cmdline parsing **/
    pptype poola, poolb;
    if (num_of_threads > 1) {
        thread parse_file_A(ParseBed, file_a, ref(poola));
        if (file_a != file_b)
            ParseBed(file_b, poolb);
        if (parse_file_A.joinable())
            parse_file_A.join();
    } else {
        ParseBed(file_a, poola);
        if (file_a != file_b)
            ParseBed(file_b, poolb);
    }
    multi_threading_queue ppt{upper_limit};
    double* z_score = new double[upper_limit];
    for (int i = 0; i < upper_limit; ++i) z_score[i] = 0.0;
    vector<thread> threads;
    for (int i = 0; i < num_of_threads; ++i) {
        if (file_a == file_b) {
            threads.emplace_back(PingPongPlayer<pptype>{poola, poola, ppt, z_score, &Phasing});
        } else {
            threads.emplace_back(PingPongPlayer<pptype>{poola, poolb, ppt, z_score, &Phasing});
        }
    }
    for (auto& thread: threads) {
        if (thread.joinable()) thread.join();
    }
    // output
    for (int k = 0; k < upper_limit; ++k) {
        cout << k + 1 << '\t' << z_score[k] << '\n';
    }
    delete[] z_score;
}
