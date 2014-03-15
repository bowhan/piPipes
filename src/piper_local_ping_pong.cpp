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
#include <vector>
#include <array>
#include <memory>
#include <unordered_map>
#include <boost/thread.hpp>
#include <boost/program_options.hpp>
#include <boost/ref.hpp>

using namespace std;
#define UPPERLIMIT 30
void do_ping_pong (
		const unordered_map<string, unordered_map<char, unordered_map<uint64_t, double>>>& ,
		const unordered_map<string, unordered_map<char, unordered_map<uint64_t, double>>>& ,
		int ,
		double* const
);
void read_file_into_unordered_map (const string& file_name, unordered_map<string, unordered_map<char, unordered_map<uint64_t, double>>>& umap)
{
	ifstream in (file_name);
	if (!in)
	{
		cerr << "file " << file_name << " cannot be opened\n";
		exit (1);
	}
	char buffer [1024];
	char* saveptr;
	string chr;
	uint64_t num_of_reads, num_of_origin;
	uint64_t start, end;
	char strand;
	while (in.getline (buffer, 1024, '\n'))
	{
		char* saveptr;
		chr    =       strtok_r (buffer, "\t", &saveptr) ;
		start  = strtoul (strtok_r (NULL,   "\t", &saveptr), NULL, 0);
		end    = strtoul (strtok_r (NULL,   "\t", &saveptr), NULL, 0);
		num_of_reads   = strtoul (strtok_r (NULL,   "\t", &saveptr), NULL, 0) ;
		num_of_origin  = strtoul (strtok_r (NULL,   "\t", &saveptr), NULL, 0) ;
		strand =       strtok_r (NULL,   "\t", &saveptr)[0];
		if (strand == '+')
		{
			umap[chr][strand][start] += (double(num_of_reads))/double(num_of_origin);
		}
		else
		{
			umap[chr][strand][end-1] += (double(num_of_reads))/double(num_of_origin);
		}
	}
}

int main (int argc, char* argv[])
{
	boost::program_options::options_description opts (R"(
This program calculate local ping pong score from bed2 format.
BED2 format modifies BED format in two ways: 
  Firstly, it replace the name field with number of times this read is been sequenced in the library. 
  Secondly, the score field has been replaced by the number of times this read been mapped to the genome. a.k.a, how many other loci, including itself, has the same sequence. 
Please contact bo.han@Umassmed.edu for any questions.
)" 	);
	/** options **/
	string file_a;
	string file_b;
	int num_of_threads;
	try {
			opts.add_options ()
				("help,h", "display this help message and exit")
				("file_a,a", boost::program_options::value<string>(&file_a)->required(), "first file used for ping-pong in BED2 format")
				("file_b,b", boost::program_options::value<string>(&file_b)->required(), "second file used for ping-pong in BED2 format")
				("thread,p", boost::program_options::value<int>(&num_of_threads)->default_value(1), "the number of threads to use")
			;
			boost::program_options::variables_map vm;
			boost::program_options::store (boost::program_options::parse_command_line(argc, argv, opts), vm);
			boost::program_options::notify (vm);
			if (vm.count("help"))	{ std::cerr << opts << std::endl; exit (1); }
		} catch (std::exception& e) {
			std::cerr << "Error: " << e.what() << std::endl;
			std::cerr << opts << std::endl;
			exit (1);
		} catch (...) {
			std::cerr << "Unknown error!" << std::endl;
			std::cerr << opts << std::endl;
			exit (1);
		} /** end of cmdline parsing **/

	unordered_map<string, unordered_map<char, unordered_map<uint64_t, double>>> seq2read1, seq2read2;

	if (num_of_threads > 1) {
		boost::thread parse_file_A (read_file_into_unordered_map, file_a, boost::ref (seq2read1));
		if ( file_a !=  file_b )
			read_file_into_unordered_map ( file_b , seq2read2);
		if ( parse_file_A.joinable () )
			parse_file_A.join ();
	} else {
		read_file_into_unordered_map ( file_a , seq2read1);
		if ( file_a !=  file_b )
			read_file_into_unordered_map ( file_b , seq2read2);
	}

	double  z_score [UPPERLIMIT];
	for (int i = 0; i < UPPERLIMIT; ++i)
		z_score[i]=0.0;
	if (num_of_threads == 1) {
		if ( file_a !=  file_b ) {
			for (int k = 0; k < UPPERLIMIT; ++k)
				do_ping_pong (seq2read1, seq2read2, k, z_score + k);
		} else {
			for (int k = 0; k < UPPERLIMIT; ++k)
				do_ping_pong (seq2read1, seq2read1, k, z_score + k);
		}
	} else { /** TODO: rewrite this part with threadpool **/
		int nthread = 0;
		int k = 0 ;
		boost::thread** threads = new boost::thread* [num_of_threads];
		while (k < UPPERLIMIT) {
			for ( nthread = 0 ;nthread < num_of_threads && k < UPPERLIMIT; ++nthread, ++k) {
				if ( file_a ==  file_b )
					threads[nthread] = new boost::thread ( do_ping_pong, boost::ref (seq2read1), boost::ref (seq2read1), k, z_score + k );
				else
					threads[nthread] = new boost::thread ( do_ping_pong, boost::ref (seq2read1), boost::ref (seq2read2), k, z_score + k );
			}
			for (int i = 0; i < nthread; ++ i ) {
				if (threads[i]->joinable ()) {
					threads[i]->join ();
					delete threads[i];
				}
			}
		}
		delete[] threads;
	}
	for (int k = 0; k < UPPERLIMIT; ++k)
	{
		cout <<  k+1 << '\t' << z_score[k] << '\n';
	}
}

void do_ping_pong (
		const unordered_map<string, unordered_map<char, unordered_map<uint64_t, double>>>& umapA,
		const unordered_map<string, unordered_map<char, unordered_map<uint64_t, double>>>& umapB,
		int k,
		double* const z_score
)
{
	for (auto const & chr_rest : umapA)
	{
		for (auto const & strand_rest : chr_rest.second)
		{
			for (auto const & position : strand_rest.second)
			{
				if (strand_rest.first == '+')
				{
					auto has_chr = umapB.find (chr_rest.first);
					if (has_chr != umapB.end ())
					{
						auto has_strand = has_chr->second.find ('-');
						if (has_strand != has_chr->second.end())
						{
							auto has_pos = has_strand->second.find (position.first+k);
							if (has_pos != has_strand->second.end())
							{
								*z_score += (position.second) * (has_pos->second);
							}
						}
						else
							continue;
					}
					else
						continue;
				}
				else
				{
					auto has_chr = umapB.find (chr_rest.first);
					if (has_chr != umapB.end ())
					{
						auto has_strand = has_chr->second.find ('+');
						if (has_strand != has_chr->second.end())
						{
							auto has_pos = has_strand->second.find (position.first-k);
							if (has_pos != has_strand->second.end())
								*z_score += (position.second) * (has_pos->second);
						}
						else
							continue;
					}
					else
						continue;
				}
			}
		}
	}
}
