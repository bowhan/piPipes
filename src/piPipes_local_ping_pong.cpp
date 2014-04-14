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
#include <deque>
#include <unordered_map>
#include <thread>
#include <mutex>
#include <boost/program_options.hpp>
#include <boost/optional.hpp>

#define UPPERLIMIT 30
using namespace std;

template <typename T>
class ping_pong_table {
	deque<T> positions;
	mutex mx;
public:
	explicit ping_pong_table (T N): positions {} {
		for (T i {} ; i < N; ++i) {
			positions.push_back (i);
		}
	}
	boost::optional<T> pop () {
		unique_lock<mutex> locker (mx);
		boost::optional<T> ret {boost::none};
		if (!positions.empty ()) {
			ret = positions.back ();
			positions.pop_back ();
		}
		return ret;
	}
};

void do_ping_pong (
		const unordered_map<string, unordered_map<char, unordered_map<uint64_t, double>>>& ,
		const unordered_map<string, unordered_map<char, unordered_map<uint64_t, double>>>& ,
		int ,
		double* const
);

template <typename C, typename T>
class ping_pong_player {
private:
	C*  _A;
	C*  _B;
	ping_pong_table<T>* _tasks;
	double* _answers;
public:
	ping_pong_player<C,T> ( C& A,  C& B, ping_pong_table<T>& tasks, double* answers):
	_A {&A}, _B {&B}, _tasks {&tasks}, _answers {answers}
	{ }
	void operator () () {
		while (1) {
			auto n = _tasks->pop ();
			if (n) {
				do_ping_pong (*_A, *_B, *n, _answers + *n);
			} else {
				break;
			}
		}
	}
};

void read_file_into_unordered_map (const string& file_name, unordered_map<string, unordered_map<char, unordered_map<uint64_t, double>>>& umap) {
	ifstream in (file_name);
	if (!in) {
		cerr << "file " << file_name << " cannot be opened\n";
		exit (1);
	}
	char buffer [1024];
	char* saveptr;
	string chr;
	uint64_t num_of_reads, num_of_origin;
	uint64_t start, end;
	char strand;
	while (in.getline (buffer, 1024, '\n')) {
		char* saveptr;
		chr = strtok_r (buffer, "\t", &saveptr) ;
		start = strtoul (strtok_r (NULL,   "\t", &saveptr), NULL, 0);
		end = strtoul (strtok_r (NULL,   "\t", &saveptr), NULL, 0);
		num_of_reads = strtoul (strtok_r (NULL,   "\t", &saveptr), NULL, 0) ;
		num_of_origin = strtoul (strtok_r (NULL,   "\t", &saveptr), NULL, 0) ;
		strand = strtok_r (NULL,   "\t", &saveptr)[0];
		if (strand == '+') {
			umap[chr][strand][start] += (double(num_of_reads))/double(num_of_origin);
		} else {
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
	int upper_limit;
	try {
			opts.add_options ()
				("help,h", "display this help message and exit")
				("file_a,a", boost::program_options::value<string>(&file_a)->required(), "first file used for ping-pong in BED2 format")
				("file_b,b", boost::program_options::value<string>(&file_b)->required(), "second file used for ping-pong in BED2 format")
				("thread,p", boost::program_options::value<int>(&num_of_threads)->default_value(1), "the number of threads to use")
				("upper_limit,u", boost::program_options::value<int>(&upper_limit)->default_value(UPPERLIMIT), "the maximal overlap to calculate")
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
		thread parse_file_A (read_file_into_unordered_map, file_a, boost::ref (seq2read1));
		if ( file_a !=  file_b )
			read_file_into_unordered_map ( file_b , seq2read2);
		if ( parse_file_A.joinable () )
			parse_file_A.join ();
	} else {
		read_file_into_unordered_map ( file_a , seq2read1);
		if ( file_a !=  file_b )
			read_file_into_unordered_map ( file_b , seq2read2);
	}
	ping_pong_table<int> ppt {upper_limit};
	double* z_score = new double [upper_limit];
	for (int i = 0; i < upper_limit; ++i)
		z_score[i]=0.0;
	thread** threads = new thread* [num_of_threads];
	for (int i = 0; i < num_of_threads; ++i) {
		if ( file_a ==  file_b ) {
			threads[i] = new thread ( ping_pong_player<decltype(seq2read1), int> {seq2read1, seq2read1, ppt, z_score} );
		} else {
			threads[i] = new thread ( ping_pong_player<decltype(seq2read1), int> {seq2read1, seq2read2, ppt, z_score} );
		}
	}
	for (int i = 0; i < num_of_threads; ++ i ) {
		if (threads[i]->joinable ()) {
			threads[i]->join ();
			delete threads[i];
		}
	}
	delete[] threads;
	for (int k = 0; k < upper_limit; ++k) {
		cout <<  k+1 << '\t' << z_score[k] << '\n';
	}
	delete[] z_score;
}

void do_ping_pong (
		const unordered_map<string, unordered_map<char, unordered_map<uint64_t, double>>>& umapA,
		const unordered_map<string, unordered_map<char, unordered_map<uint64_t, double>>>& umapB,
		int k,
		double* const z_score )  {
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

