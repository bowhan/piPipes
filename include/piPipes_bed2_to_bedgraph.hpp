#ifndef PIPIPES_PIPIPES_BED2_TO_BEDGRAPH_HPP
#define PIPIPES_PIPIPES_BED2_TO_BEDGRAPH_HPP

#include <iostream>
#include <string>
#include <unordered_map>
#include <boost/algorithm/string.hpp>

template <typename T, char S>
class bed2_bedgraph_converter {
private:
    std::string infile_;
    std::ostream *outtream_;
    unsigned long max_size_;
    std::unordered_map<std::string, unsigned long> *sz_table_;
public:
    bed2_bedgraph_converter(std::string in_file, std::ostream *out, unsigned long size
                            , std::unordered_map<std::string, unsigned long> *size_table
                           )
        :
        infile_{in_file}, outtream_{out}, max_size_{size}, sz_table_{size_table} {
    }

    void operator()();
};

template <typename T, char S>
void
bed2_bedgraph_converter<T, S>::operator()() {
    using namespace std;
    *outtream_ << std::setprecision(5);
    *outtream_ << std::fixed;
    ifstream *_in = new ifstream{infile_};
    T *_container = new T[max_size_];
    T *p = _container;
    T previous_value{};
    unsigned long start_point{};
    for (unsigned long i = 0; i < max_size_; ++i)
        *p = T{};
    string line;
    vector<string> fields;
    string chr;
    unsigned long pos;
    T num_of_reads;
    T num_of_mapping;
    getline(*_in, line);
    boost::split(fields, line, boost::is_any_of("\t"));
    chr = fields[0];
    _in->seekg(ios_base::beg); /// C++11 says it also clear eofbit flag
    while (getline(*_in, line)) {
        boost::split(fields, line, boost::is_any_of("\t"));
        if (fields[5][0] == S) { /// if this is the strand to do
            if (!chr.empty() && chr != fields[0]) { /// new chromosome
                /* begin to output */
                p = _container;
                previous_value = *p++;
                start_point = 0;
                unsigned long i{1};
                for (; i < (*sz_table_)[chr]; ++i) {
                    if (*p != previous_value) { /// new value, need to output previous value
                        if (previous_value) {
                            if (S == '+') /// compiling time decision
                                *outtream_ << chr << '\t' << start_point << '\t' << i << '\t' << previous_value
                                           << '\n';
                            if (S == '-')
                                *outtream_ << chr << '\t' << start_point << '\t' << i << "\t-" << previous_value
                                           << '\n';
                        }
                        start_point = i;
                        previous_value = *p;
                    }
                    *p++ = T{};
                }
                if (previous_value) {
                    if (S == '+') /// compiling time decision
                        *outtream_ << chr << '\t' << start_point << '\t' << i << '\t' << previous_value << '\n';
                    if (S == '-')
                        *outtream_ << chr << '\t' << start_point << '\t' << i << "\t-" << previous_value << '\n';
                }
                /* end of output */
                for (i = 0, p = _container; i < (*sz_table_)[fields[0]]; ++i) /// continue to blank container
                    *p++ = T{};
                chr = fields[0];
            }
            num_of_reads = boost::lexical_cast<T>(fields[3]);
            num_of_mapping = boost::lexical_cast<T>(fields[4]);
            if (fields[5] == "+")
                pos = stoul(fields[1]);
            else if (fields[5] == "-")
                pos = stoul(fields[2]) - 1;
            else {
                cerr << "Error: unknown strand " << fields[5] << endl;
                exit(1);
            }
            if (pos < (*sz_table_)[chr])
                _container[pos] += num_of_reads / num_of_mapping;
            else
                cerr << "Warning: chr size overflowing \n" << line << endl;
        }
    }
    p = _container;
    previous_value = *p++;
    start_point = 0;
    unsigned long i{1};
    for (; i < (*sz_table_)[chr]; ++i, ++p) {
        if (*p != previous_value) { /// new value, need to output previous value
            if (previous_value) {
                if (S == '+') /// compiling time decision
                    *outtream_ << chr << '\t' << start_point << '\t' << i << '\t' << previous_value << '\n';
                if (S == '-')
                    *outtream_ << chr << '\t' << start_point << '\t' << i << "\t-" << previous_value << '\n';
            }
            start_point = i;
            previous_value = *p;
        }
        //*p++ = T{}; no need to re-initialize
    }
    if (previous_value) {
        if (S == '+') /// compiling time decision
            *outtream_ << chr << '\t' << start_point << '\t' << i << '\t' << previous_value << '\n';
        if (S == '-')
            *outtream_ << chr << '\t' << start_point << '\t' << i << "\t-" << previous_value << '\n';
    }
    delete[] _container;
    _in->close();
    delete _in;
}


#endif //PIPIPES_PIPIPES_BED2_TO_BEDGRAPH_HPP
