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

#ifndef bowhan_chiSquared_h
#define bowhan_chiSquared_h

#include <iostream>
#include <fstream>
#include <string>
#include <memory>
#include <vector>
#include <cmath>
#include <iomanip>
#include "boost/algorithm/string.hpp"
#include "boost/lexical_cast.hpp"

template <typename T>
class chiSquared {
private:
    uint32_t _N = 0;
    std::vector<T> *_data = nullptr;
    T *_count = nullptr;
    std::vector<std::string> _names;
public:
    chiSquared(uint32_t N)
        :
        _N(N), _data{new std::vector<T>[N]}, _count{new T[N]} {
        std::cerr << "data initialization" << std::endl;
    }

    chiSquared(const chiSquared&) = delete;

    chiSquared& operator=(const chiSquared&) = delete;

    chiSquared(chiSquared&& other)
        :
        _N{other._N}, _data{other._data}, _count{other._count}, _names{other._names} {
        other._data = nullptr;
        other._count = nullptr;
        other.~chiSquared();
    }

    chiSquared& operator=(chiSquared&& other) {
        _N = other._N;
        _data = other._data;
        _count = other._count;
        _names = other._names;
        other._data = nullptr;
        other._count = nullptr;
        other.~chiSquared();
    }

    ~chiSquared() {
        delete[] _data;
        delete[] _count;
    }

    void readFile(std::ifstream *in) {
        std::cout << std::setprecision(2) << std::fixed;
        std::string line;
        std::vector<std::string> tokens(_N + 1);
        T tempValue{};
        T tableTotal{};
        while (getline(*in, line)) {
            boost::split(tokens, line, boost::is_any_of("\t"));
            _names.push_back(tokens[0]);
            for (uint32_t i = 1; i <= _N; ++i) {
                try {
                    tempValue = boost::lexical_cast<T>(tokens[i]);
                } catch (const boost::bad_lexical_cast& e) {
                    tempValue = T {};
                    std::cerr << "Warning: illegal character " << tokens[i] << std::endl;
                }
                _data[i - 1].push_back(tempValue);
                _count[i - 1] += tempValue;
                tableTotal += tempValue;
            }
        }
        std::cerr << "finishing reading file" << std::endl;
        T results{};
        T expected{};
        T rowTotal{};
        T yets{};
        std::vector<T> rowTwo(_N);
        for (size_t i = 0; i < _names.size(); ++i) {
            rowTotal = T{};
            yets = 0.0; /// initialize yets
            /** first iterator to get row 2 and total */
            for (uint32_t j = 0; j < _N; ++j) {
                rowTotal += _data[j][i];
                if (_data[j][i] < 5.0)
                    yets = 0.5; /// if any cell is smaller than 5, use yets coorection
            }
            results = T {};
            std::cout << _names[i];
            for (uint32_t j = 0; j < _N; ++j) {
                expected = rowTotal * _count[j] / tableTotal;
                std::cout << '\t' << _data[j][i] << '\t' << expected;
                results +=
                    std::pow((_data[j][i] - expected > 0.0 ? (_data[j][i] - expected) : (expected - _data[j][i])) - yets
                             , 2) / expected;

                expected = (tableTotal - rowTotal) * _count[j] / tableTotal;
                results += std::pow(
                    (_count[j] - _data[j][i] - expected > 0 ? (_count[j] - _data[j][i] - expected) : (_data[j][i]
                        + expected - _count[j])) - yets, 2) / expected;
            }
            std::cout << '\t' << results << '\n';
        }
    }
};


#endif