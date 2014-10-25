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

#ifndef bowhan_insertMerge_h
#define bowhan_insertMerge_h

#include <string>
#include <fstream>
#include <vector>
#include <unordered_map>
#include "boost/algorithm/string.hpp"

void readFile(const std::string& fileName, std::unordered_map<std::string,std::string>* storage, int commonField, int mergeField) {
	std::ifstream in1(fileName);
	std::string line{};
	std::vector<std::string> tokens;
	if (mergeField != -1) {
		while(getline(in1, line)) {
			boost::split(tokens, line, boost::is_any_of("\t"));
			assert(tokens.size() > commonField && tokens.size() > mergeField);
			storage->insert(std::pair<std::string, std::string> {tokens[commonField], tokens[mergeField]});
		}
	} else {
		while(getline(in1,line)) {
			boost::split(tokens, line, boost::is_any_of("\t"));
			assert(tokens.size() > commonField && tokens.size() > mergeField);
			storage->insert(std::pair<std::string, std::string> {tokens[commonField], line});
		}
	}
}

#endif
