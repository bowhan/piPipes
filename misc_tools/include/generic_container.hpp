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

#ifndef bowhan_genericcontainer_h
#define bowhan_genericcontainer_h

#include <deque>
#include <queue>
#include <vector>
#include <stack>

/* policies for generic containers */
template <typename T, template <typename...> class C>
struct containerPolicy {};

/* Vector */
template <typename T>
struct containerPolicy<T, std::vector> {

    template <typename U>
    static void push(std::vector<T>& vec, U&& data) {
        vec.push_back(data);
    }

    static void pop(std::vector<T>& vec) {
        vec.pop_back();
    }

    static T& top(std::vector<T>& vec) {
        return vec.back();
    }
};

/* Deque */
template <typename T>
struct containerPolicy<T, std::deque> {
    template <typename U>
    static void push(std::deque<T>& deq, U&& data) {
        deq.push_back(data);
    }

    static void pop(std::deque<T>& deq) {
        deq.pop_back();
    }

    static T& top(std::deque<T>& deq) {
        return deq.back();
    }
};

/* Queue */
template <typename T>
struct containerPolicy<T, std::queue> {
    template <typename U>
    static void push(std::queue<T>& que, U&& data) {
        que.push(data);
    }

    static void pop(std::queue<T>& que) {
        que.pop();
    }

    static T& top(std::queue<T>& que) {
        return que.front();
    }

};

/* Stack */
template <typename T>
struct containerPolicy<T, std::stack> {
    template <typename U>
    static void push(std::stack<T>& stk, U&& data) {
        stk.push(data);
    }

    static void pop(std::stack<T>& stk) {
        stk.pop();
    }

    static T& top(std::stack<T>& stk) {
        return stk.top();
    }
};

#endif