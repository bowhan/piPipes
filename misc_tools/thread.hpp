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

#ifndef bowhan_thread_h
#define bowhan_thread_h

#include <thread>
#include <condition_variable>
#include "boost/optional.hpp"
#include "generic_container.hpp"

/* implementation of thread safe queue */
template<typename T, template <typename, typename...> class C>
class threadSafeQueue {
public:
	threadSafeQueue<T,C> () { }
	
	bool empty () const {
		return _quene.empty ();
	}
	
	bool done () const {
		return _finished;
	}
	
	void setDone() noexcept {
		std::unique_lock<std::mutex> lock(_mx);
		_finished = true;
		_cond.notify_all();
	}
	
	void push(T&& data) {
		std::unique_lock<std::mutex> lock(_mx);
		bool const was_empty = _quene.empty ();
		containerPolicy<T,C>::push(_quene, std::forward(data));
		lock.unlock();
		if (was_empty)
			_cond.notify_one();
	}
	
	void push(char* data) {
		std::unique_lock<std::mutex> lock(_mx);
		bool const was_empty = _quene.empty ();
		containerPolicy<T,C>::push(_quene, T{data});
		lock.unlock();
		if (was_empty){
			_cond.notify_one();
		}
	}
	
	bool try_pop (T& data) {
		std::unique_lock<std::mutex> lock(_mx);
		_cond.wait(lock, [this]()->bool { return !_quene.empty() || _finished; });
		if (_quene.empty()) return false;
		else {
			data = containerPolicy<T,C>::top(_quene);
			containerPolicy<T,C>::pop(_quene);
			return true;
		}
	}
private:
	C<T>  _quene;
	bool _finished = false;
	mutable std::mutex _mx;
	std::condition_variable _cond;
	threadSafeQueue<T,C> (const threadSafeQueue<T,C>&) = delete;
	threadSafeQueue<T,C> (threadSafeQueue<T,C>&&) = delete;
	threadSafeQueue<T,C> operator=(const threadSafeQueue<T,C>&) = delete;
	threadSafeQueue<T,C> operator=(threadSafeQueue<T,C>&&) = delete;
};

#endif