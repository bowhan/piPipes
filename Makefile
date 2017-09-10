all:
	mkdir build && cd build && cmake -DCMAKE_INSTALL_PREFIX=.. .. && make && make install

clean:
	rm -r build 

.PHONY: clean