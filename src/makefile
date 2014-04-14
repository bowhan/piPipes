CXX = g++
CXXFLAGS = -Wall -O3 -std=c++11 -I.
LIBS = -lboost_system -lboost_program_options
PROG = piPipes_bed2Summary piPipes_bed2_to_bedGraph piPipes_calculate_miRNA_heterogeneity piPipes_fastq_to_insert piPipes_filter_CIGAR piPipes_insertBed_to_bed2 piPipes_local_ping_pong 

all:	$(PROG)

piPipes_bed2Summary:	piPipes_bed2Summary.cpp
	$(CXX) $(CXXFLAGS) piPipes_bed2Summary.cpp $(LIBS) -o $@

piPipes_bed2_to_bedGraph:	piPipes_bed2_to_bedGraph.cpp
	$(CXX) $(CXXFLAGS) piPipes_bed2_to_bedGraph.cpp $(LIBS) -lboost_thread -o $@

piPipes_calculate_miRNA_heterogeneity:	piPipes_calculate_miRNA_heterogeneity.cpp piPipes_hairpin.hpp
	$(CXX) $(CXXFLAGS) piPipes_calculate_miRNA_heterogeneity.cpp $(LIBS) -o $@

piPipes_fastq_to_insert:	piPipes_fastq_to_insert.cpp gzstream.h 
	$(CXX) $(CXXFLAGS) piPipes_fastq_to_insert.cpp gzstream.C -lz $(LIBS) -o $@

piPipes_filter_CIGAR:	piPipes_filter_CIGAR.cpp piPipes_sam.hpp
	$(CXX) $(CXXFLAGS) piPipes_filter_CIGAR.cpp $(LIBS) -o $@

piPipes_insertBed_to_bed2:	piPipes_insertBed_to_bed2.cpp
	$(CXX) $(CXXFLAGS) piPipes_insertBed_to_bed2.cpp $(LIBS) -o $@

piPipes_local_ping_pong:	piPipes_local_ping_pong.cpp
	$(CXX) $(CXXFLAGS) piPipes_local_ping_pong.cpp $(LIBS) -o $@

clean:
	rm -rf $(PROG)
