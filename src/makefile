CXX = g++
CXXFLAGS = -Wall -O3 -std=c++11 -I.
LIBS = -lz -pthread -lboost_system -lboost_program_options
PROG = piPipes_bed2Summary piPipes_bed2_to_bedGraph piPipes_calculate_miRNA_heterogeneity piPipes_fastq_to_insert piPipes_filter_CIGAR piPipes_insertBed_to_bed2 piPipes_local_ping_pong 

all:	$(PROG)

piPipes_bed2Summary:	piPipes_bed2Summary.cpp
	$(CC) $(CXXFLAGS) piPipes_bed2Summary.cpp $(LIBS) -o $@

piPipes_bed2_to_bedGraph:	piPipes_bed2_to_bedGraph.cpp
	$(CC) $(CXXFLAGS) piPipes_bed2_to_bedGraph.cpp $(LIBS) -lboost_thread -o $@

piPipes_calculate_miRNA_heterogeneity:	piPipes_calculate_miRNA_heterogeneity.cpp piPipes_hairpin.hpp
	$(CC) $(CXXFLAGS) piPipes_calculate_miRNA_heterogeneity.cpp $(LIBS) -o $@

piPipes_fastq_to_insert:	piPipes_fastq_to_insert.cpp gzstream.h 
	$(CC) $(CXXFLAGS) piPipes_fastq_to_insert.cpp gzstream.C $(LIBS) -o $@

piPipes_filter_CIGAR:	piPipes_filter_CIGAR.cpp piPipes_sam.hpp
	$(CC) $(CXXFLAGS) piPipes_filter_CIGAR.cpp $(LIBS) -o $@

piPipes_insertBed_to_bed2:	piPipes_insertBed_to_bed2.cpp
	$(CC) $(CXXFLAGS) piPipes_insertBed_to_bed2.cpp $(LIBS) -o $@

piPipes_local_ping_pong:	piPipes_local_ping_pong.cpp
	$(CC) $(CXXFLAGS) piPipes_local_ping_pong.cpp $(LIBS) -o $@

clean:
	rm -rf $(PROG)
