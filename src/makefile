CC = gcc
CXX = g++
CXXFLAGS = -Wall -O3 -std=c++11 -static -I.
LIBS = -lboost_system -lboost_program_options
HTSLIB = $(HOME)/htslib/
PROG = piPipes_bed2Summary piPipes_bed2_to_bedGraph piPipes_calculate_miRNA_heterogeneity piPipes_fastq_to_insert piPipes_filter_CIGAR piPipes_insertBed_to_bed2 piPipes_local_ping_pong piPipes_phasing piPipes_bam_ZW_filter


all:	$(PROG)

piPipes_bed2Summary:	piPipes_bed2Summary.cpp
	$(CXX) $(CXXFLAGS) piPipes_bed2Summary.cpp $(LIBS) -o $@

piPipes_bed2_to_bedGraph:	piPipes_bed2_to_bedGraph.cpp
	$(CXX) $(CXXFLAGS) piPipes_bed2_to_bedGraph.cpp $(LIBS) -lboost_thread -pthread -lrt -o $@

piPipes_calculate_miRNA_heterogeneity:	piPipes_calculate_miRNA_heterogeneity.cpp piPipes_hairpin.hpp
	$(CXX) $(CXXFLAGS) piPipes_calculate_miRNA_heterogeneity.cpp $(LIBS) -o $@

piPipes_fastq_to_insert:	piPipes_fastq_to_insert.cpp gzstream.h 
	$(CXX) $(CXXFLAGS) piPipes_fastq_to_insert.cpp gzstream.C -lz $(LIBS) -o $@

piPipes_filter_CIGAR:	piPipes_filter_CIGAR.cpp piPipes_sam.hpp
	$(CXX) $(CXXFLAGS) piPipes_filter_CIGAR.cpp $(LIBS) -o $@

piPipes_insertBed_to_bed2:	piPipes_insertBed_to_bed2.cpp
	$(CXX) $(CXXFLAGS) piPipes_insertBed_to_bed2.cpp $(LIBS) -o $@

piPipes_local_ping_pong:	piPipes_local_ping_pong.cpp
	$(CXX) $(CXXFLAGS) piPipes_local_ping_pong.cpp -pthread -lrt $(LIBS) -o $@
piPipes_phasing:	piPipes_phasing.cpp
	$(CXX) $(CXXFLAGS) piPipes_phasing.cpp -pthread -lrt $(LIBS) -o $@

piPipes_bam_ZW_filter: piPipes_bam_ZW_filter.c
	$(CC) piPipes_bam_ZW_filter.c -lhts -pthread -lz -O3 -o $@ -I $(HTSLIB)  -L $(HTSLIB)
	
clean:
	rm -rf $(PROG)
