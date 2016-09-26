cmake_minimum_required(VERSION 3.0)
project(piPipes_src)

find_package(Threads)
find_package(Boost COMPONENTS system iostreams program_options REQUIRED)
include_directories(${Boost_INCLUDE_DIRS})
find_package(ZLIB REQUIRED)
find_package(BZip2 REQUIRED)

add_executable(piPipes_fastq_to_insert piPipes_fastq_to_insert.cpp)
target_link_libraries(piPipes_fastq_to_insert
        boost_iostreams z bz2)

add_executable(piPipes_bed2_to_bedGraph
        ${piPipes_INCLUDEDIR}/piPipes_bed2_to_bedgraph.hpp
        piPipes_bed2_to_bedGraph.cpp)
target_link_libraries(piPipes_bed2_to_bedGraph
        boost_system boost_program_options ${CMAKE_THREAD_LIBS_INIT})

add_executable(piPipes_bed2Summary piPipes_bed2Summary.cpp)
target_link_libraries(piPipes_bed2Summary
        boost_system boost_program_options)

add_executable(piPipes_calculate_miRNA_heterogeneity
        ${piPipes_INCLUDEDIR}/piPipes_hairpin.hpp
        piPipes_hairpin.cpp
        piPipes_calculate_miRNA_heterogeneity.cpp)

add_executable(piPipes_filter_CIGAR
        ${piPipes_INCLUDEDIR}/piPipes_sam.hpp
        piPipes_sam.cpp
        piPipes_filter_CIGAR.cpp)
target_link_libraries(piPipes_filter_CIGAR
        boost_system boost_program_options)

add_executable(piPipes_insertBed_to_bed2 piPipes_insertBed_to_bed2.cpp)

add_executable(piPipes_local_ping_pong
        ${piPipes_INCLUDEDIR}/piPipes_ping_pong.hpp
        ${piPipes_INCLUDEDIR}/piPipes_common.hpp
        piPipes_common.cpp
        piPipes_ping_pong.cpp
        piPipes_local_ping_pong.cpp)
target_link_libraries(piPipes_local_ping_pong
        boost_system boost_program_options ${CMAKE_THREAD_LIBS_INIT})

add_executable(piPipes_phasing
        ${piPipes_INCLUDEDIR}/piPipes_ping_pong.hpp
        ${piPipes_INCLUDEDIR}/piPipes_common.hpp
        piPipes_common.cpp
        piPipes_ping_pong.cpp
        piPipes_phasing.cpp)
target_link_libraries(piPipes_phasing
        boost_system boost_program_options ${CMAKE_THREAD_LIBS_INIT})

set_target_properties(
        piPipes_fastq_to_insert
        piPipes_bed2_to_bedGraph
        piPipes_bed2Summary
        piPipes_calculate_miRNA_heterogeneity
        piPipes_insertBed_to_bed2
        piPipes_local_ping_pong
        piPipes_phasing
        PROPERTIES
        RUNTIME_OUTPUT_DIRECTORY ${piPipes_BINDIR}
)