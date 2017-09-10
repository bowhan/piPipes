include(files.cmake)

add_library(piPipesLib SHARED ${HEADER_FILES} ${SOURCE_FILES})
target_include_directories(piPipesLib PUBLIC ${piPipes_INCLUDEDIR})

target_link_libraries(piPipesLib
        ${ZLIB_LIBRARIES}
        ${Bzlib2_LIBRARIES}
        ${CMAKE_THREAD_LIBS_INIT}
        ${Boost_LIBRARIES}
        )



add_executable(piPipes_fastq_to_insert piPipes_fastq_to_insert.cpp)
target_link_libraries(piPipes_fastq_to_insert PRIVATE piPipesLib)

add_executable(piPipes_bed2_to_bedGraph piPipes_bed2_to_bedGraph.cpp)
target_link_libraries(piPipes_bed2_to_bedGraph PRIVATE piPipesLib)

add_executable(piPipes_bed2Summary piPipes_bed2Summary.cpp)
target_link_libraries(piPipes_bed2Summary PRIVATE piPipesLib)

add_executable(piPipes_calculate_miRNA_heterogeneity piPipes_calculate_miRNA_heterogeneity.cpp)
target_link_libraries(piPipes_calculate_miRNA_heterogeneity PRIVATE piPipesLib)

add_executable(piPipes_filter_CIGAR piPipes_filter_CIGAR.cpp)
target_link_libraries(piPipes_filter_CIGAR PRIVATE piPipesLib)


add_executable(piPipes_insertBed_to_bed2 piPipes_insertBed_to_bed2.cpp)

add_executable(piPipes_local_ping_pong piPipes_local_ping_pong.cpp)
target_link_libraries(piPipes_local_ping_pong PRIVATE piPipesLib)

add_executable(piPipes_phasing piPipes_phasing.cpp)
target_link_libraries(piPipes_phasing PRIVATE piPipesLib)

install(TARGETS
        piPipesLib
        piPipes_fastq_to_insert
        piPipes_bed2_to_bedGraph
        piPipes_bed2Summary
        piPipes_calculate_miRNA_heterogeneity
        piPipes_insertBed_to_bed2
        piPipes_local_ping_pong
        piPipes_phasing
        PUBLIC_HEADER DESTINATION include
        ARCHIVE DESTINATION lib
        LIBRARY DESTINATION lib
        RUNTIME DESTINATION bin
        )
