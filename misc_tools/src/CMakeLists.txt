find_package(Boost COMPONENTS system program_options iostreams REQUIRED)
find_package(ZLIB REQUIRED)
find_package(BZip2 REQUIRED)
find_package(Threads)

include_directories(${BOOST_INCLUDE_DIRS})

add_executable(chiSquared
        ${MISC_INCLUDE_DIR}/chiSquared.hpp
        chiSquared.cpp
        )

add_executable(insertMerge
        ${MISC_INCLUDE_DIR}/insertMerge.hpp
        insertMerge.cpp
        )
target_link_libraries(insertMerge
        ${Boost_LIBRARIES}
        ${CMAKE_THREAD_LIBS_INIT}
        )

add_executable(kmer_mean_len_from_insert
        kmer_mean_len_from_insert.cpp
        )

add_executable(splitBarcode
        splitBarcode.cpp
        )
target_link_libraries(splitBarcode
        ${Boost_LIBRARIES}
        ${ZLIB_LIBRARIES}
        ${BZip2_LIBRARIES}
        ${CMAKE_THREAD_LIBS_INIT}
        )

add_executable(pp_detail
        pp_detail.cpp
        )
target_link_libraries(pp_detail
        ${Boost_LIBRARIES}
        )

set_target_properties(chiSquared insertMerge kmer_mean_len_from_insert splitBarcode pp_detail
        PROPERTIES
        RUNTIME_OUTPUT_DIRECTORY ${MISC_BIN_DIR}
        )