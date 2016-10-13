cmake_minimum_required(VERSION 3.0)
project(piPipes_misc_tools)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")

set(MISC_INCLUDE_DIR ${PROJECT_SOURCE_DIR}/include)
set(MISC_BIN_DIR ${PROJECT_SOURCE_DIR}/bin)

include_directories(${MISC_INCLUDE_DIR})
add_subdirectory(src)
