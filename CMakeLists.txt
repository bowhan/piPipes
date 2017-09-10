project(piPipes C CXX)
cmake_minimum_required(VERSION 3.0)

foreach (level MAJOR MINOR PATCH)
    file(STRINGS piPipes
            _define_${level}
            REGEX "PROG_${level}_VERSION=")
    string(REGEX MATCH "([0-9]+)" _${level} "${_define_${level}}")
endforeach ()

set(piPipes_VERSION "${_MAJOR}.${_MINOR}.${_PATCH}")
set(piPipes_ROOTDIR ${PROJECT_SOURCE_DIR})
set(piPipes_INCLUDEDIR ${piPipes_ROOTDIR}/include)
set(piPipes_SOURCEDIR ${piPipes_ROOTDIR}/src)
set(piPipes_LIBDIR ${piPipes_ROOTDIR}/lib)
set(piPipes_BINDIR ${piPipes_ROOTDIR}/bin)

SET(CMAKE_INSTALL_RPATH ${CMAKE_INSTALL_PREFIX}/lib)
SET(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)

find_package(Threads)
find_package(Boost COMPONENTS system iostreams program_options REQUIRED)
find_package(ZLIB REQUIRED)
find_package(BZip2 REQUIRED)

#include_directories(${piPipes_INCLUDEDIR})
#file(MAKE_DIRECTORY ${piPipes_LIBDIR})
set(CMAKE_CXX_FLAGS "-g -std=c++11 -Wall")
#set(CMAKE_EXE_LINKER_FLAGS "-static")
add_subdirectory(${piPipes_SOURCEDIR})
