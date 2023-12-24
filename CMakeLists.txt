cmake_minimum_required(VERSION 3.2)

project(qlibs_cpp_test)

add_subdirectory( src )
include_directories( src )

set(CMAKE_VERBOSE_MAKEFILE ON)
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/build/lib)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/build/lib)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/build/bin)

add_executable ( qlibs_cpp_test check/qlibs_cpp_test.cpp )
target_include_directories( qlibs_cpp_test PUBLIC ${PROJECT_SOURCE_DIR}/src )
target_link_libraries( qlibs_cpp_test qlibs-cpp )