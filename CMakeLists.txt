cmake_minimum_required (VERSION 3.2)
project (PhyloCSF++ CXX)

find_package (OpenMP)

#include_directories ("include/libBigWig")
#link_directories ("include/libBigWig")

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O0 -g -march=native -Wall -pedantic -Wno-write-strings -Wno-pointer-arith -fopenmp") # -Wextra -DNOCURL -static -fsanitize=address -fno-omit-frame-pointer

#set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS}")

add_executable (phylocsf++ phylocsfpp.cpp)
add_executable (memtest phylo_memtest.cpp)
add_executable (tests automatic_tests.cpp)
#add_executable (reader parallel_file_reader.hpp)

target_link_libraries (phylocsf++ gsl)
target_link_libraries (phylocsf++ gslcblas)

target_link_libraries (memtest gsl)
target_link_libraries (memtest gslcblas)

target_link_libraries (tests gsl)
target_link_libraries (tests gslcblas)

#target_link_libraries (phylocsf++ BigWig)
#target_link_libraries (phylocsf++ z) # needed by libbigwig
#target_link_libraries (phylocsf++ pthread)
