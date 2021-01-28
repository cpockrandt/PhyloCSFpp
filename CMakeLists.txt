cmake_minimum_required (VERSION 3.2)
project (PhyloCSF++ CXX)

find_package (OpenMP)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -march=native -Wall -pedantic -Wno-write-strings -Wno-pointer-arith -fopenmp") # -Wextra -DNOCURL -static -fsanitize=address -fno-omit-frame-pointer

add_executable (phylocsf++ src/phylocsf++.cpp)
target_link_libraries (phylocsf++ gsl)
target_link_libraries (phylocsf++ gslcblas)
#target_link_libraries (phylocsf++ pthread)

#add_executable (memtest phylo_memtest.cpp)
#target_link_libraries (memtest gsl)
#target_link_libraries (memtest gslcblas)

#add_executable (tests automatic_tests.cpp)
#target_link_libraries (tests gsl)
#target_link_libraries (tests gslcblas)
