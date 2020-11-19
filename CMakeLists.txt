cmake_minimum_required(VERSION 3.14)
project(phylocsf)

set (CXX11_FOUND "1")
set (CMAKE_PREFIX_PATH "/home/chris/dev-uni/seqan/util/cmake")
set (SEQAN_INCLUDE_PATH "/home/chris/dev-uni/seqan/include")

set (CMAKE_BUILD_TYPE Debug)

# set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -pedantic")

find_package (SeqAn REQUIRED)
# find_package (OpenMP REQUIRED)
# if (OPENMP_FOUND)
#    set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
#    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
#    set (CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
#endif()

include_directories (${SEQAN_INCLUDE_DIRS})
add_definitions (${SEQAN_DEFINITIONS})

set(CMAKE_CXX_STANDARD 14)

add_executable(phylocsf prototype.cpp src/newick.hpp src/ecm.hpp)
target_link_libraries (phylocsf ${SEQAN_LIBRARIES})
target_link_libraries (phylocsf gsl gslcblas) #pthread
