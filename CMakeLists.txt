cmake_minimum_required(VERSION 2.8)
project(phylocsf)

set (CMAKE_BUILD_TYPE Debug)

set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O0 -g -fsanitize=address -fno-omit-frame-pointer") #  -static-libasan  -Wall -pedantic

#set (CMAKE_STATIC_LINKER_FLAGS "${CMAKE_STATIC_LINKER_FLAGS} -fsanitize=address -static-libasan")

# find_package (OpenMP REQUIRED)
# if (OPENMP_FOUND)
#    set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
#    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
#    set (CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
#endif()

add_executable(phylocsf prototype.cpp src/newick.hpp src/ecm.hpp)
target_link_libraries (phylocsf gsl gslcblas) # pthread
