cmake_minimum_required(VERSION 2.8.12)
project(PhyloCSF++ CXX)

include(ExternalProject)

find_package(OpenMP REQUIRED)

if (OPENMP_FOUND)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
endif()

# Pass git hash and git date to C++ compiler as macros to include them into the help string.
# Source: https://stackoverflow.com/a/21028226 (by Drew Noakes)
find_package(Git)

if (GIT_FOUND)
    execute_process(COMMAND
        "${GIT_EXECUTABLE}" describe --match=NeVeRmAtCh --always --abbrev=8 --dirty
        WORKING_DIRECTORY "${CMAKE_SOURCE_DIR}"
        OUTPUT_VARIABLE GIT_HASH
        ERROR_QUIET OUTPUT_STRIP_TRAILING_WHITESPACE)
    execute_process(COMMAND
        "${GIT_EXECUTABLE}" log -1 --date=short --format=%cd
        WORKING_DIRECTORY "${CMAKE_SOURCE_DIR}"
        OUTPUT_VARIABLE GIT_DATE
        ERROR_QUIET OUTPUT_STRIP_TRAILING_WHITESPACE)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DGIT_HASH=\"${GIT_HASH}\" -DGIT_DATE=\"${GIT_DATE}\"")
endif()

ExternalProject_Add(libBigWig
    GIT_REPOSITORY "https://github.com/cpockrandt/libBigWig"
    UPDATE_COMMAND ""
    DOWNLOAD_DIR "${CMAKE_CURRENT_SOURCE_DIR}/external/"
    SOURCE_DIR "${CMAKE_CURRENT_SOURCE_DIR}/external/libBigWig"
    BUILD_IN_SOURCE 1
    CONFIGURE_COMMAND ""
    BUILD_COMMAND ${MAKE}
    INSTALL_COMMAND ""
)

include_directories(${CMAKE_CURRENT_SOURCE_DIR}/external/libBigWig)
link_directories(${CMAKE_CURRENT_SOURCE_DIR}/external/libBigWig)

# -static -fsanitize=address -fno-omit-frame-pointer

# TODO: fix errors
# -Werror
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wextra -pedantic -DNOCURL -std=c++11 -fpermissive") # -Wconversion

add_executable(phylocsf++ src/phylocsf++.cpp)
add_dependencies(phylocsf++ libBigWig)
target_link_libraries(phylocsf++ gsl gslcblas BigWig z)

#add_executable(tests automatic_tests.cpp)
#target_link_libraries(tests gsl gslcblas)
