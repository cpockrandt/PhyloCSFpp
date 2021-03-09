cmake_minimum_required(VERSION 3.2)
string(ASCII 27 Esc)
set(ColourBold  "${Esc}[1m")
set(ColourReset "${Esc}[m")
set(ColourRed   "${Esc}[31m")

include(ExternalProject)

message("${ColourBold}Compiler Detection${ColourReset}") # This needs to go before "project ()".

project(PhyloCSF++ CXX)

# ----------------------------------------------------------------------------
# Compiler detection
# ----------------------------------------------------------------------------
# TODO: check whether older compilers work as well
if (CMAKE_CXX_COMPILER_ID MATCHES "GNU")
    if (CMAKE_CXX_COMPILER_VERSION VERSION_LESS 4.9.1)
        message(FATAL_ERROR "Your GCC version is too old. Minimum version is GCC-4.9.1!")
        return()
    endif ()
elseif (CMAKE_CXX_COMPILER_ID MATCHES "Clang")
    if (CMAKE_CXX_COMPILER_VERSION VERSION_LESS 3.6)
        message(FATAL_ERROR "Your Clang version seems too old. Please upgrade to 3.8.0 or use GCC.")
        return()
    elseif (CMAKE_CXX_COMPILER_VERSION VERSION_LESS 3.8)
        message(WARNING "Your Clang version is too old, you will not have parallelism! Upgrade to 3.8.0 or newer.")
    endif ()
elseif (CMAKE_CXX_COMPILER_ID MATCHES "Intel")
    if (CMAKE_CXX_COMPILER_VERSION VERSION_LESS 17)
        message(WARNING "Your Intel Compiler version is too old. Please upgrade to 17.0.0 or newer!")
    endif ()
else ()
    message(WARNING "Unknown compiler, you are on your own!")
endif ()

message(STATUS "Compiling with ${CMAKE_CXX_COMPILER_ID} ${CMAKE_CXX_COMPILER_VERSION}.")
message(STATUS "The requirements were met.")

# ----------------------------------------------------------------------------
# Dependency detection
# ----------------------------------------------------------------------------

message("\n${ColourBold}Dependency detection${ColourReset}")

find_package(OpenMP QUIET)
find_package(Git)
find_package(GSL REQUIRED)

option(PHYLOCSFPP_BUILD_LIBBIGWIG "Turn on/off building step of libBigWig (used for bioconda)." ON)
if (NOT PHYLOCSFPP_BUILD_LIBBIGWIG)
    add_definitions(-DPHYLOCSFPP_BUILD_LIBBIGWIG=0)
    set(LIBBIGWIG_MAKE_CMD "")
else ()
    add_definitions(-DPHYLOCSFPP_BUILD_LIBBIGWIG=1)
    set(LIBBIGWIG_MAKE_CMD "make")
endif ()

# TODO: separate source and binary dir?
ExternalProject_Add(libBigWig
    GIT_REPOSITORY "https://github.com/cpockrandt/libBigWig"
    UPDATE_COMMAND ""
    DOWNLOAD_DIR "${CMAKE_SOURCE_DIR}/external/"
    SOURCE_DIR "${CMAKE_BINARY_DIR}/libBigWig"
    BUILD_IN_SOURCE 1
    CONFIGURE_COMMAND ""
    BUILD_COMMAND "${LIBBIGWIG_MAKE_CMD}"
    INSTALL_COMMAND ""
)

# Pass git hash and git date to C++ compiler as macros to include them into the help string.
# Source: https://stackoverflow.com/a/21028226 (by Drew Noakes)
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

message(STATUS "These dependencies where found:")
message(   "     GSL        ${GSL_FOUND}      ${GSL_VERSION}")
message(   "     OPENMP     ${OPENMP_FOUND}      ${OpenMP_CXX_FLAGS}")
message(   "     GIT        ${GIT_FOUND}      (optional)")

if (OPENMP_FOUND)
    add_definitions(-DENABLE_OPENMP)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
else ()
    message(WARNING "WARNING WARNING WARNING\nWARNING: OpenMP not found. PhyloCSF++ will be built without multi-threading! "
            "This is probably not what you want! Use GCC >= 4.9.1 or Clang >= 3.8.0\nWARNING WARNING WARNING")
endif ()

include_directories(${GSL_INCLUDE_DIRS})

message(STATUS "The requirements were met.")

# ----------------------------------------------------------------------------
# Build configuration
# ----------------------------------------------------------------------------

message("\n${ColourBold}Build configuration${ColourReset}")

# ----------------------------------------------------------------------------
# Make "Release" the default cmake build type
# ----------------------------------------------------------------------------

if (NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release CACHE STRING
            "Choose the type of build, options are: Debug Release RelWithDebInfo"
            FORCE)
endif ()

# ----------------------------------------------------------------------------
# Warn if cmake build type is not "Release"
# ----------------------------------------------------------------------------

if (NOT CMAKE_BUILD_TYPE STREQUAL Release)
    message(STATUS "${ColourRed}CMAKE_BUILD_TYPE is not \"Release\", your binaries will be slow.${ColourReset}")
endif ()

option(PHYLOCSFPP_NATIVE_BUILD "Architecture-specific optimizations, i.e., g++ -march=native." OFF)
if (PHYLOCSFPP_NATIVE_BUILD)
    add_definitions(-DPHYLOCSFPP_NATIVE_BUILD=1)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -march=native")
    if (CMAKE_CXX_COMPILER_ID MATCHES "Intel")
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -xHOST -ipo -no-prec-div -fp-model fast=2")
    endif ()
endif ()

option(PHYLOCSFPP_STATIC_BUILD "Static build." OFF)
if (PHYLOCSFPP_STATIC_BUILD)
    add_definitions(-DPHYLOCSFPP_STATIC_BUILD=1)
    set(CMAKE_FIND_LIBRARY_SUFFIXES ".a")
    set(BUILD_SHARED_LIBS OFF)
    set (CMAKE_EXE_LINKER_FLAGS "-static -static-libgcc -static-libstdc++")
    # apple does not support fully static builds, but at least libgcc and libstdc++
    if (APPLE)
        message(WARNING "WARNING: Builds on Mac are never fully static.")
    endif ()
    # on linux cmake adds -rdynamic automatically which clang can't handle in static builds
    # if (CMAKE_SYSTEM_NAME MATCHES "Linux")
    #     set(CMAKE_SHARED_LIBRARY_LINK_CXX_FLAGS "")
    # endif ()
endif ()

include_directories(SYSTEM ${CMAKE_BINARY_DIR}/libBigWig)
link_directories(${CMAKE_BINARY_DIR}/libBigWig)

message(STATUS "The following options are selected for the build:")
message(   "     PHYLOCSFPP_NATIVE_BUILD             ${PHYLOCSFPP_NATIVE_BUILD}")
message(   "     PHYLOCSFPP_STATIC_BUILD             ${PHYLOCSFPP_STATIC_BUILD}")
message(   "     PHYLOCSFPP_BUILD_LIBBIGWIG          ${PHYLOCSFPP_BUILD_LIBBIGWIG}")
message(STATUS "Run 'cmake -LH' to get a comment on each option.")
message(STATUS "Remove CMakeCache.txt and re-run cmake with -DOPTIONNAME=ON|OFF to change an option.")

# ----------------------------------------------------------------------------
# Build Setup
# ----------------------------------------------------------------------------

add_definitions(-DCMAKE_BUILD_TYPE="${CMAKE_BUILD_TYPE}")

# TODO: -Wconversion
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wextra -pedantic -DNOCURL -std=c++11 -fpermissive -DHAVE_INLINE")

add_executable(phylocsf++ src/phylocsf++.cpp)
add_dependencies(phylocsf++ libBigWig)
target_link_libraries(phylocsf++ z)
target_link_libraries(phylocsf++ gsl)
target_link_libraries(phylocsf++ gslcblas)
target_link_libraries(phylocsf++ BigWig)

install(TARGETS phylocsf++ DESTINATION bin)

message("\n${ColourBold}Finishing up${ColourReset}")
