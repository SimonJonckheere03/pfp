cmake_minimum_required(VERSION 3.5.0)

project(PFP++)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

################################################################################
# Set a default build type if none was specified
if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release CACHE STRING 
      "Choose the type of build, options are: Release|Debug|RelWithDebInfo (for distros)." FORCE)
endif()
message(STATUS "Build type: ${CMAKE_BUILD_TYPE}")
message(STATUS "Install directory: ${CMAKE_INSTALL_PREFIX}")
################################################################################
# Version number
execute_process(
        COMMAND git rev-parse --abbrev-ref HEAD
        WORKING_DIRECTORY "${CMAKE_CURRENT_LIST_DIR}"
        OUTPUT_VARIABLE GIT_BRANCH
        OUTPUT_STRIP_TRAILING_WHITESPACE)

execute_process(
        COMMAND git log -1 --format=%h
        WORKING_DIRECTORY "${CMAKE_CURRENT_LIST_DIR}"
        OUTPUT_VARIABLE GIT_COMMIT_HASH
        OUTPUT_STRIP_TRAILING_WHITESPACE)

message(STATUS "Git current branch: ${GIT_BRANCH}")
message(STATUS "Git commit hash: ${GIT_COMMIT_HASH}")

message(STATUS "Generating version.hpp")
configure_file(
        ${CMAKE_CURRENT_LIST_DIR}/include/version.hpp.in
        ${CMAKE_CURRENT_BINARY_DIR}/generated/version.hpp)

include_directories(${CMAKE_CURRENT_BINARY_DIR}/generated)

################################################################################
# Compiler and linker flags

if (CMAKE_BUILD_TYPE STREQUAL "Release")
    set(BUILD_TYPE_FLAGS "-O3")
elseif(CMAKE_BUILD_TYPE STREQUAL "Debug")
    set(BUILD_TYPE_FLAGS "-g")
elseif(CMAKE_BUILD_TYPE STREQUAL "RelWithDebInfo")
    set(BUILD_TYPE_FLAGS "-O3 -g -DNDEBUG")
endif()

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${BUILD_TYPE_FLAGS}")

################################################################################
# Check dependencies

find_package(OpenMP)
if(OPENMP_FOUND)
    include_directories(${OpenMP_CXX_INCLUDE_DIRS})
    link_libraries(${OpenMP_CXX_LIBRARIES})
    if(${CMAKE_SYSTEM_NAME} MATCHES "Linux")
        set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS} -fopenmp")
        set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS} -fopenmp")
        set (CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
    elseif(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
        # assumes clang build
        # we can't reliably detect when we're using clang, so for the time being we assume
        set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS} -Xpreprocessor -fopenmp")
        set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS} -Xpreprocessor -fopenmp")
        set (CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
    endif()
endif()

# HTS
list(APPEND CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/cmake")
find_package(ZLIB REQUIRED)

if(NOT HTSlib_LIBRARIES)
    find_package(HTSlib REQUIRED)
else()
    message(STATUS "Bypassing FindHTSlib.cmake. Using injected HTSlib: ${HTSlib_LIBRARIES}")
endif()

if(HTSlib_FOUND OR HTSlib_LIBRARIES)
    # Added the /htslib subfolder to prevent the <hts.h> missing error
    include_directories(${HTSlib_INCLUDE_DIRS} ${HTSlib_INCLUDE_DIRS}/htslib)
    message(STATUS "HTSlib linked successfully.")
endif()

if(NOT TARGET htslib)
    add_library(htslib UNKNOWN IMPORTED)
    set_target_properties(htslib PROPERTIES
        IMPORTED_LOCATION "${HTSlib_LIBRARIES}"
        INTERFACE_INCLUDE_DIRECTORIES "${HTSlib_INCLUDE_DIRS}"
    )
endif()

option(PFP_USE_INTERNAL_LEVIOSAM "Build/fetch LevioSAM inside PFP" ON)

if(PFP_USE_INTERNAL_LEVIOSAM)
    if(TARGET lvsam)
        message(STATUS "Using caller-provided LevioSAM target: lvsam")
    else()
        find_library(LEVIOSAM_LIB lvsam)
        if(LEVIOSAM_LIB)
            message(STATUS "Found external leviosam library: ${LEVIOSAM_LIB}")
            add_library(lvsam INTERFACE IMPORTED)
            set_target_properties(lvsam PROPERTIES
                INTERFACE_LINK_LIBRARIES "${LEVIOSAM_LIB}"
            )
        else()
            message(STATUS "leviosam library not found. Building as an external content.")
            include(FetchContent)
            
            FetchContent_Declare(
                leviosam
                GIT_REPOSITORY https://github.com/maxrossi91/levioSAM.git
                GIT_TAG master
            )
            
            FetchContent_GetProperties(leviosam)
            if(NOT leviosam_POPULATED)
                FetchContent_Populate(leviosam)
                set(BUILD_TESTS OFF CACHE BOOL "Avoid building LevioSAM tests")

                add_subdirectory(${leviosam_SOURCE_DIR} ${leviosam_BINARY_DIR} EXCLUDE_FROM_ALL)
                cmake_policy(SET CMP0079 NEW)
                target_include_directories(lvsam PUBLIC ${leviosam_SOURCE_DIR}/src ${CMAKE_INCLUDE_PATH})
                if(SDSL_LIB)
                    target_link_libraries(lvsam INTERFACE htslib ${SDSL_LIB})
                else()
                    target_link_libraries(lvsam INTERFACE htslib)
                endif()
            endif()
        endif()
    endif()
else()
    message(STATUS "Skipping internal LevioSAM integration in PFP")
    if(NOT TARGET lvsam)
        add_library(lvsam INTERFACE)
    endif()
endif()
################################################################################
# Lib

include_directories(${CMAKE_CURRENT_LIST_DIR}/thirdparty)
include_directories(${CMAKE_CURRENT_LIST_DIR}/include)
include_directories(${CMAKE_CURRENT_BINARY_DIR}/generated)

file(GLOB SOURCES ${CMAKE_CURRENT_LIST_DIR}/src/*.cpp)
file(GLOB MURMUR_SOURCES ${CMAKE_CURRENT_LIST_DIR}/thirdparty/murmur/*.cpp)

# 32 and 64 bits parse
add_library(pfp ${SOURCES} ${MURMUR_SOURCES})
add_library(pfp64 ${SOURCES} ${MURMUR_SOURCES})
target_link_libraries(pfp lvsam)
target_link_libraries(pfp64 lvsam)
target_compile_options(pfp64 PRIVATE "-DPFP_LONG_TYPE=ON")

set(VCF_LIB_DEPS htslib ${ZLIB_LIBRARIES} dl pthread lvsam)
################################################################################
# Tests
# 
# add_executable(unit_tests tests/unit_tests.cpp)
# target_link_libraries(unit_tests pfp ${VCF_LIB_DEPS})
# 
# add_executable(unit_tests64 tests/unit_tests.cpp)
# target_compile_options(unit_tests64 PRIVATE "-DPFP_LONG_TYPE=ON")
# target_link_libraries(unit_tests64 pfp64 ${VCF_LIB_DEPS})
# 
################################################################################
# Tools

# PFP++

# 32 and 64 bits parse
add_executable(pfp++ pfp++.cpp)
target_link_libraries(pfp++ pfp ${VCF_LIB_DEPS})
add_executable(pfp++64 pfp++.cpp)
target_compile_options(pfp++64 PRIVATE "-DPFP_LONG_TYPE=ON")
target_link_libraries(pfp++64 pfp64 ${VCF_LIB_DEPS})

# AuPair
add_executable(AuPair au_pair.cpp)
target_link_libraries(AuPair pfp ${VCF_LIB_DEPS})
add_executable(AuPair64 au_pair.cpp)
target_link_libraries(AuPair64 pfp64 ${VCF_LIB_DEPS})

# PFP++ merge
add_executable(mpfp++ merge.cpp)
target_link_libraries(mpfp++ pfp ${VCF_LIB_DEPS})
add_executable(mpfp++64 merge.cpp)
target_compile_options(mpfp++64 PRIVATE "-DPFP_LONG_TYPE=ON")
target_link_libraries(mpfp++64 pfp64 ${VCF_LIB_DEPS})

# VCT to Fasta
add_executable(vcf_to_fa vcf_to_fa.cpp)
target_link_libraries(vcf_to_fa pfp ${VCF_LIB_DEPS})

# Unparse
add_executable(pfp++_unparse unparse.cpp)
target_link_libraries(pfp++_unparse pfp ${VCF_LIB_DEPS})
add_executable(pfp++64_unparse unparse.cpp)
target_compile_options(pfp++64_unparse PRIVATE "-DPFP_LONG_TYPE=ON")
target_link_libraries(pfp++64_unparse pfp64 ${VCF_LIB_DEPS})

# Integrity Check
add_executable(check check_integrity.cpp)
target_link_libraries(check pfp ${VCF_LIB_DEPS})
add_executable(check64 check_integrity.cpp)
target_compile_options(check64 PRIVATE "-DPFP_LONG_TYPE=ON")
target_link_libraries(check64 pfp64 ${VCF_LIB_DEPS})

install(TARGETS pfp++ pfp++64 AuPair AuPair64 mpfp++ mpfp++64 pfp++_unparse pfp++64_unparse check check64 TYPE RUNTIME)