cmake_minimum_required(VERSION 2.8)

list(APPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_LIST_DIR}/cmake)

project(LBC_lib)
set(CMAKE_CXX_STANDARD 11)

#TODO: remove mmio and use the sparse_io in smp repo
add_compile_options(-Wno-write-strings)
add_compile_options(-Wunused-result)

set(SPARSE_UTIL_INC "${CMAKE_CURRENT_SOURCE_DIR}/utils/includes")
set(SPARSE_UTIL_SRC "${CMAKE_CURRENT_SOURCE_DIR}/utils")

#AUX_SOURCE_DIRECTORY(${SPARSE_UTIL_SRC} SPARSE_UTIL_SRC_FILES)
set(SPARSE_UTIL_SRC_FILES
        "${CMAKE_CURRENT_SOURCE_DIR}/utils/sparse_inspector.cpp"
        "${CMAKE_CURRENT_SOURCE_DIR}/utils/sparse_io.cpp"
        "${CMAKE_CURRENT_SOURCE_DIR}/utils/sparse_utilities.cpp"
        "${CMAKE_CURRENT_SOURCE_DIR}/utils/test_utils.cpp"
        "${CMAKE_CURRENT_SOURCE_DIR}/utils/utils.cpp")


set(SPARSE_EXTN_INC "${SPARSE_UTIL_SRC}/external/includes")
set(SPARSE_EXTN_SRC "${SPARSE_UTIL_SRC}/external")
AUX_SOURCE_DIRECTORY(${SPARSE_EXTN_SRC} SPARSE_EXTN_SRC_FILES)

set(LBC_INC "${CMAKE_CURRENT_SOURCE_DIR}/lbc/includes")
set(LBC_SRC "${CMAKE_CURRENT_SOURCE_DIR}/lbc")

set(Group_INC "${CMAKE_CURRENT_SOURCE_DIR}/Group/includes/")

AUX_SOURCE_DIRECTORY(${LBC_SRC} LBC_SRC_FILES)


find_package(METIS OPTIONAL_COMPONENTS)

if (APPLE)
    SET(CMAKE_C_COMPILER /usr/local/Cellar/gcc/9.1.0/bin/gcc-9)
    SET(CMAKE_CXX_COMPILER /usr/local/Cellar/gcc/9.1.0/bin/g++-9)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fopenmp ")
endif ()
if (UNIX AND NOT APPLE)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fopenmp -std=c++11 ")
endif ()

add_subdirectory(lbc)
if(EXISTS "${METIS_INCLUDES}")
    add_definitions(-DMETIS)
    set(METIS_SOURCE_FILES "${SPARSE_UTIL_SRC}/metis_interface.cpp")
    add_subdirectory(example)
endif()
