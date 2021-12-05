cmake_minimum_required(VERSION 3.12)

list(APPEND CMAKE_MODULE_PATH
        ${CMAKE_CURRENT_LIST_DIR}/cmake
        "${CMAKE_CURRENT_LIST_DIR}/cmake/third_party"
        )

include(DownloadProject)


set(LBC_ROOT_DIR ${CMAKE_CURRENT_LIST_DIR})
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_CURRENT_LIST_DIR}/lib/)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_CURRENT_LIST_DIR}/lib/)
set(LBC_EXTERNAL ${CMAKE_CURRENT_LIST_DIR}/external/)


project(LBC_lib CXX C)
set(CMAKE_CXX_STANDARD 11)

get_directory_property(HAS_PARENT PARENT_DIRECTORY)
if(HAS_PARENT)
    set(LBC_IS_TOPLEVEL OFF)
else()
    set(LBC_IS_TOPLEVEL ON)
endif()

option(LBC_BUILD_EXAMPLES "Build LBC Examples" ${LBC_IS_TOPLEVEL})
option(LBC_BUILD_TESTS "Build LBC TESTs" ${LBC_IS_TOPLEVEL})
option(LBC_PROFILING "Build LBC Profiling Routines" ${LBC_IS_TOPLEVEL})

#TODO: remove mmio and use the sparse_io in smp repo
add_compile_options(-Wno-write-strings)
add_compile_options(-Wunused-result)

set(SPARSE_UTIL_SRC "${CMAKE_CURRENT_SOURCE_DIR}/utils")

set(Sparse_BLAS_INC "${CMAKE_CURRENT_SOURCE_DIR}/sparse_blas/includes/")
set(Sparse_BLAS_SRC "${CMAKE_CURRENT_SOURCE_DIR}/sparse_blas/")
AUX_SOURCE_DIRECTORY(${Sparse_BLAS_SRC} Sparse_BLAS_SRC_FILES)

set(LBC_INC "${CMAKE_CURRENT_SOURCE_DIR}/include")
set(LBC_SRC "${CMAKE_CURRENT_SOURCE_DIR}/src")
AUX_SOURCE_DIRECTORY(${LBC_SRC} LBC_SRC_FILES)

find_package(OpenMP)
if(OpenMP_FOUND)
    if(APPLE) #TODO: there might be a better support
        set(OpenMP_CXX_INCLUDE_DIRS  "/usr/local/include/")
    endif()
    include_directories(
            ${OpenMP_CXX_INCLUDE_DIRS}
    )
    link_libraries(
            OpenMP::OpenMP_CXX
    )
    add_definitions(-DENABLE_OPENMP)
    add_definitions(-D__OPENMP__)
endif()


include(metis)
if(METIS_FOUND)
    add_definitions(-DMETIS)
    set(METIS_SOURCE_FILES "${SPARSE_UTIL_SRC}/metis_interface.cpp")
endif()

### Find PAPI for Profiling
### Can be cloned directly from https://bitbucket.org/icl/papi.git
##TODO: create a cmake header
#find_package(PAPI)
if (PAPI_FOUND)
    message(STATUS "PAPI found")
    include_directories(
            ${PAPI_INCLUDE_DIRS}
    )
    add_definitions(-DPAPI)
    if(LBC_PROFILING)
        add_definitions(-DPROFILE)
        message(STATUS "Profiling is enabled.")
    endif()
endif()

file(GLOB_RECURSE EXT_SRC ${CMAKE_CURRENT_SOURCE_DIR}/utils/external/*.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/utils/external/*.c)
add_library(sym_amd STATIC
        ${EXT_SRC})
target_include_directories(sym_amd PRIVATE
        ${LBC_INC}
        )
set(SPARSE_UTIL_SRC_FILES
        "${CMAKE_CURRENT_SOURCE_DIR}/utils/sparse_inspector.cpp"
        "${CMAKE_CURRENT_SOURCE_DIR}/utils/sparse_io.cpp"
        "${CMAKE_CURRENT_SOURCE_DIR}/utils/sparse_utilities.cpp"
        "${CMAKE_CURRENT_SOURCE_DIR}/utils/test_utils.cpp"
        "${CMAKE_CURRENT_SOURCE_DIR}/utils/utils.cpp"
        "${CMAKE_CURRENT_SOURCE_DIR}/utils/lfactor_creation.cpp"
        "${CMAKE_CURRENT_SOURCE_DIR}/utils/etree.cpp"
        "${CMAKE_CURRENT_SOURCE_DIR}/utils/FusionDemo.cpp"
        "${CMAKE_CURRENT_SOURCE_DIR}/utils/BCSCMatrix.cpp"
        "${CMAKE_CURRENT_SOURCE_DIR}/utils/external/mmio.cpp"
        "${METIS_SOURCE_FILES}"
        )

add_subdirectory(src)


#if(EXISTS "$ENV{PAPIROOT}")
#    if(EXISTS "$ENV{PAPIROOT}/src")
#        set(PAPI_LIB "$ENV{PAPIROOT}/src")
#        set(PAPI_LIB2 "${PAPI_LIB}/libpfm4/lib/")
#    elseif(EXISTS "$ENV{PAPIROOT}/lib")
#        set(PAPI_LIB "$ENV{PAPIROOT}/lib")
#    endif()
#    if(EXISTS "$ENV{PAPIROOT}/include")
#        set(PAPI_INC "$ENV{PAPIROOT}/include")
#    else()
#        set(PAPI_INC "$ENV{PAPIROOT}/src")
#    endif()
#
#    add_definitions(-DPAPI)
#    set(PAPI_LIB_NAME "papi")
#    set(PAPI_LIB_NAME2 "pfm")
#    message(STATUS "PAPI found")
#else()
#    message(STATUS "PAPI not found, building without PAPI.")
#    set(PAPI_LIB_NAME "")
#endif()

if(LBC_PROFILING AND PAPI_FOUND)
    add_subdirectory(profiler)
endif()

if(LBC_BUILD_EXAMPLES)
    add_subdirectory(example)
endif()

if(LBC_BUILD_TESTS)
    add_subdirectory(Catch_tests)
endif()
