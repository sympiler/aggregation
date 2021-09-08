cmake_minimum_required(VERSION 3.0)

list(APPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_LIST_DIR}/cmake)

include(DownloadProject)


set(LBC_ROOT_DIR ${CMAKE_CURRENT_LIST_DIR})
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_CURRENT_LIST_DIR}/lib/)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_CURRENT_LIST_DIR}/lib/)
set(LBC_EXTERNAL ${CMAKE_CURRENT_LIST_DIR}/external)


project(LBC_lib)
set(CMAKE_CXX_STANDARD 11)

#TODO: remove mmio and use the sparse_io in smp repo
add_compile_options(-Wno-write-strings)
add_compile_options(-Wunused-result)

set(SPARSE_UTIL_INC "${CMAKE_CURRENT_SOURCE_DIR}/includes")
set(SPARSE_UTIL_SRC "${CMAKE_CURRENT_SOURCE_DIR}/utils")

#AUX_SOURCE_DIRECTORY(${SPARSE_UTIL_SRC} SPARSE_UTIL_SRC_FILES)


set(SPARSE_EXTN_INC "${SPARSE_UTIL_SRC}/external/includes")
set(SPARSE_EXTN_SRC "${SPARSE_UTIL_SRC}/external")
AUX_SOURCE_DIRECTORY(${SPARSE_EXTN_SRC} SPARSE_EXTN_SRC_FILES)

set(LBC_INC "${CMAKE_CURRENT_SOURCE_DIR}/includes")
set(LBC_SRC "${CMAKE_CURRENT_SOURCE_DIR}/src")
AUX_SOURCE_DIRECTORY(${LBC_SRC} LBC_SRC_FILES)


find_package(METIS OPTIONAL_COMPONENTS)
if(NOT METIS_FOUND)
    message(STATUS  "METIS not detected, downloading it..." )
#    find_program(MAKE_EXE NAMES gmake nmake make)
#    set(metis_INSTALL_DIR "${CMAKE_CURRENT_BINARY_DIR}/metis")
#    set(metis_CMAKE_ARGS -DCMAKE_INSTALL_PREFIX=${metis_INSTALL_DIR})
#set(EXTERNAL_INSTALL_LOCATION ${CMAKE_CURRENT_LIST_DIR}/external)
# set(GKLIB_PATH ${SYMPILER_EXTERNAL}/metis/metis-src/GKLIB/ CACHE PATH "d")
download_project(PROJ metis_shared
        GIT_REPOSITORY https://github.com/cheshmi/METIS.git  # https://github.com/KarypisLab/METIS.git
        # BUILD_IN_SOURCE 1
        #PREFIX  ${metis_BINARY_DIR}
        #     URL http://glaros.dtc.umn.edu/gkhome/fetch/sw/metis/metis-5.1.0.tar.gz
        PREFIX  ${LBC_EXTERNAL}/metis/
        #SOURCE_DIR ${SYMPILER_EXTERNAL}/metis/metis-src/
        BINARY_DIR ${CMAKE_CURRENT_BINARY_DIR}/metis/
        UPDATE_DISCONNECTED 1
        #CMAKE_ARGS  -DGKLIB_PATH=${metis_SOURCE_DIR}/GKlib/
        # BUILD_COMMAND ${MAKE_EXE}
        )

    set(SHARED  TRUE CACHE BOOL "shared library" FORCE)
    add_subdirectory(${metis_shared_SOURCE_DIR})
    set(METIS_INCLUDE_DIRS ${metis_shared_SOURCE_DIR}/include)
    set(METIS_INCLUDES ${metis_shared_SOURCE_DIR}/include/)
    set(METIS_LIBRARY_DIRS ${LBC_ROOT_DIR}/lib/)
    set(METIS_LIBRARIES ${LBC_ROOT_DIR}/lib/libmetis.so)
    set(METIS_FOUND TRUE)
    message(STATUS "METIS library is found at ${METIS_LIBRARY}")
endif()

if(METIS_FOUND)
    add_definitions(-DMETIS)
    set(METIS_SOURCE_FILES "${SPARSE_UTIL_SRC}/metis_interface.cpp")
endif()

file(GLOB_RECURSE EXT_SRC ${CMAKE_CURRENT_SOURCE_DIR}/utils/external/*.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/utils/external/*.c)
add_library(sym_amd STATIC
        ${EXT_SRC})

set(SPARSE_UTIL_SRC_FILES
        "${CMAKE_CURRENT_SOURCE_DIR}/utils/sparse_inspector.cpp"
        "${CMAKE_CURRENT_SOURCE_DIR}/utils/sparse_io.cpp"
        "${CMAKE_CURRENT_SOURCE_DIR}/utils/sparse_utilities.cpp"
        "${CMAKE_CURRENT_SOURCE_DIR}/utils/test_utils.cpp"
        "${CMAKE_CURRENT_SOURCE_DIR}/utils/utils.cpp"
        "${CMAKE_CURRENT_SOURCE_DIR}/utils/lfactor_creation.cpp"
        "${CMAKE_CURRENT_SOURCE_DIR}/utils/etree.cpp"
        "${CMAKE_CURRENT_SOURCE_DIR}/utils/BCSCMatrix.cpp"
        "${CMAKE_CURRENT_SOURCE_DIR}/utils/external/mmio.cpp"
        "${METIS_SOURCE_FILES}"
        )

add_subdirectory(src)
add_subdirectory(example)
