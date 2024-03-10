cmake_minimum_required( VERSION 3.0)
project(lbc)

add_library (lbc STATIC
        ${LBC_SRC_FILES}
        ${SPARSE_UTIL_SRC_FILES}
        )



target_include_directories(lbc PRIVATE
        ${LBC_INC}
)

target_link_libraries(lbc PRIVATE
        sym_amd
        metis::metis
        )

add_library (lbcShared SHARED
        ${LBC_SRC_FILES}
        ${SPARSE_UTIL_SRC_FILES}
        )
target_include_directories(lbcShared PRIVATE
        ${LBC_INC}
        )

target_link_libraries(lbcShared PRIVATE
        sym_amd
        metis::metis
        )
########## HDAGG
add_library(aggregation
        ${HDAGG_SRCFILES}
        )
target_include_directories(aggregation PRIVATE
        ${LBC_INC}
        )
target_link_libraries(aggregation PRIVATE
        lbc
        sym_amd
        metis::metis
        )


add_library(aggregationShared
        ${HDAGG_SRCFILES}
        )
target_include_directories(aggregationShared PRIVATE
        ${LBC_INC}
        )
target_link_libraries(aggregationShared PRIVATE
        lbcShared
        sym_amd
        metis::metis
        )
### Find Intel MKL if available - it works with oneAPI
if(HDAGG_WITH_MKL)
    find_package(MKL OPTIONAL_COMPONENTS)
    if(EXISTS "${MKL_INCLUDE_DIRS}")
        target_include_directories(aggregation PRIVATE
                ${MKL_INCLUDE_DIRS}
                )
        target_link_libraries(aggregation PRIVATE
                ${MKL_GNU_THREAD}
                ${PTHREAD}
                ${IOMP5}
                ${MKL_SEQUENTIAL}
                ${MKL_INTEL_THREAD}
                ${MKL_CORE}
                ${M}
                ${DL}
                )
        add_definitions(-DMKL)

    else()
        message(STATUS "Intel MKL not found, building without MKL.")
        set(MKL 0)
    endif()
endif()

### Check whether DAGP root exists - Relative works
if(HDAGG_WITH_DAGP)
    if(EXISTS "$ENV{DAGPROOT}")
        message(STATUS "The DAGP is found - $ENV{DAGPROOT}.")
        target_include_directories(aggregation PRIVATE
                "$ENV{DAGPROOT}/src/recBisection/"
                "$ENV{DAGPROOT}/src/common/"
                "${CMAKE_CURRENT_SOURCE_DIR}/utils/external/DAGPInterface/includes"
                )
        target_link_libraries(aggregation PRIVATE
                "$ENV{DAGPROOT}/lib/libdagp.a"
                "${CMAKE_CURRENT_SOURCE_DIR}/utils/external/DAGPInterface/*.cpp"
                )
        add_definitions(-DDAGP)
    else()
        message(STATUS "Please defined the DAGP directory using DAGP_ROOT variable")
    endif()
endif()

### Check whether SPMP root exists - Relative works
if(HDAGG_WITH_SPMP)
    if(EXISTS "$ENV{SPMPROOT}")
        message(STATUS "Found SpMP - $ENV{SPMPROOT}.")
        target_include_directories(aggregation PRIVATE
                $ENV{SPMPROOT}
                "${CMAKE_CURRENT_SOURCE_DIR}/utils/external/SpMPInterface/includes"
                )
        target_link_libraries(aggregation PRIVATE
                "$ENV{SPMPROOT}/libspmp.a"
                )
        add_definitions(-DSPMP)
    else()
        message(STATUS "Please defined the SpMP directory using SPMP_ROOT variable - not found")
    endif()
endif()
if(OpenMP_FOUND)
target_link_libraries(lbc PRIVATE OpenMP::OpenMP_CXX)
target_link_libraries(aggregation PRIVATE OpenMP::OpenMP_CXX)
endif()