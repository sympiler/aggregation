cmake_minimum_required(VERSION 2.8)

include_directories(includes
        ${SPARSE_UTIL_INC}
        ${SPARSE_EXTN_INC}
        ${LBC_INC}
        ${METIS_INCLUDE_DIR} ${METIS_INCLUDES}
        ${Group_INC}
        ${Sparse_BLAS_INC}
        ${Profiler_INC}
        )

project(sptrsv_demo)
set(SPTRSVD_SOURCE_FILES
        ${SPARSE_UTIL_SRC}/sparse_utilities.cpp
        ${SPARSE_UTIL_SRC}/sparse_inspector.cpp
        ${SPARSE_UTIL_SRC}/sparse_io.cpp
        ${SPARSE_UTIL_SRC}/test_utils.cpp
        ${SPARSE_UTIL_SRC}/utils.cpp
        ${SPARSE_EXTN_SRC}/mmio.cpp
        ${Profiler_SRC}/StatSpMat.cpp
        ${Sparse_BLAS_SRC}/executor.cpp
        ${Sparse_BLAS_SRC}/sptrsv.cpp
        ${Group_SRC}/Group.cpp
        ${Group_SRC}/Utils.cpp
        ${METIS_SOURCE_FILES}
        FusionDemo.cpp
        )

add_executable(sptrsv_demo ${SPTRSVD_SOURCE_FILES}
        sptrsv_demo.cpp)
target_link_libraries(sptrsv_demo
        lbc
        ${METIS_LIBRARY} ${METIS_LIBRARIES}         ${PAPI_LIB}
        ${PAPI_LIB2} )


if(NOT ${PAPI})
    message(STATUS "PAPI not available, skipping ")
else()
add_executable(sptrsvprofiler_demo ${SPTRSVD_SOURCE_FILES}
        sptrsvprofiler_demo.cpp)
target_link_libraries(sptrsvprofiler_demo
        lbc
        ${METIS_LIBRARY} ${METIS_LIBRARIES}         ${PAPI_LIB}
        ${PAPI_LIB2})
endif()



if(NOT ${PAPI})
    message(STATUS "PAPI not available, skipping ")
else()
add_executable(sptrsv_csr_profile ${SPTRSVD_SOURCE_FILES}
        sptrsv_csr_profile.cpp)
target_link_libraries( sptrsv_csr_profile lbc
        ${METIS_LIBRARY} ${METIS_LIBRARIES} ${PAPI_LIB}
        ${PAPI_LIB2})
endif()

add_executable(sptrsv_dag_test ${SPTRSVD_SOURCE_FILES}
        sptrsv_dag_test.cpp)
target_link_libraries(sptrsv_dag_test
        lbc
        ${METIS_LIBRARY} ${METIS_LIBRARIES})


project(spic0_demo)
set(SPIC0_SOURCE_FILES
        ${SPARSE_UTIL_SRC}/sparse_utilities.cpp
        ${SPARSE_UTIL_SRC}/sparse_inspector.cpp
        ${SPARSE_UTIL_SRC}/sparse_io.cpp
        ${SPARSE_UTIL_SRC}/test_utils.cpp
        ${SPARSE_UTIL_SRC}/utils.cpp
        ${SPARSE_EXTN_SRC}/mmio.cpp
        ${Profiler_SRC}/StatSpMat_v1.cpp
        ${Sparse_BLAS_SRC}/executor.cpp
        ${Sparse_BLAS_SRC}/spic0.cpp
        ${Group_SRC}/Group.cpp
        ${Group_SRC}/Utils.cpp
        ${METIS_SOURCE_FILES}
        FusionDemo.cpp
        )

add_executable(spic0_demo ${SPIC0_SOURCE_FILES}
        spic0_demo.cpp)
target_link_libraries(spic0_demo
        lbc
        ${METIS_LIBRARY} ${METIS_LIBRARIES})

add_executable(spic0_csc_demo ${SPIC0_SOURCE_FILES}
        spico_csc_demo.cpp)

target_link_libraries(spic0_csc_demo
        lbc
        ${METIS_LIBRARY} ${METIS_LIBRARIES})

if(NOT ${PAPI})
    message(STATUS "PAPI not available, skipping ")
else()
add_executable(spic0_csc_profile  ${SPIC0_SOURCE_FILES}
                        spic0_csc_profile.cpp)
target_link_libraries(spic0_csc_profile
                lbc
        ${METIS_LIBRARY} ${METIS_LIBRARIES})
endif()