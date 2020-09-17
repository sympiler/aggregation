cmake_minimum_required(VERSION 2.8)

include_directories(includes
        ${SPARSE_UTIL_INC}
        ${SPARSE_EXTN_INC}
        ${LBC_INC}
        ${METIS_INCLUDE_DIR} ${METIS_INCLUDES}
        ${Group_INC}
        )

project(sptrsv_demo)
set(SPTRSVD_SOURCE_FILES
        sptrsv.cpp
        ${SPARSE_UTIL_SRC}/sparse_utilities.cpp
        ${SPARSE_UTIL_SRC}/sparse_inspector.cpp
        ${SPARSE_UTIL_SRC}/sparse_io.cpp
        ${SPARSE_UTIL_SRC}/test_utils.cpp
        ${SPARSE_UTIL_SRC}/utils.cpp
        ${SPARSE_EXTN_SRC}/mmio.cpp
        ${Group_SRC}/StatSpMat.cpp
        ${Group_SRC}/executor.cpp
        ${Group_SRC}/Group.cpp
        ${Group_SRC}/Utils.cpp
        ${METIS_SOURCE_FILES}
        FusionDemo.cpp
        )

add_executable(sptrsv_demo ${SPTRSVD_SOURCE_FILES}
        sptrsv_demo.cpp )

target_link_libraries(sptrsv_demo
        lbc
        ${METIS_LIBRARY} ${METIS_LIBRARIES})

add_executable(sptrsv_csr_profile ${SPTRSVD_SOURCE_FILES}
        sptrsv_csr_profile.cpp)
target_link_libraries( sptrsv_csr_profile lbc
        ${METIS_LIBRARY} ${METIS_LIBRARIES})