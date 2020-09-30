include_directories(includes
        ${SPARSE_UTIL_INC}
        ${SPARSE_EXTN_INC}
        ${LBC_INC}
        ${Sparse_BLAS_INC}
        ${METIS_INCLUDE_DIR} ${METIS_INCLUDES}
        )

link_directories(
        ${SPARSE_LIB_LIB}
        ${SPARSE_UTIL_LIB}
        ${METIS_LIB}
        ${MKL_LIB}
        ${MKL_LIB2}
        ${PAPI_LIB}
        ${PAPI_LIB2})

set(SOURCE_FILES
        ${Sparse_BLAS_SRC}/sptrsv.cpp
        ${SPARSE_UTIL_SRC}/sparse_utilities.cpp
        ${SPARSE_UTIL_SRC}/sparse_inspector.cpp
        ${SPARSE_UTIL_SRC}/sparse_io.cpp
        ${SPARSE_UTIL_SRC}/test_utils.cpp
        ${SPARSE_UTIL_SRC}/utils.cpp
        ${SPARSE_EXTN_SRC}/mmio.cpp
        ${METIS_SOURCE_FILES}
        )

add_executable(Catch_tests ${SOURCE_FILES}
        sptrsv_tests.cpp
        )

target_link_libraries(Catch_tests
        lbc
        ${METIS_LIBRARY} ${METIS_LIBRARIES})






