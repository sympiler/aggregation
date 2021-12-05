
include_directories(
        ${LBC_INC}
        )

add_executable(Catch_tests
        ${Sparse_BLAS_SRC_FILES}
        sptrsv_tests.cpp
        )

target_include_directories(Catch_tests PRIVATE
        ${Sparse_BLAS_INC}
        )

target_link_libraries(Catch_tests
        lbc
        metis::metis
        )

