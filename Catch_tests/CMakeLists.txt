
include_directories(
        ${LBC_INC}
        )

add_executable(Catch_tests
        ${Sparse_TRSV_SRC_FILES}
        sptrsv_tests.cpp
        )

target_include_directories(Catch_tests PRIVATE
        ${Sparse_TRSV_INC}
        )

target_link_libraries(Catch_tests
        lbc
        metis::metis
        )

add_executable(hdagg_tests
        ${Sparse_TRSV_SRC_FILES}
        hdagg_tests.cpp
        )

target_include_directories(hdagg_tests PRIVATE
        ${Sparse_TRSV_INC}
        )

target_link_libraries(hdagg_tests
        aggregation
        metis::metis
        )
