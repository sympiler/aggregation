
include_directories(
        ${LBC_INC}
)

project(sptrsv_profiler_demo)

add_executable(sptrsv_profiler_demo
        sptrsv_profiler.cpp
        ${Sparse_BLAS_SRC_FILES}
        )

target_include_directories(sptrsv_profiler_demo PRIVATE
        ${Sparse_BLAS_INC}
        )
target_link_libraries(sptrsv_profiler_demo PRIVATE
        metis::metis
        ${PAPI_LIBRARIES}
        lbc
        )
if(OpenMP_FOUND)
    target_link_libraries(sptrsv_profiler_demo PRIVATE OpenMP::OpenMP_CXX)
endif()
