cmake_minimum_required( VERSION 3.0)
project(lbc)

add_library (lbc STATIC
        ${LBC_SRC_FILES}
        ${SPARSE_UTIL_SRC_FILES}
        )

target_include_directories( lbc PRIVATE
        ${LBC_INC}
)

target_link_libraries(lbc PRIVATE
        sym_amd
        metis::metis
        )

if(OpenMP_FOUND)
target_link_libraries(lbc PRIVATE OpenMP::OpenMP_CXX)
endif()