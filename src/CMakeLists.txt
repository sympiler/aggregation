cmake_minimum_required( VERSION 3.0)
project(lbc)

include_directories(includes
        ${SPARSE_UTIL_INC}
        ${OpenMP_CXX_INCLUDE_DIRS}
        ${METIS_INCLUDE_DIRS}
        )

add_library (lbc STATIC
        lbc_tree.cpp
        lbc_csc.cpp
        lbc_utils.cpp
        ${SPARSE_UTIL_SRC_FILES}
        )
target_link_libraries(lbc PRIVATE
        sym_amd)