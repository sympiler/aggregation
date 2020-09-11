cmake_minimum_required( VERSION 2.5)
project(lbc)

include_directories(includes
        ${SPARSE_UTIL_INC}
        )

add_library (lbc STATIC
        lbc_csc.cpp
        lbc_utils.cpp
        lbc_tree.cpp
        ${SPARSE_UTIL_SRC}/sparse_utilities.cpp
        ${SPARSE_UTIL_SRC}/etree.cpp
        ${SPARSE_UTIL_SRC}/utils.cpp)