cmake_minimum_required(VERSION 2.8)



if(EXISTS "${METIS_INCLUDES}")
include_directories(includes
        ${SPARSE_UTIL_INC}
        ${SPARSE_EXTN_INC}
        ${LBC_INC}
        ${METIS_INCLUDE_DIR} ${METIS_INCLUDES}
        )
else()
include_directories(includes
        ${SPARSE_UTIL_INC}
        ${SPARSE_EXTN_INC}
        ${LBC_INC}
        )
endif()


project(lbc_demo)
set(LBCD_SOURCE_FILES
        ${SPARSE_UTIL_SRC}/sparse_utilities.cpp
        ${SPARSE_UTIL_SRC}/sparse_inspector.cpp
        ${SPARSE_UTIL_SRC}/sparse_io.cpp
        ${SPARSE_UTIL_SRC}/test_utils.cpp
        ${SPARSE_UTIL_SRC}/utils.cpp
        ${SPARSE_UTIL_SRC}/BCSCMatrix.cpp
        ${SPARSE_EXTN_SRC}/mmio.cpp

        ${METIS_SOURCE_FILES}
        )

add_executable(lbc_demo ${LBCD_SOURCE_FILES}
        lbc_demo.cpp )
if (NOT METIS_FOUND)
add_dependencies(lbc_demo metis)
endif ()
if(METIS_FOUND)
target_link_libraries(lbc_demo
        lbc
        ${METIS_LIBRARY}
        )
else()
target_link_libraries(lbc_demo
        lbc)
endif()
