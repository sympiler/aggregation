cmake_minimum_required(VERSION 3.0)




if(METIS_FOUND)

include_directories(includes
        ${SPARSE_UTIL_INC}
        ${SPARSE_EXTN_INC}
        ${LBC_INC}
        ${METIS_INCLUDE_DIRS} ${METIS_INCLUDES}
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

link_directories(${METIS_LIBRARY})
if(METIS_FOUND)
target_link_libraries(lbc_demo
        lbc
        ${METIS_LIBRARIES}
        )
else()
target_link_libraries(lbc_demo
        lbc)
endif()
