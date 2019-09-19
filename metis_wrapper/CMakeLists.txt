cmake_minimum_required(VERSION 2.8)
project(metisTest)

set(SOURCE_FILES metis_testing.cpp)

include_directories(${METIS_INC} )
link_directories(${METIS_LIB})

add_executable(metisTest ${SOURCE_FILES} metis_testing.cpp)

target_link_libraries(metisTest
        #gomp
        m
        dl
        metis
        )