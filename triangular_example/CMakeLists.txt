cmake_minimum_required(VERSION 2.8)
project(LBCTest)


include_directories(${MKL_INC}  ${METIS_INC}
        "../src/")
link_directories(${MKL_LIB} ${METIS_LIB})
add_executable(LBCTest lbc_test.cpp )

target_link_libraries(LBCTest
        #MKL BLAS
        mkl_intel_lp64
        mkl_sequential
        mkl_core
        mkl_intel_ilp64
        #mkl_gnu_thread
        pthread
        m
        dl
        metis
        )