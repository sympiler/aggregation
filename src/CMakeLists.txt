cmake_minimum_required(VERSION 2.8)
project(LBC_LIB)

add_executable(lbc_test main.cpp csc.cpp csr.cpp matrix_format.cpp)

#TODO: build this as a static library
