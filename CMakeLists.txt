cmake_minimum_required(VERSION 2.8)
project(LBC)
set(CMAKE_CXX_STANDARD 11)

set(METIS_INC "$ENV{METISROOT}/include")
set(METIS_LIB "$ENV{METISROOT}/libmetis")

if (APPLE)
    set(MKL_INC "$ENV{MKLROOT}/include")
    set(MKL_LIB "$ENV{MKLROOT}/lib")
    set(MKL_LIB2 "$ENV{MKLROOT}/../compiler/lib")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fopenmp -O3 -g")
endif ()
if (UNIX AND NOT APPLE)
    set(MKL_INC "$ENV{MKLROOT}/include")
    set(MKL_LIB "$ENV{MKLROOT}/lib/intel64")
    # for iomp lib
    set(MKL_LIB2 "$ENV{MKLROOT}/../compiler/lib/intel64")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fopenmp -std=c++11 -O3 -g")
endif ()

#add_subdirectory(src)
add_subdirectory(triangular_example)