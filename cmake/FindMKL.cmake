# Try to find MKL headers and libraries.
#
# Usage of this module as follows:
#
#     find_package(MKL)
#
# Variables used by this module, they can change the default behaviour and need
# to be set before calling find_package:
#
#  MKL_PREFIX         Set this variable to the root installation of
#                      libpapi if the module has problems finding the
#                      proper installation path.
#
# Variables defined by this module:
#
#  MKL_FOUND              System has MKL libraries and headers
#  MKL_LIBRARIES          The MKL library
#  MKL_INCLUDE_DIRS       The location of PAPI headers

find_path(MKL_PERFIX
    NAMES include/mkl.h
    PATHS $ENV{MKLROOT}
)

#compiler/2021.1.1/linux/compiler/lib/intel64_lin
find_library(MKL_GNU_THREAD NAMES mkl_intel_lp64 HINTS ${MKL_PERFIX}/lib/intel64  ${HILTIDEPS}/lib)
find_library(PTHREAD NAMES pthread HINTS ${MKL_PERFIX}/lib/intel64  ${HILTIDEPS}/lib)
find_library(MKL_INTEL_THREAD NAMES mkl_intel_thread HINTS ${MKL_PERFIX}/lib/intel64  ${HILTIDEPS}/lib)
find_library(MKL_SEQUENTIAL NAMES mkl_sequential HINTS ${MKL_PERFIX}/lib/intel64  ${HILTIDEPS}/lib)
find_library(MKL_CORE NAMES mkl_core HINTS ${MKL_PERFIX}/lib/intel64  ${HILTIDEPS}/lib)
find_library(M NAMES m HINTS ${MKL_PERFIX}/lib/intel64  ${HILTIDEPS}/lib)
find_library(DL NAMES dl HINTS ${MKL_PERFIX}/lib/intel64  ${HILTIDEPS}/lib)

find_library(IOMP5 NAMES iomp5 HINTS ${MKL_PERFIX}/../../compiler/latest/linux/compiler/lib/intel64 ${MKL_PERFIX}/../compilers_and_libraries/linux/lib/intel64 ${HILTIDEPS}/lib)

find_path(MKL_INCLUDE_DIRS
        NAMES mkl.h
        HINTS ${MKL_PERFIX}/include ${HILTIDEPS}/include
        )

add_definitions(-DMKL)

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(MKL DEFAULT_MSG
        MKL_GNU_THREAD
        PTHREAD
        MKL_INTEL_THREAD
        IOMP5
        MKL_SEQUENTIAL
        MKL_CORE
        M
        DL
        MKL_INCLUDE_DIRS
        )

mark_as_advanced(
        MKL_PREFIX_DIRS
        MKL_GNU_THREAD
        PTHREAD
        MKL_INTEL_THREAD
        IOMP5
        MKL_SEQUENTIAL
        MKL_CORE
        M
        DL
        MKL_INCLUDE_DIRS
)
