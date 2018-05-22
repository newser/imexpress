# Copyright (C) 2017 haniu (niuhao.cn@gmail.com)
#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 3 of the License, or (at
# your option) any later version.
#
# This program is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
# General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software
# Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301,
# USA.

project(imexpress)

cmake_minimum_required(VERSION 3.1.3)

#
# definition
#

set(ROOT_PATH ${CMAKE_CURRENT_LIST_DIR})
set(OUT_PATH ${CMAKE_CURRENT_BINARY_DIR})

set(BUILD_PATH ${ROOT_PATH}/build)
set(LIBRARY_PATH ${ROOT_PATH}/library)

include(${BUILD_PATH}/environment.cmake)
include(${BUILD_PATH}/${ENV}.cmake)
include(${BUILD_PATH}/util.cmake)

option(DEVELOP_BUILD "build mode during developing" TRUE)
mark_as_advanced(DEVELOP_BUILD)

#
# source file
#

function(add_group group inc src f_list)
    # make group
    file(GLOB f_inc "${inc}/*.h")
    #message(STATUS "${group} header: ${f_inc}")
    aux_source_directory(${src} f_src)
    #message(STATUS "${group} source: ${f_src}")
    source_group(${group} FILES ${f_inc} ${f_src})

    # add to source list
    set(${f_list} ${${f_list}} ${f_inc} ${f_src} PARENT_SCOPE)

endfunction(add_group)

set(IEXP_SOURCE)
add_group(common ${ROOT_PATH}/include/common ${ROOT_PATH}/source/common IEXP_SOURCE)
add_group(math ${ROOT_PATH}/include/math ${ROOT_PATH}/source/math IEXP_SOURCE)
add_group(complex ${ROOT_PATH}/include/complex ${ROOT_PATH}/source/complex IEXP_SOURCE)
add_group(poly ${ROOT_PATH}/include/poly ${ROOT_PATH}/source/poly IEXP_SOURCE)
add_group(special ${ROOT_PATH}/include/special ${ROOT_PATH}/source/special IEXP_SOURCE)
add_group(combination ${ROOT_PATH}/include/combination ${ROOT_PATH}/source/combination IEXP_SOURCE)
add_group(multiset ${ROOT_PATH}/include/multiset ${ROOT_PATH}/source/multiset IEXP_SOURCE)
add_group(sort ${ROOT_PATH}/include/sort ${ROOT_PATH}/source/sort IEXP_SOURCE)
add_group(fft ${ROOT_PATH}/include/fft ${ROOT_PATH}/source/fft IEXP_SOURCE)
add_group(fft ${ROOT_PATH}/include/fft/fftw ${ROOT_PATH}/source/fft/fftw IEXP_SOURCE)
add_group(integral ${ROOT_PATH}/include/integral ${ROOT_PATH}/source/integral IEXP_SOURCE)
add_group(rand ${ROOT_PATH}/include/rand ${ROOT_PATH}/source/rand IEXP_SOURCE)
add_group(randist ${ROOT_PATH}/include/randist ${ROOT_PATH}/source/randist IEXP_SOURCE)
add_group(stats ${ROOT_PATH}/include/stats ${ROOT_PATH}/source/stats IEXP_SOURCE)
add_group(rstat ${ROOT_PATH}/include/rstat ${ROOT_PATH}/source/rstat IEXP_SOURCE)
add_group(histogram ${ROOT_PATH}/include/histogram ${ROOT_PATH}/source/histogram IEXP_SOURCE)
add_group(siman ${ROOT_PATH}/include/siman ${ROOT_PATH}/source/siman IEXP_SOURCE)
add_group(dae ${ROOT_PATH}/include/dae ${ROOT_PATH}/source/dae IEXP_SOURCE)

#
# build
#

add_library(imexpress STATIC ${IEXP_SOURCE})

#
# header file path
#

target_include_directories(imexpress PRIVATE ${ROOT_PATH}/include)

# gsl
target_include_directories(imexpress PRIVATE ${OUT_PATH}/gsl)

# eigen
target_include_directories(imexpress PRIVATE ${LIBRARY_PATH}/eigen)

# fft
target_include_directories(imexpress PRIVATE ${LIBRARY_PATH}/fft/fftw/api)

# dae
target_include_directories(imexpress PRIVATE ${LIBRARY_PATH}/dae/sundials/include)
target_include_directories(imexpress PRIVATE ${OUT_PATH}/dae/sundials/include)

#
# link
#

# gsl
check_exist(gsl_exist ${CMAKE_CURRENT_BINARY_DIR}/gsl "Debug;Release" "libgsl.a")
check_exist(gslcblas_exist ${CMAKE_CURRENT_BINARY_DIR}/gsl "Debug;Release" "libgslcblas.a")
if (gsl_exist AND gslcblas_exist)
    add_library(gsl STATIC IMPORTED GLOBAL)
    set_target_properties(gsl PROPERTIES IMPORTED_LOCATION ${gsl_exist})

    add_library(gslcblas STATIC IMPORTED GLOBAL)
    set_target_properties(gslcblas PROPERTIES IMPORTED_LOCATION ${gslcblas_exist})
else ()
    add_subdirectory(${LIBRARY_PATH}/gsl gsl)
endif ()
target_link_libraries(imexpress gsl gslcblas)

# fft
add_subdirectory(${LIBRARY_PATH}/fft fft)
if (FFT_FLOAT_LIB)
    target_link_libraries(imexpress ${FFT_FLOAT_LIB})
endif ()
target_link_libraries(imexpress ${FFT_LIB})

# dae
add_subdirectory(${LIBRARY_PATH}/dae dae)
target_link_libraries(imexpress ${DAE_LIB})

# test
add_subdirectory(${ROOT_PATH}/test test)

#
# package
#
