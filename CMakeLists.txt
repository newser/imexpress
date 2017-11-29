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
add_group(polynomial ${ROOT_PATH}/include/polynomial ${ROOT_PATH}/source/polynomial IEXP_SOURCE)
add_group(special ${ROOT_PATH}/include/special ${ROOT_PATH}/source/special IEXP_SOURCE)
add_group(combination ${ROOT_PATH}/include/combination ${ROOT_PATH}/source/combination IEXP_SOURCE)
add_group(multiset ${ROOT_PATH}/include/multiset ${ROOT_PATH}/source/multiset IEXP_SOURCE)
add_group(sort ${ROOT_PATH}/include/sort ${ROOT_PATH}/source/sort IEXP_SOURCE)

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

#
# link
#

# gsl
add_subdirectory(${LIBRARY_PATH}/gsl gsl)
target_link_libraries(imexpress gsl)

# fft
add_subdirectory(${LIBRARY_PATH}/fft fft)
target_link_libraries(imexpress ${FFT_LIB})

# test
add_subdirectory(${ROOT_PATH}/test test)

#
# package
#
