# Copyright (C) 2017 haniu (niuhao.cn@gmail.com)
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.

project(imexpress)

cmake_minimum_required(VERSION 3.1.3)

#
# definition
#

set(ROOT_PATH ${CMAKE_CURRENT_LIST_DIR})
set(OUT_PATH ${CMAKE_CURRENT_BINARY_DIR})

set(BUILD_PATH ${ROOT_PATH}/build)
set(LIBRARY_PATH ${ROOT_PATH}/library)

#
# source file
#

function(add_group group inc src f_list)
    # make group
    file(GLOB f_inc "${inc}/*.h")
	#message(STATUS "${group} header: ${f_inc}")	
    aux_source_directory(${src} f_src)
	message(STATUS "${group} source: ${f_src}")
    source_group(${group} FILES ${f_inc} ${f_src})

    # add to source list
    set(${f_list} ${${f_list}} ${f_inc} ${f_src} PARENT_SCOPE)

endfunction(add_group)

set(IEP_SOURCE)

add_group(common 
          ${ROOT_PATH}/include/common
          ${ROOT_PATH}/source/common
          IEP_SOURCE)

#
# build
#

add_library(imexpress STATIC ${IEP_SOURCE})

#
# header file path
#

# gsl
target_include_directories(imexpress PRIVATE ${OUT_PATH}/gsl)

# eigen
target_include_directories(imexpress PRIVATE ${LIBRARY_PATH}/eigen)

#
# link
#

# gsl
add_subdirectory(${LIBRARY_PATH}/gsl gsl)
target_link_libraries(imexpress gsl)

# test
add_subdirectory(${ROOT_PATH}/test test)

#
# package
#
