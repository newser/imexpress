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

#
# definition
#

#
# source file
#

set(test_src test_imexpress.cpp)

aux_source_directory(case test_src)
#message(STATUS "test source: ${test_src}")

#
# build
#

add_executable(test_imexpress ${test_src})

#
# header file path
#

target_include_directories(test_imexpress PRIVATE ${ROOT_PATH}/include)
target_include_directories(test_imexpress PRIVATE ${ROOT_PATH}/test)

# gsl
target_include_directories(test_imexpress PRIVATE ${OUT_PATH}/gsl)

# eigen
target_include_directories(test_imexpress PRIVATE ${LIBRARY_PATH}/eigen)

# fft
target_include_directories(test_imexpress PRIVATE ${LIBRARY_PATH}/fft/fftw/api)

# dae
target_include_directories(test_imexpress PRIVATE ${LIBRARY_PATH}/dae/sundials/include)
target_include_directories(test_imexpress PRIVATE ${OUT_PATH}/dae/sundials/include)

# catch
target_include_directories(test_imexpress PRIVATE ${LIBRARY_PATH}/Catch2/single_include)

# plot
find_package(MathGL2 QUIET)
if (MATHGL2_FOUND)
    target_include_directories(test_imexpress PRIVATE ${MATHGL2_INCLUDE_DIRS})
    target_compile_definitions(test_imexpress PRIVATE -DIEXP_MGL2)
endif ()

#
# link
#

target_link_libraries(test_imexpress imexpress)

if (MATHGL2_FOUND)
    target_link_libraries(test_imexpress ${MATHGL2_LIBRARIES})
endif ()

#
# package
#
