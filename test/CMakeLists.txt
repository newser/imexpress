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

set(test_src test_run.cpp)

aux_source_directory(case test_src)
message(STATUS "test source: ${test_src}")

#
# build
#

add_executable(test_run ${test_src})

#
# header file path
#

# gsl
target_include_directories(test_run PRIVATE ${OUT_PATH}/gsl)

# eigen
target_include_directories(test_run PRIVATE ${LIBRARY_PATH}/eigen)

# catch
target_include_directories(test_run PRIVATE ${ROOT_PATH}/test/catch/include)

#
# link
#

target_link_libraries(test_run imexpress)

#
# package
#
