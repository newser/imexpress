# Licensed to the Apache Software Foundation (ASF) under one or more
# contributor license agreements.  See the NOTICE file distributed with
# this work for additional information regarding copyright ownership.
# The ASF licenses this file to You under the Apache License, Version 2.0
# (the "License"); you may not use this file except in compliance with
# the License.  You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

#
# definition
#

set(FFTW_PATH ${CMAKE_CURRENT_SOURCE_DIR}/fftw)

#
# customization
#

hide(BUILD_SHARED_LIBS BOOL OFF)
hide(BUILD_TESTS BOOL OFF)

hide(ENABLE_OPENMP BOOL OFF) #todo
hide(ENABLE_THREADS BOOL OFF) #todo
hide(WITH_COMBINED_THREADS BOOL OFF) #todo

hide(ENABLE_FLOAT BOOL ON)
hide(ENABLE_LONG_DOUBLE BOOL OFF)
hide(ENABLE_QUAD_PRECISION BOOL OFF)

hide(ENABLE_SSE BOOL ON)
hide(ENABLE_SSE2 BOOL ON)
hide(ENABLE_AVX BOOL ON)
hide(ENABLE_AVX2 BOOL ON)

add_subdirectory(${FFTW_PATH})

if (ENABLE_QUAD_PRECISION)
    set(fft_name fftw3q)
elseif (ENABLE_LONG_DOUBLE)
    set(fft_name fftw3l)
elseif (ENABLE_FLOAT)
    set(fft_name fftw3f)
endif ()

# no need to show fftw compile warnings
target_compile_options(${fft_name} PRIVATE ${FLAG_DISABLE_WARNING})

#
# result
#

set(FFT_LIB ${fft_name} PARENT_SCOPE)
