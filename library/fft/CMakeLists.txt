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

hide(ENABLE_LONG_DOUBLE BOOL OFF)
hide(ENABLE_QUAD_PRECISION BOOL OFF)

hide(ENABLE_SSE BOOL ON)
hide(ENABLE_SSE2 BOOL ON)
hide(ENABLE_AVX BOOL ON)
hide(ENABLE_AVX2 BOOL ON)

# todo: if enabled float
if (1)
    foreach(t "Debug" "Release")
        foreach(n "libfftw3f.a")
            if (EXISTS ${CMAKE_CURRENT_BINARY_DIR}/fftw/float/${t}/${n})
                set(fft3wf_path ${CMAKE_CURRENT_BINARY_DIR}/fftw/float/${t}/${n})
                break()
            endif ()
        endforeach ()
    endforeach ()

    if (fft3wf_path)
        add_library(fftw3f STATIC IMPORTED GLOBAL)
        set_target_properties(fftw3f PROPERTIES IMPORTED_LOCATION ${fft3wf_path})
    else ()
        hide(ENABLE_FLOAT BOOL ON)
        add_subdirectory(${FFTW_PATH} ${CMAKE_CURRENT_BINARY_DIR}/fftw/float)

        # no need to show fftw compile warnings
        target_compile_options(fftw3f PRIVATE ${FLAG_DISABLE_WARNING})
    endif ()

    set(FFT_FLOAT_LIB fftw3f PARENT_SCOPE)
endif ()

# todo: if enabled double
if (1)
    foreach(t "Debug" "Release")
        foreach(n "libfftw3.a")
            if (EXISTS ${CMAKE_CURRENT_BINARY_DIR}/fftw/double/${t}/${n})
                set(fft3w_path ${CMAKE_CURRENT_BINARY_DIR}/fftw/double/${t}/${n})
                break()
            endif ()
        endforeach ()
    endforeach ()

    if (fft3w_path)
        add_library(fftw3 STATIC IMPORTED GLOBAL)
        set_target_properties(fftw3 PROPERTIES IMPORTED_LOCATION ${fft3w_path})
    else ()
        hide(ENABLE_FLOAT BOOL OFF)
        add_subdirectory(${FFTW_PATH} ${CMAKE_CURRENT_BINARY_DIR}/fftw/double)

        # no need to show fftw compile warnings
        target_compile_options(fftw3 PRIVATE ${FLAG_DISABLE_WARNING})
    endif ()

    set(FFT_LIB fftw3 PARENT_SCOPE)
endif ()

#
# result
#
