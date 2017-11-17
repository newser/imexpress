# Licensed to the Apache Software Foundation (ASF) under one or more
# contributor license agreements.    See the NOTICE file distributed with
# this work for additional information regarding copyright ownership.
# The ASF licenses this file to You under the Apache License, Version 2.0
# (the "License"); you may not use this file except in compliance with
# the License.    You may obtain a copy of the License at
#
#         http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

#
# load environment info
#

if (CMAKE_SYSTEM_NAME STREQUAL Android)
    if (NOT ANDROID_STUDIO)
        # cmake support ndk build from 3.7 but android studio may not have latest cmake
        cmake_minimum_required(VERSION 3.7.2)
    endif ()

    message(STATUS "CMAKE_SYSTEM_NAME: ${CMAKE_SYSTEM_NAME}")
    message(STATUS "CMAKE_SYSTEM_VERSION: ${CMAKE_SYSTEM_VERSION}")
    message(STATUS "CMAKE_ANDROID_ARCH_ABI: ${CMAKE_ANDROID_ARCH_ABI}")
    message(STATUS "CMAKE_ANDROID_NDK: ${CMAKE_ANDROID_NDK}")

elseif (CMAKE_HOST_SYSTEM_NAME STREQUAL Darwin)
    # available only on mac
    set(IOS 0 CACHE BOOL "use predefined ios-xcode toolchain")

    if (IOS)
        include(${BUILD_PATH}/ios-xcode.toolchain.cmake)
    endif ()

endif ()

# get application os
function(get_os os)
    if (CMAKE_SYSTEM_NAME STREQUAL Windows)
        set(${os} windows PARENT_SCOPE)
    elseif (CMAKE_SYSTEM_NAME STREQUAL Linux)
        set(${os} linux PARENT_SCOPE)
    elseif (CMAKE_SYSTEM_NAME STREQUAL Darwin)
        set(${os} macos PARENT_SCOPE)
    elseif (CMAKE_SYSTEM_NAME STREQUAL iOS)
        set(${os} ios PARENT_SCOPE)
    elseif (CMAKE_SYSTEM_NAME STREQUAL Android)
        set(${os} android PARENT_SCOPE)
    else ()
        # more to be added
        message(FATAL_ERROR "unsupported os: ${${CMAKE_SYSTEM_NAME}}")
    endif ()
endfunction(get_os)

# get application cpu
function(get_cpu cpu)
    set(__cpu ${CMAKE_SYSTEM_PROCESSOR})
    string(TOLOWER ${__cpu} __cpu)
    # CMAKE_SYSTEM_PROCESSOR is exactly the cpu even in cross build mode
    
    if (__cpu MATCHES "(.*intel.*|.*amd.*|.*x64.*|.*x86.*|.*i.86)")
        set(${cpu} x86 PARENT_SCOPE)
    elseif (__cpu MATCHES "(.*arm.*)" OR __cpu MATCHES "(.*aarch64.*)")
        set(${cpu} arm PARENT_SCOPE)
    elseif (__cpu MATCHES "(.*mips.*)")
        set(${cpu} mips PARENT_SCOPE)
    else ()
        # more to be add
        message(FATAL_ERROR "unsupported cpu: ${__cpu}")
    endif ()

endfunction(get_cpu)

# get application toolchain
function(get_toolchain toolchain)
    if (CMAKE_SYSTEM_NAME STREQUAL Android)
        set(${toolchain} ndk PARENT_SCOPE)
    elseif (CMAKE_GENERATOR MATCHES ".*Visual Studio.*")
        # vs must use msvc as compiler
        set(${toolchain} msvc PARENT_SCOPE)
    elseif (CMAKE_GENERATOR MATCHES ".*Unix Makefiles.*" AND 
            CMAKE_C_COMPILER_ID MATCHES GNU AND
            CMAKE_CXX_COMPILER_ID MATCHES GNU)
        set(${toolchain} gnu PARENT_SCOPE)
    elseif (CMAKE_GENERATOR MATCHES ".*Xcode.*" AND
            CMAKE_C_COMPILER_ID MATCHES Clang AND
            CMAKE_CXX_COMPILER_ID MATCHES Clang)
        set(${toolchain} xcode PARENT_SCOPE)
    else ()
        message(FATAL_ERROR "unsupported toolchain: ${CMAKE_GENERATOR}, ${CMAKE_C_COMPILER_ID}")
    endif ()
endfunction(get_toolchain)

# construct environment name
message(STATUS "DETECTING ENVIRONMENT ...")

# os
message(STATUS "detecting os")
get_os(ENV_OS)
set(ENV_OS ${ENV_OS} CACHE INTERNAL "os")
message(STATUS "os: ${ENV_OS}")

# cpu
message(STATUS "detecting cpu")
get_cpu(ENV_CPU)
if (IOS_SIMULATOR)
    set(ENV_CPU x86)
endif ()
set(ENV_CPU ${ENV_CPU} CACHE INTERNAL "cpu")
message(STATUS "cpu: ${ENV_CPU}")

# toolchain
message(STATUS "detecting toochain")
get_toolchain(ENV_TOOLCHAIN)
set(ENV_TOOLCHAIN ${ENV_TOOLCHAIN} CACHE INTERNAL "toochain")
message(STATUS "toochain: ${ENV_TOOLCHAIN}")

set(ENV "${ENV_OS}-${ENV_TOOLCHAIN}")

include(TestBigEndian)
TEST_BIG_ENDIAN(ENV_BIG_ENDIAN)

message(STATUS "DETECTING ENVIRONMENT DONE: ${ENV}")

#
# configuration
#

# generate header file
configure_file(${ROOT_PATH}/include/common/environment.h.in
               ${ROOT_PATH}/include/common/environment.h
               @ONLY)
