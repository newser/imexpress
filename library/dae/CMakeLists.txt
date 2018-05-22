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

set(DAE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/sundials)

#
# customization
#

hide(BUILD_ARKODE BOOL ON)
hide(BUILD_CVODE BOOL ON)
hide(BUILD_CVODES BOOL ON)
hide(BUILD_IDA BOOL ON)
hide(BUILD_IDAS BOOL ON)
hide(BUILD_KINSOL BOOL ON)

hide(USE_XSDK_DEFAULTS BOOL OFF)

hide(SUNDIALS_PRECISION STRING double)

#if (CMAKE_SIZEOF_VOID_P EQUAL 8)
#    hide(SUNDIALS_INDEX_TYPE STRING "int64_t")
#else()
#    hide(SUNDIALS_INDEX_TYPE STRING "int32_t")
#endif()
hide(SUNDIALS_INDEX_TYPE STRING "int32_t")

hide(FCMIX_ENABLE BOOL OFF)

hide(BUILD_STATIC_LIBS BOOL ON)
hide(BUILD_SHARED_LIBS BOOL OFF)

hide(USE_GENERIC_MATH BOOL ON)

# may enable some of below in future
hide(MPI_ENABLE BOOL OFF)
hide(OPENMP_ENABLE BOOL OFF)
hide(PTHREAD_ENABLE BOOL OFF)
hide(CUDA_ENABLE BOOL OFF)
hide(RAJA_ENABLE BOOL OFF)
hide(BLAS_ENABLE BOOL OFF)
hide(LAPACK_ENABLE BOOL OFF)
hide(SUPERLUMT_ENABLE BOOL OFF)
hide(KLU_ENABLE BOOL OFF)
hide(HYPRE_ENABLE BOOL OFF)
hide(PETSC_ENABLE BOOL OFF)

hide(EXAMPLES_ENABLE_C BOOL OFF)
hide(EXAMPLES_ENABLE_CXX BOOL OFF)

set(DAE_LIB
    sundials_nvecserial_static
    sundials_sunmatrixdense_static
    sundials_sunmatrixband_static
    sundials_sunmatrixsparse_static
    sundials_sunlinsoldense_static
    sundials_sunlinsolband_static
    sundials_sunlinsolspgmr_static
    sundials_sunlinsolspfgmr_static
    sundials_sunlinsolspbcgs_static
    sundials_sunlinsolsptfqmr_static
    sundials_arkode_static
    sundials_cvode_static
    sundials_cvodes_static
    sundials_ida_static
    sundials_idas_static
    sundials_kinsol_static
    PARENT_SCOPE)

set(build_sundials TRUE)
set(mode)
foreach(t "Debug" "Release")
    if (EXISTS ${CMAKE_CURRENT_BINARY_DIR}/sundials/src/nvec_ser/${t}/libsundials_nvecserial.a)
        # we assume all libs are built if one of them is found
        set(build_sundials FALSE)
        set(mode ${t})
        break()
    endif ()
endforeach ()

macro(import_sundials name dir)
    add_library(${name}_static STATIC IMPORTED GLOBAL)
    set_target_properties(${name}_static PROPERTIES IMPORTED_LOCATION ${CMAKE_CURRENT_BINARY_DIR}/sundials/src/${dir}/${mode}/lib${name}.a)
endmacro()

if (build_sundials)
    add_subdirectory(${DAE_PATH} ${CMAKE_CURRENT_BINARY_DIR}/sundials)

    # no need to show compile warnings
    foreach (lib ${DAE_LIB})
        target_compile_options(${lib} PRIVATE ${FLAG_DISABLE_WARNING})
    endforeach()
else ()
    import_sundials(sundials_nvecserial nvec_ser)
    import_sundials(sundials_sunmatrixdense sunmat_dense)
    import_sundials(sundials_sunmatrixband sunmat_band)
    import_sundials(sundials_sunmatrixsparse sunmat_sparse)
    import_sundials(sundials_sunlinsoldense sunlinsol_dense)
    import_sundials(sundials_sunlinsolband sunlinsol_band)
    import_sundials(sundials_sunlinsolspgmr sunlinsol_spgmr)
    import_sundials(sundials_sunlinsolspfgmr sunlinsol_spfgmr)
    import_sundials(sundials_sunlinsolspbcgs sunlinsol_spbcgs)
    import_sundials(sundials_sunlinsolsptfqmr sunlinsol_sptfqmr)
    import_sundials(sundials_arkode arkode)
    import_sundials(sundials_cvode cvode)
    import_sundials(sundials_cvodes cvodes)
    import_sundials(sundials_ida ida)
    import_sundials(sundials_idas idas)
    import_sundials(sundials_kinsol kinsol)
endif ()

#
# result
#
