/* Copyright (C) 2017 haniu (niuhao.cn@gmail.com)
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or (at
 * your option) any later version.
 *
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301,
 * USA.
 */

////////////////////////////////////////////////////////////
// import header files
////////////////////////////////////////////////////////////

#include <dae/error.h>

#include <cvode/cvode.h>
#include <sundials/sundials_linearsolver.h>

IEXP_NS_BEGIN

namespace dae {

////////////////////////////////////////////////////////////
// internal macro
////////////////////////////////////////////////////////////

#define throw_if(e, s)                                                         \
    case e:                                                                    \
        throw std::runtime_error(s);                                           \
        break

////////////////////////////////////////////////////////////
// internal type
////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////
// extern declaration
////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////
// global variant
////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////
// interface declaration
////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////
// interface implementation
////////////////////////////////////////////////////////////

void cv_throw(int e)
{
    if (e < 0) {
        switch (e) {
            throw_if(CV_TOO_MUCH_WORK, "CV_TOO_MUCH_WORK");
            throw_if(CV_TOO_MUCH_ACC, "CV_TOO_MUCH_ACC");
            throw_if(CV_ERR_FAILURE, "CV_ERR_FAILURE");
            throw_if(CV_CONV_FAILURE, "CV_CONV_FAILURE");

            throw_if(CV_LINIT_FAIL, "CV_LINIT_FAIL");
            throw_if(CV_LSETUP_FAIL, "CV_LSETUP_FAIL");
            throw_if(CV_LSOLVE_FAIL, "CV_TOO_MUCH_WORK");
            throw_if(CV_RHSFUNC_FAIL, "CV_TOO_MUCH_WORK");
            throw_if(CV_FIRST_RHSFUNC_ERR, "CV_LINIT_FAIL");
            throw_if(CV_REPTD_RHSFUNC_ERR, "CV_LSETUP_FAIL");
            throw_if(CV_UNREC_RHSFUNC_ERR, "CV_TOO_MUCH_WORK");
            throw_if(CV_RTFUNC_FAIL, "CV_TOO_MUCH_WORK");

            throw_if(CV_MEM_FAIL, "CV_LINIT_FAIL");
            throw_if(CV_MEM_NULL, "CV_LSETUP_FAIL");
            throw_if(CV_ILL_INPUT, "CV_TOO_MUCH_WORK");
            throw_if(CV_NO_MALLOC, "CV_TOO_MUCH_WORK");
            throw_if(CV_BAD_K, "CV_LINIT_FAIL");
            throw_if(CV_BAD_T, "CV_LSETUP_FAIL");
            throw_if(CV_BAD_DKY, "CV_TOO_MUCH_WORK");
            throw_if(CV_TOO_CLOSE, "CV_TOO_MUCH_WORK");

            default: {
                throw std::runtime_error("unknown cv error");
            } break;
        }
    }
}

void ls_throw(int e)
{
    if (e < 0) {
        switch (e) {
            throw_if(SUNLS_MEM_NULL,
                     "the memory argument to the function is NULL");
            throw_if(SUNLS_ILL_INPUT,
                     "an illegal input has been provided to the function");
            throw_if(SUNLS_MEM_FAIL, "failed memory access or allocation");
            throw_if(SUNLS_ATIMES_FAIL_UNREC,
                     "an unrecoverable failure occurred in the ATimes routine");
            throw_if(SUNLS_PSET_FAIL_UNREC,
                     "an unrecoverable failure occurred in the Pset routine");
            throw_if(SUNLS_PSOLVE_FAIL_UNREC,
                     "an unrecoverable failure occurred in the Psolve routine");
            throw_if(SUNLS_PACKAGE_FAIL_UNREC,
                     "an unrecoverable failure occurred in an external linear "
                     "base_linsol package");
            throw_if(SUNLS_GS_FAIL,
                     "a failure occurred during Gram-Schmidt "
                     "orthogonalization");
            throw_if(SUNLS_QRSOL_FAIL,
                     "a singular R matrix was encountered in a QR "
                     "factorization");
            default: {
                throw std::runtime_error("unknown cv error");
            } break;
        }
    }
}
}

IEXP_NS_END
