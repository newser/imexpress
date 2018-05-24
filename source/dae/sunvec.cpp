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

#include <dae/sunvec.h>

IEXP_NS_BEGIN

namespace dae {

////////////////////////////////////////////////////////////
// internal macro
////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////
// internal type
////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////
// extern declaration
////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////
// global variant
////////////////////////////////////////////////////////////

struct _generic_N_Vector_Ops sunvec_ops = {
    N_VGetVectorID_Serial,
    N_VClone_Serial,
    N_VCloneEmpty_Serial,
    N_VDestroy_Serial,
    N_VSpace_Serial,
    N_VGetArrayPointer_Serial,
    N_VSetArrayPointer_Serial,
    N_VLinearSum_Serial,
    N_VConst_Serial,
    N_VProd_Serial,
    N_VDiv_Serial,
    N_VScale_Serial,
    N_VAbs_Serial,
    N_VInv_Serial,
    N_VAddConst_Serial,
    N_VDotProd_Serial,
    N_VMaxNorm_Serial,
    N_VWrmsNorm_Serial,
    N_VWrmsNormMask_Serial,
    N_VMin_Serial,
    N_VWL2Norm_Serial,
    N_VL1Norm_Serial,
    N_VCompare_Serial,
    N_VInvTest_Serial,
    N_VConstrMask_Serial,
    N_VMinQuotient_Serial,
};

////////////////////////////////////////////////////////////
// interface declaration
////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////
// interface implementation
////////////////////////////////////////////////////////////
}

IEXP_NS_END
