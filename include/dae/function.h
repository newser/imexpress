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

#ifndef __IEXP_DAE_FUNCTION__
#define __IEXP_DAE_FUNCTION__

////////////////////////////////////////////////////////////
// import header files
////////////////////////////////////////////////////////////

#include <common/common.h>

#include <cvode/cvode.h>
#include <nvector/nvector_serial.h>
#include <sundials/sundials_matrix.h>
#include <sunmatrix/sunmatrix_dense.h>

IEXP_NS_BEGIN

namespace dae {

////////////////////////////////////////////////////////////
// macro definition
////////////////////////////////////////////////////////////

#define DY_TYPE typename dy_func<void>::type

#define JAC_TYPE typename jac_func<void>::type

#define WEIGHT_TYPE typename weight_func<void>::type

#define PSETUP_TYPE typename psetup_func<void>::type

#define PSOLVE_TYPE typename psolve_func<void>::type

////////////////////////////////////////////////////////////
// type definition
////////////////////////////////////////////////////////////

template <typename T>
class dy_func
{
  public:
    using type = std::function<
        int(double t, Map<const VectorXd> &y, Map<VectorXd> &dy, void *opaque)>;

    static int s_dy(realtype t, N_Vector y, N_Vector ydot, void *user_data)
    {
        static_assert(TYPE_IS(realtype, double), "only support double scalar");

        eigen_assert(N_VGetVectorID(y) == SUNDIALS_NVEC_SERIAL);
        eigen_assert(N_VGetVectorID(ydot) == SUNDIALS_NVEC_SERIAL);
        eigen_assert(N_VGetLength_Serial(y) == N_VGetLength_Serial(ydot));

        Map<const VectorXd> mapped_y(N_VGetArrayPointer(y),
                                     N_VGetLength_Serial(y));
        Map<VectorXd> mapped_dy(N_VGetArrayPointer(ydot),
                                N_VGetLength_Serial(ydot));
        return ((T *)user_data)
            ->compute_dy(t, mapped_y, mapped_dy, ((T *)user_data)->opaque());
    }
};

template <typename T>
class jac_func
{
  public:
    using type = std::function<int(double t,
                                   Map<const VectorXd> &y,
                                   Map<const VectorXd> &fy,
                                   Map<MatrixXd> &jac,
                                   void *opaque)>;

    static int s_jac(realtype t,
                     N_Vector y,
                     N_Vector fy,
                     SUNMatrix jac,
                     void *user_data,
                     N_Vector tmp1,
                     N_Vector tmp2,
                     N_Vector tmp3)
    {
        static_assert(TYPE_IS(realtype, double), "only support double scalar");

        eigen_assert(N_VGetVectorID(y) == SUNDIALS_NVEC_SERIAL);
        eigen_assert(N_VGetVectorID(fy) == SUNDIALS_NVEC_SERIAL);
        eigen_assert(N_VGetLength_Serial(y) == N_VGetLength_Serial(fy));

        eigen_assert(SUNMatGetID(jac) == SUNMATRIX_DENSE);

        Map<const VectorXd> mapped_y(N_VGetArrayPointer(y),
                                     N_VGetLength_Serial(y));
        Map<const VectorXd> mapped_fy(N_VGetArrayPointer(fy),
                                      N_VGetLength_Serial(fy));
        Map<MatrixXd> mapped_jac(SM_DATA_D(jac),
                                 SM_ROWS_D(jac),
                                 SM_COLUMNS_D(jac));
        return ((T *)user_data)
            ->compute_jac(t,
                          mapped_y,
                          mapped_fy,
                          mapped_jac,
                          ((T *)user_data)->opaque());
    }
};

template <typename T>
class weight_func
{
  public:
    using type = std::function<
        int(Map<const VectorXd> &y, Map<VectorXd> &ewt, void *opaque)>;

    static int s_weight(N_Vector y, N_Vector ewt, void *user_data)
    {
        eigen_assert(N_VGetVectorID(y) == SUNDIALS_NVEC_SERIAL);
        eigen_assert(N_VGetVectorID(ewt) == SUNDIALS_NVEC_SERIAL);

        Map<const VectorXd> mapped_y(N_VGetArrayPointer(y),
                                     N_VGetLength_Serial(y));
        Map<VectorXd> mapped_ewt(N_VGetArrayPointer(ewt),
                                 N_VGetLength_Serial(ewt));
        return ((T *)user_data)
            ->compute_weight(mapped_y, mapped_ewt, ((T *)user_data)->opaque());
    }
};

template <typename T>
class psetup_func
{
  public:
    using type = std::function<int(double t,
                                   Map<const VectorXd> &y,
                                   Map<const VectorXd> &dy,
                                   bool jac_ok,
                                   bool &jac_updated,
                                   double gamma,
                                   void *opaque)>;

    static int s_psetup(realtype t,
                        N_Vector y,
                        N_Vector fy,
                        booleantype jok,
                        booleantype *jcurPtr,
                        realtype gamma,
                        void *user_data)
    {
        static_assert(TYPE_IS(realtype, double), "only support double scalar");

        eigen_assert(N_VGetVectorID(y) == SUNDIALS_NVEC_SERIAL);
        eigen_assert(N_VGetVectorID(fy) == SUNDIALS_NVEC_SERIAL);
        eigen_assert(N_VGetLength_Serial(y) == N_VGetLength_Serial(fy));

        Map<const VectorXd> mapped_y(N_VGetArrayPointer(y),
                                     N_VGetLength_Serial(y));
        Map<const VectorXd> mapped_dy(N_VGetArrayPointer(fy),
                                      N_VGetLength_Serial(fy));
        bool jac_updated;
        int ret = ((T *)user_data)
                      ->psetup(t,
                               mapped_y,
                               mapped_dy,
                               jok,
                               jac_updated,
                               gamma,
                               ((T *)user_data)->opaque());
        *jcurPtr = jac_updated;
        return ret;
    }
};

template <typename T>
class psolve_func
{
  public:
    using type = std::function<int(double t,
                                   Map<const VectorXd> &y,
                                   Map<const VectorXd> &dy,
                                   Map<const VectorXd> &r,
                                   Map<VectorXd> &z,
                                   double gamma,
                                   double delta,
                                   precondition pretype,
                                   void *opaque)>;

    static int s_psolve(realtype t,
                        N_Vector y,
                        N_Vector fy,
                        N_Vector r,
                        N_Vector z,
                        realtype gamma,
                        realtype delta,
                        int lr,
                        void *user_data)
    {
        static_assert(TYPE_IS(realtype, double), "only support double scalar");

        eigen_assert(N_VGetVectorID(y) == SUNDIALS_NVEC_SERIAL);
        eigen_assert(N_VGetVectorID(fy) == SUNDIALS_NVEC_SERIAL);
        eigen_assert(N_VGetVectorID(r) == SUNDIALS_NVEC_SERIAL);
        eigen_assert(N_VGetVectorID(z) == SUNDIALS_NVEC_SERIAL);
        eigen_assert(N_VGetLength_Serial(y) == N_VGetLength_Serial(fy));
        eigen_assert(N_VGetLength_Serial(y) == N_VGetLength_Serial(r));
        eigen_assert(N_VGetLength_Serial(y) == N_VGetLength_Serial(z));

        Map<const VectorXd> mapped_y(N_VGetArrayPointer(y),
                                     N_VGetLength_Serial(y));
        Map<const VectorXd> mapped_dy(N_VGetArrayPointer(fy),
                                      N_VGetLength_Serial(fy));
        Map<const VectorXd> mapped_r(N_VGetArrayPointer(r),
                                     N_VGetLength_Serial(r));
        Map<VectorXd> mapped_z(N_VGetArrayPointer(z), N_VGetLength_Serial(z));
        return ((T *)user_data)
            ->psolve(t,
                     mapped_y,
                     mapped_dy,
                     mapped_r,
                     mapped_z,
                     gamma,
                     delta,
                     lr == 1 ? precondition::LEFT : precondition::RIGHT,
                     ((T *)user_data)->opaque());
    }
};

////////////////////////////////////////////////////////////
// global variants
////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////
// interface declaration
////////////////////////////////////////////////////////////
}

IEXP_NS_END

#endif /* __IEXP_DAE_FUNCTION__ */
