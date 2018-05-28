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

#ifndef __IEXP_DAE_DENSE_ODE__
#define __IEXP_DAE_DENSE_ODE__

////////////////////////////////////////////////////////////
// import header files
////////////////////////////////////////////////////////////

#include <common/common.h>

#include <dae/linsol.h>
#include <dae/ode.h>

#include <cvode/cvode_direct.h>

IEXP_NS_BEGIN

namespace dae {

////////////////////////////////////////////////////////////
// macro definition
////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////
// type definition
////////////////////////////////////////////////////////////

class dense_ode : public ode<dense_ode>
{
  public:
    template <typename T>
    dense_ode(multistep lmm,
              const DY_TYPE &dy,
              double t0,
              const DenseBase<T> &y0)
        : ode<dense_ode>(lmm, iteration::NEWTON, dy, t0, y0)
        , m_A(SUNDenseMatrix(y0.size(), y0.size()))
        , m_ls(m_y0.n_vector(), m_A)
    {
        IEXP_NOT_NULLPTR(m_A);
    }

    dense_linsol &linsol()
    {
        return m_ls;
    }

    dense_ode &jac(const JAC_TYPE &jac)
    {
        m_jac = jac;
        return *this;
    }

    void prepare()
    {
        cv_check(CVDlsSetLinearSolver(m_cvode, m_ls.sunls(), m_A));

        if (m_jac) {
            cv_check(CVDlsSetJacFn(m_cvode, jac_func<dense_ode>::s_jac));
        } else {
            cv_check(CVDlsSetJacFn(m_cvode, nullptr));
        }
    }

    int compute_jac(double t,
                    Map<const VectorXd> &y,
                    Map<const VectorXd> &fy,
                    Map<MatrixXd> &jac)
    {
        eigen_assert(m_jac);
        return m_jac(t, y, fy, jac);
    }

  private:
    SUNMatrix m_A;
    dense_linsol m_ls;
    JAC_TYPE m_jac;
};

////////////////////////////////////////////////////////////
// global variants
////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////
// interface declaration
////////////////////////////////////////////////////////////
}

IEXP_NS_END

#endif /* __IEXP_DAE_DENSE_ODE__ */
