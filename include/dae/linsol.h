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

#ifndef __IEXP_DAE_LINSOL__
#define __IEXP_DAE_LINSOL__

////////////////////////////////////////////////////////////
// import header files
////////////////////////////////////////////////////////////

#include <common/common.h>

#include <dae/sunmat.h>
#include <dae/sunvec.h>

#include <sundials/sundials_linearsolver.h>
#include <sunlinsol/sunlinsol_dense.h>

IEXP_NS_BEGIN

namespace dae {

namespace linsol {
////////////////////////////////////////////////////////////
// macro definition
////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////
// type definition
////////////////////////////////////////////////////////////

// ========================================
// base
// ========================================

template <typename Derived>
class solver
{
    DEFINE_DERIVED

  public:
    SUNLinearSolver sunls()
    {
        return m_ls;
    }

  protected:
    static const char *s_err_desc[9];

    solver()
        : m_ls(nullptr)
    {
    }

    void check(int e)
    {
        if (e < 0) {
            throw std::runtime_error(s_err_desc[-e - 1]);
        }
    }

    SUNLinearSolver m_ls;
};

template <typename Derived>
const char *solver<Derived>::s_err_desc[9] = {
    "the memory argument to the function is NULL",
    "an illegal input has been provided to the function",
    "failed memory access or allocation",
    "an unrecoverable failure occurred in the ATimes routine",
    "an unrecoverable failure occurred in the Pset routine",
    "an unrecoverable failure occurred in the Psolve routine",
    "an unrecoverable failure occurred in an external linear solver package",
    "a failure occurred during Gram-Schmidt orthogonalization",
    "a singular R matrix was encountered in a QR factorization",
};

// ========================================
// dense
// ========================================

class dense : public solver<dense>
{
  public:
    template <typename T, typename U>
    dense(const DenseBase<T> &y, const DenseBase<T> &a)
    {
        sunvec_serial sun_y(y, false);
        sunmat_dense sun_a(a, false);
        m_ls = SUNDenseLinearSolver(sun_y.n_vector(), sun_a.sunmatrix());
        IEXP_NOT_NULLPTR(m_ls);
    }
};

////////////////////////////////////////////////////////////
// global variants
////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////
// interface declaration
////////////////////////////////////////////////////////////
}
}

IEXP_NS_END

#endif /* __IEXP_DAE_LINSOL__ */