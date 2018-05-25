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
    template <typename T>
    solver()
        : m_ls(nullptr)
    {
    }

    SUNLinearSolver m_ls;
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
        m_ls = SUNDenseLinearSolver(sun_y.sunvec(), sun_a.sunmat());
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
