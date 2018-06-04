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

#ifndef __IEXP_DAE_KRYLOV_ODE__
#define __IEXP_DAE_KRYLOV_ODE__

////////////////////////////////////////////////////////////
// import header files
////////////////////////////////////////////////////////////

#include <common/common.h>

#include <dae/linsol.h>
#include <dae/ode.h>

#include <cvode/cvode_bandpre.h>
#include <cvode/cvode_spils.h>

IEXP_NS_BEGIN

namespace dae {

////////////////////////////////////////////////////////////
// macro definition
////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////
// type definition
////////////////////////////////////////////////////////////

template <typename LS>
class krylov_ode : public ode<krylov_ode<LS>>
{
  public:
    template <typename T>
    krylov_ode(multistep lmm,
               const DY_TYPE &dy,
               double t0,
               const DenseBase<T> &y0,
               int max_kdim = 0)
        : ode<krylov_ode<LS>>(lmm, iteration::NEWTON, dy, t0, y0)
        , m_ls(this->m_y0.n_vector(), precondition::NONE, max_kdim)
        , m_dim(0)
        , m_up_dim(0)
        , m_low_dim(0)
    {
    }

    template <typename T>
    krylov_ode(multistep lmm,
               const DY_TYPE &dy,
               double t0,
               const DenseBase<T> &y0,
               precondition pretype,
               int dim,
               int up_dim,
               int low_dim,
               int max_kdim = 0)
        : ode<krylov_ode<LS>>(lmm, iteration::NEWTON, dy, t0, y0)
        , m_ls(this->m_y0.n_vector(), pretype, max_kdim)
        , m_dim(dim)
        , m_up_dim(up_dim)
        , m_low_dim(low_dim)
    {
    }

    LS &linsol()
    {
        return m_ls;
    }

    void prepare()
    {
        cv_check(CVSpilsSetLinearSolver(this->m_cvode, m_ls.sunls()));

        if (m_dim != 0) {
            cv_check(CVBandPrecInit(this->m_cvode, m_dim, m_up_dim, m_low_dim));
        }
    }

  protected:
    LS m_ls;
    sunindextype m_dim, m_up_dim, m_low_dim;
};

using spgmr_ode = krylov_ode<spgmr_linsol>;

using spfgmr_ode = krylov_ode<spfgmr_linsol>;

using spbcgs_ode = krylov_ode<spbcgs_linsol>;

using sptfqmr_ode = krylov_ode<sptfqmr_linsol>;

////////////////////////////////////////////////////////////
// global variants
////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////
// interface declaration
////////////////////////////////////////////////////////////
}

IEXP_NS_END

#endif /* __IEXP_DAE_KRYLOV_ODE__ */
