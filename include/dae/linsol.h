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

#include <dae/def.h>
#include <dae/sunmat.h>
#include <dae/sunvec.h>

#include <sundials/sundials_linearsolver.h>
#include <sunlinsol/sunlinsol_band.h>
#include <sunlinsol/sunlinsol_dense.h>
#include <sunlinsol/sunlinsol_spbcgs.h>
#include <sunlinsol/sunlinsol_spfgmr.h>
#include <sunlinsol/sunlinsol_spgmr.h>
#include <sunlinsol/sunlinsol_sptfqmr.h>

IEXP_NS_BEGIN

namespace dae {

////////////////////////////////////////////////////////////
// macro definition
////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////
// type definition
////////////////////////////////////////////////////////////

// ========================================
// base linear solver
// ========================================

template <typename Derived>
class linsol
{
    DEFINE_DERIVED

  public:
    SUNLinearSolver sunls()
    {
        return m_ls;
    }

  protected:
    linsol()
        : m_ls(nullptr)
    {
    }

    SUNLinearSolver m_ls;
};

// ========================================
// dense linear solver
// ========================================

class dense_linsol : public linsol<dense_linsol>
{
  public:
    dense_linsol(N_Vector y, SUNMatrix a)
    {
        m_ls = SUNDenseLinearSolver(y, a);
        IEXP_NOT_NULLPTR(m_ls);
    }
};

// ========================================
// band linear solver
// ========================================

class band_linsol : public linsol<band_linsol>
{
  public:
    band_linsol(N_Vector y, SUNMatrix a)
    {
        m_ls = SUNBandLinearSolver(y, a);
        IEXP_NOT_NULLPTR(m_ls);
    }
};

// ========================================
// Scaled, Preconditioned, Generalized Minimum Residual
// ========================================

class spgmr_linsol : public linsol<spgmr_linsol>
{
  public:
    spgmr_linsol(N_Vector y, precondition pretype, int maxl)
    {
        m_ls = SUNSPGMR(y, (int)pretype, maxl);
        IEXP_NOT_NULLPTR(m_ls);
    }

    spgmr_linsol &precondition(precondition pretype)
    {
        SUNSPGMRSetPrecType(m_ls, (int)pretype);
        return *this;
    }

    spgmr_linsol &gram_schmidt(gram_schmidt gstype)
    {
        SUNSPGMRSetGSType(m_ls, (int)gstype);
        return *this;
    }

    spgmr_linsol &max_restart(int maxrs)
    {
        SUNSPGMRSetMaxRestarts(m_ls, maxrs);
        return *this;
    }
};

// ========================================
// Scaled, Preconditioned, Flexible, Generalized Minimum Residual
// ========================================

class spfgmr_linsol : public linsol<spfgmr_linsol>
{
  public:
    spfgmr_linsol(N_Vector y, precondition pretype, int maxl)
    {
        m_ls = SUNSPFGMR(y, (int)pretype, maxl);
        IEXP_NOT_NULLPTR(m_ls);
    }

    spfgmr_linsol &precondition(precondition pretype)
    {
        SUNSPFGMRSetPrecType(m_ls, (int)pretype);
        return *this;
    }

    spfgmr_linsol &gram_schmidt(gram_schmidt gstype)
    {
        SUNSPFGMRSetGSType(m_ls, (int)gstype);
        return *this;
    }

    spfgmr_linsol &max_restart(int maxrs)
    {
        SUNSPFGMRSetMaxRestarts(m_ls, maxrs);
        return *this;
    }
};

// ========================================
// Scaled, Preconditioned,  Bi-Conjugate Gradient, Stabilized
// ========================================

class spbcgs_linsol : public linsol<spbcgs_linsol>
{
  public:
    spbcgs_linsol(N_Vector y, precondition pretype, int maxl)
    {
        m_ls = SUNSPBCGS(y, (int)pretype, maxl);
        IEXP_NOT_NULLPTR(m_ls);
    }

    spbcgs_linsol &precondition(precondition pretype)
    {
        SUNSPBCGSSetPrecType(m_ls, (int)pretype);
        return *this;
    }

    spbcgs_linsol &maxl(int maxl)
    {
        SUNSPBCGSSetMaxl(m_ls, maxl);
        return *this;
    }
};

// ========================================
//  Scaled, Preconditioned, Transpose-Free Quasi-Minimum Residual
// ========================================

class sptfqmr_linsol : public linsol<sptfqmr_linsol>
{
  public:
    sptfqmr_linsol(N_Vector y, precondition pretype, int maxl)
    {
        m_ls = SUNSPTFQMR(y, (int)pretype, maxl);
        IEXP_NOT_NULLPTR(m_ls);
    }

    sptfqmr_linsol &precondition(precondition pretype)
    {
        SUNSPTFQMRSetPrecType(m_ls, (int)pretype);
        return *this;
    }

    sptfqmr_linsol &maxl(int maxl)
    {
        SUNSPTFQMRSetMaxl(m_ls, maxl);
        return *this;
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

#endif /* __IEXP_DAE_LINSOL__ */
