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

#ifndef __IEXP_ELLIPTIC_JACOBI__
#define __IEXP_ELLIPTIC_JACOBI__

////////////////////////////////////////////////////////////
// import header files
////////////////////////////////////////////////////////////

#include <common/common.h>

#include <gsl/gsl_errno.h>
#include <gsl/gsl_sf_elljac.h>

IEXP_NS_BEGIN

namespace sf {

////////////////////////////////////////////////////////////
// macro definition
////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////
// type definition
////////////////////////////////////////////////////////////

template <typename T>
class elljac_functor
    : public functor_m2vnum_2d<elljac_functor<T>,
                               T,
                               std::tuple<double, double, double>>
{
  public:
    elljac_functor(const T &u_m)
        : functor_m2vnum_2d<elljac_functor<T>,
                            T,
                            std::tuple<double, double, double>>(u_m)
    {
    }

    std::tuple<double, double, double> m2vnum_impl(double u, double m) const
    {
        double sn, cn, dn;
        gsl_sf_elljac_e(u, m, &sn, &cn, &dn);
        return std::make_tuple(sn, cn, dn);
    }
};

template <typename T>
inline CwiseNullaryOp<elljac_functor<T>, typename elljac_functor<T>::ResultType>
elljac(const DenseBase<T> &u_m)
{
    using ResultType = typename elljac_functor<T>::ResultType;
    return ResultType::NullaryExpr(M2VNUM_ROW(T, u_m),
                                   M2VNUM_COL(T, u_m),
                                   elljac_functor<T>(u_m.derived()));
}

////////////////////////////////////////////////////////////
// global variants
////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////
// interface declaration
////////////////////////////////////////////////////////////
}

IEXP_NS_END

#endif /* __IEXP_ELLIPTIC_JACOBI__ */
