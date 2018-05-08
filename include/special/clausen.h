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

#ifndef __IEXP_CLAUSEN__
#define __IEXP_CLAUSEN__

////////////////////////////////////////////////////////////
// import header files
////////////////////////////////////////////////////////////

#include <common/common.h>

#include <gsl/gsl_errno.h>
#include <gsl/gsl_sf_clausen.h>

IEXP_NS_BEGIN

namespace sf {

////////////////////////////////////////////////////////////
// macro definition
////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////
// type definition
////////////////////////////////////////////////////////////

template <typename T>
class clausen_functor : public functor_foreach<clausen_functor<T>, T, double>
{
  public:
    clausen_functor(const T &x)
        : functor_foreach<clausen_functor<T>, T, double>(x)
    {
    }

    double foreach_impl(double x) const
    {
        return gsl_sf_clausen(x);
    }
};

template <typename T>
inline CwiseNullaryOp<clausen_functor<T>,
                      typename clausen_functor<T>::ResultType>
clausen(const DenseBase<T> &x)
{
    using ResultType = typename clausen_functor<T>::ResultType;
    return ResultType::NullaryExpr(x.rows(),
                                   x.cols(),
                                   clausen_functor<T>(x.derived()));
}

template <typename T, typename U>
class clausen_e_functor
    : public functor_foreach_e<clausen_e_functor<T, U>, T, U, double>
{
  public:
    clausen_e_functor(const T &x, U &e)
        : functor_foreach_e<clausen_e_functor<T, U>, T, U, double>(x, e)
    {
    }

    double foreach_e_impl(double x, double &e) const
    {
        gsl_sf_result r;
        gsl_sf_clausen_e(x, &r);
        e = r.err;
        return r.val;
    }
};

template <typename T, typename U>
inline CwiseNullaryOp<clausen_e_functor<T, U>,
                      typename clausen_e_functor<T, U>::ResultType>
clausen(const DenseBase<T> &x, DenseBase<U> &e)
{
    using ResultType = typename clausen_e_functor<T, U>::ResultType;
    return ResultType::NullaryExpr(x.rows(),
                                   x.cols(),
                                   clausen_e_functor<T, U>(x.derived(),
                                                           e.derived()));
}

////////////////////////////////////////////////////////////
// global variants
////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////
// interface declaration
////////////////////////////////////////////////////////////
}

IEXP_NS_END

#endif /* __IEXP_CLAUSEN__ */
