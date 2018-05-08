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

#ifndef __IEXP_DAWSON__
#define __IEXP_DAWSON__

////////////////////////////////////////////////////////////
// import header files
////////////////////////////////////////////////////////////

#include <common/common.h>

#include <gsl/gsl_errno.h>
#include <gsl/gsl_sf_dawson.h>

IEXP_NS_BEGIN

namespace sf {

////////////////////////////////////////////////////////////
// macro definition
////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////
// type definition
////////////////////////////////////////////////////////////

template <typename T>
class dawson_functor : public functor_foreach<dawson_functor<T>, T, double>
{
  public:
    dawson_functor(const T &x)
        : functor_foreach<dawson_functor<T>, T, double>(x)
    {
    }

    double foreach_impl(double x) const
    {
        return gsl_sf_dawson(x);
    }
};

template <typename T>
inline CwiseNullaryOp<dawson_functor<T>, typename dawson_functor<T>::ResultType>
dawson(const DenseBase<T> &x)
{
    using ResultType = typename dawson_functor<T>::ResultType;
    return ResultType::NullaryExpr(x.rows(),
                                   x.cols(),
                                   dawson_functor<T>(x.derived()));
}

template <typename T, typename U>
class dawson_e_functor
    : public functor_foreach_e<dawson_e_functor<T, U>, T, U, double>
{
  public:
    dawson_e_functor(const T &x, U &e)
        : functor_foreach_e<dawson_e_functor<T, U>, T, U, double>(x, e)
    {
    }

    double foreach_e_impl(double x, double &e) const
    {
        gsl_sf_result r;
        gsl_sf_dawson_e(x, &r);
        e = r.err;
        return r.val;
    }
};

template <typename T, typename U>
inline CwiseNullaryOp<dawson_e_functor<T, U>,
                      typename dawson_e_functor<T, U>::ResultType>
dawson(const DenseBase<T> &x, DenseBase<U> &e)
{
    using ResultType = typename dawson_e_functor<T, U>::ResultType;
    return ResultType::NullaryExpr(x.rows(),
                                   x.cols(),
                                   dawson_e_functor<T, U>(x.derived(),
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

#endif /* __IEXP_DAWSON__ */
