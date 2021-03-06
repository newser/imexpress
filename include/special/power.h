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

#ifndef __IEXP_SF_POWER__
#define __IEXP_SF_POWER__

////////////////////////////////////////////////////////////
// import header files
////////////////////////////////////////////////////////////

#include <common/common.h>

#include <gsl/gsl_errno.h>
#include <gsl/gsl_sf_pow_int.h>

IEXP_NS_BEGIN

namespace sf {

////////////////////////////////////////////////////////////
// macro definition
////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////
// type definition
////////////////////////////////////////////////////////////

template <typename T>
class pow_functor : public functor_foreach<pow_functor<T>, T, double>
{
  public:
    pow_functor(int n, const T &x)
        : functor_foreach<pow_functor<T>, T, double>(x)
        , m_n(n)
    {
    }

    double foreach_impl(double x) const
    {
        return gsl_sf_pow_int(x, m_n);
    }

  private:
    int m_n;
};

template <typename T>
inline CwiseNullaryOp<pow_functor<T>, typename pow_functor<T>::ResultType> pow(
    int n, const DenseBase<T> &x)
{
    using ResultType = typename pow_functor<T>::ResultType;
    return ResultType::NullaryExpr(x.rows(),
                                   x.cols(),
                                   pow_functor<T>(n, x.derived()));
}

template <typename T, typename U>
class pow_e_functor
    : public functor_foreach_e<pow_e_functor<T, U>, T, U, double>
{
  public:
    pow_e_functor(int n, const T &x, U &e)
        : functor_foreach_e<pow_e_functor<T, U>, T, U, double>(x, e)
        , m_n(n)
    {
    }

    double foreach_e_impl(double x, double &e) const
    {
        gsl_sf_result r;
        gsl_sf_pow_int_e(x, m_n, &r);
        e = r.err;
        return r.val;
    }

  private:
    int m_n;
};

template <typename T, typename U>
inline CwiseNullaryOp<pow_e_functor<T, U>,
                      typename pow_e_functor<T, U>::ResultType>
pow(int n, const DenseBase<T> &x, DenseBase<U> &e)
{
    using ResultType = typename pow_e_functor<T, U>::ResultType;
    return ResultType::NullaryExpr(x.rows(),
                                   x.cols(),
                                   pow_e_functor<T, U>(n,
                                                       x.derived(),
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

#endif /* __IEXP_SF_POWER__ */
