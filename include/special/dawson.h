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
inline T dawson_impl(const T x)
{
    UNSUPPORTED_TYPE(T);
}

template <>
inline double dawson_impl(const double x)
{
    return gsl_sf_dawson(x);
}

template <typename T>
class dawson_functor
{
  public:
    typedef Array<typename T::Scalar,
                  T::RowsAtCompileTime,
                  T::ColsAtCompileTime,
                  T::Flags & RowMajorBit ? RowMajor : ColMajor,
                  T::MaxRowsAtCompileTime,
                  T::MaxColsAtCompileTime>
        ArrayType;

    dawson_functor(const T &x)
        : m_x(x)
    {
    }

    const typename T::Scalar operator()(Index i, Index j) const
    {
        return dawson_impl(m_x(i, j));
    }

  private:
    const T &m_x;
};

template <typename T>
inline CwiseNullaryOp<dawson_functor<T>, typename dawson_functor<T>::ArrayType>
dawson(const ArrayBase<T> &x)
{
    typedef typename dawson_functor<T>::ArrayType ArrayType;
    return ArrayType::NullaryExpr(x.rows(),
                                  x.cols(),
                                  dawson_functor<T>(x.derived()));
}

template <typename T>
inline T dawson_e_impl(const T x, T &e)
{
    UNSUPPORTED_TYPE(T);
}

template <>
inline double dawson_e_impl(const double x, double &e)
{
    gsl_sf_result r;
    if (gsl_sf_dawson_e(x, &r) == GSL_SUCCESS) {
        e = r.err;
        return r.val;
    }
    RETURN_NAN_OR_THROW(std::runtime_error("dawson"));
}

template <typename T, typename U>
class dawson_e_functor
{
  public:
    typedef Array<typename T::Scalar,
                  T::RowsAtCompileTime,
                  T::ColsAtCompileTime,
                  T::Flags & RowMajorBit ? RowMajor : ColMajor,
                  T::MaxRowsAtCompileTime,
                  T::MaxColsAtCompileTime>
        ArrayType;

    dawson_e_functor(const T &x, U &e)
        : m_x(x)
        , m_e(e)
    {
    }

    const typename T::Scalar operator()(Index i, Index j) const
    {
        return dawson_e_impl(m_x(i, j), m_e(i, j));
    }

  private:
    const T &m_x;
    U &m_e;
};

template <typename T, typename U>
inline CwiseNullaryOp<dawson_e_functor<T, U>,
                      typename dawson_e_functor<T, U>::ArrayType>
dawson(const ArrayBase<T> &x, ArrayBase<U> &e)
{
    typedef typename dawson_e_functor<T, U>::ArrayType ArrayType;
    return ArrayType::NullaryExpr(x.rows(),
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