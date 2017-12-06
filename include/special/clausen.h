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
inline T clausen_impl(const T x)
{
    UNSUPPORTED_TYPE(T);
}

template <>
inline double clausen_impl(const double x)
{
    return gsl_sf_clausen(x);
}

template <typename T>
class clausen_functor
{
  public:
    using ArrayType = Array<typename T::Scalar,
                            T::RowsAtCompileTime,
                            T::ColsAtCompileTime,
                            T::Flags & RowMajorBit ? RowMajor : ColMajor,
                            T::MaxRowsAtCompileTime,
                            T::MaxColsAtCompileTime>;

    clausen_functor(const T &x)
        : m_x(x)
    {
    }

    const typename T::Scalar operator()(Index i, Index j) const
    {
        return clausen_impl(m_x(i, j));
    }

  private:
    const T &m_x;
};

template <typename T>
inline CwiseNullaryOp<clausen_functor<T>,
                      typename clausen_functor<T>::ArrayType>
clausen(const ArrayBase<T> &x)
{
    using ArrayType = typename clausen_functor<T>::ArrayType;
    return ArrayType::NullaryExpr(x.rows(),
                                  x.cols(),
                                  clausen_functor<T>(x.derived()));
}

template <typename T>
inline T clausen_e_impl(const T x, T &e)
{
    UNSUPPORTED_TYPE(T);
}

template <>
inline double clausen_e_impl(const double x, double &e)
{
    gsl_sf_result r;
    if (gsl_sf_clausen_e(x, &r) == GSL_SUCCESS) {
        e = r.err;
        return r.val;
    }
    RETURN_NAN_OR_THROW(std::runtime_error("clausen"));
}

template <typename T, typename U>
class clausen_e_functor
{
  public:
    using ArrayType = Array<typename T::Scalar,
                            T::RowsAtCompileTime,
                            T::ColsAtCompileTime,
                            T::Flags & RowMajorBit ? RowMajor : ColMajor,
                            T::MaxRowsAtCompileTime,
                            T::MaxColsAtCompileTime>;

    clausen_e_functor(const T &x, U &e)
        : m_x(x)
        , m_e(e)
    {
    }

    const typename T::Scalar operator()(Index i, Index j) const
    {
        return clausen_e_impl(m_x(i, j), m_e(i, j));
    }

  private:
    const T &m_x;
    U &m_e;
};

template <typename T, typename U>
inline CwiseNullaryOp<clausen_e_functor<T, U>,
                      typename clausen_e_functor<T, U>::ArrayType>
clausen(const ArrayBase<T> &x, ArrayBase<U> &e)
{
    using ArrayType = typename clausen_e_functor<T, U>::ArrayType;
    return ArrayType::NullaryExpr(x.rows(),
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
