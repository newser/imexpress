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

#ifndef __IEXP_SPHERICAL_BESSEL__
#define __IEXP_SPHERICAL_BESSEL__

////////////////////////////////////////////////////////////
// import header files
////////////////////////////////////////////////////////////

#include <common/common.h>

#include <gsl/gsl_errno.h>
#include <gsl/gsl_sf_bessel.h>

IEXP_NS_BEGIN

namespace sf {

////////////////////////////////////////////////////////////
// macro definition
////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////
// type definition
////////////////////////////////////////////////////////////

// ========================================
// bessel first kind
// ========================================

template <typename T>
inline T sbessel_j_impl(const int n, const T x)
{
    UNSUPPORTED_TYPE(T);
}

template <>
inline double sbessel_j_impl(const int n, const double x)
{
    // gsl_sf_bessel_jl does not allow negative x...
    if (n == 0) {
        return gsl_sf_bessel_j0(x);
    } else if (n == 1) {
        return gsl_sf_bessel_j1(x);
    } else if (n == 2) {
        return gsl_sf_bessel_j2(x);
    } else {
        return gsl_sf_bessel_jl(n, x);
    }
}

template <typename T>
class sbessel_j_functor
{
  public:
    typedef Array<typename T::Scalar,
                  T::RowsAtCompileTime,
                  T::ColsAtCompileTime,
                  T::Flags & RowMajorBit ? RowMajor : ColMajor,
                  T::MaxRowsAtCompileTime,
                  T::MaxColsAtCompileTime>
        ArrayType;

    sbessel_j_functor(const int n, const T &x)
        : m_n(n)
        , m_x(x)
    {
    }

    const typename T::Scalar operator()(Index i, Index j) const
    {
        return sbessel_j_impl<typename T::Scalar>(m_n, m_x(i, j));
    }

  private:
    const int m_n;
    const T &m_x;
};

template <typename T>
inline CwiseNullaryOp<sbessel_j_functor<T>,
                      typename sbessel_j_functor<T>::ArrayType>
sbessel_j(const int n, const ArrayBase<T> &x)
{
    typedef typename sbessel_j_functor<T>::ArrayType ArrayType;
    return ArrayType::NullaryExpr(x.rows(),
                                  x.cols(),
                                  sbessel_j_functor<T>(n, x.derived()));
}

template <typename T>
inline T sbessel_j_e_impl(const int n, const T x, T &e)
{
    UNSUPPORTED_TYPE(T);
}

template <>
inline double sbessel_j_e_impl(const int n, const double x, double &e)
{
    gsl_sf_result r;
    int s;
    // gsl_sf_bessel_jl does not allow negative x...
    if (n == 0) {
        s = gsl_sf_bessel_j0_e(x, &r);
    } else if (n == 1) {
        s = gsl_sf_bessel_j1_e(x, &r);
    } else if (n == 2) {
        s = gsl_sf_bessel_j2_e(x, &r);
    } else {
        s = gsl_sf_bessel_jl_e(n, x, &r);
    }
    if (s == GSL_SUCCESS) {
        e = r.err;
        return r.val;
    }
    THROW_OR_RETURN_NAN(std::runtime_error("sbessel_j"));
}

template <typename T, typename U>
class sbessel_j_e_functor
{
  public:
    typedef Array<typename T::Scalar,
                  T::RowsAtCompileTime,
                  T::ColsAtCompileTime,
                  T::Flags & RowMajorBit ? RowMajor : ColMajor,
                  T::MaxRowsAtCompileTime,
                  T::MaxColsAtCompileTime>
        ArrayType;

    sbessel_j_e_functor(const int n, const T &x, U &e)
        : m_n(n)
        , m_x(x)
        , m_e(e)
    {
    }

    const typename T::Scalar operator()(Index i, Index j) const
    {
        return sbessel_j_e_impl<typename T::Scalar>(m_n, m_x(i, j), m_e(i, j));
    }

  private:
    const int m_n;
    const T &m_x;
    U &m_e;
};

template <typename T, typename U>
inline CwiseNullaryOp<sbessel_j_e_functor<T, U>,
                      typename sbessel_j_e_functor<T, U>::ArrayType>
sbessel_j(const int n, const ArrayBase<T> &x, ArrayBase<U> &e)
{
    typedef typename sbessel_j_e_functor<T, U>::ArrayType ArrayType;
    return ArrayType::NullaryExpr(x.rows(),
                                  x.cols(),
                                  sbessel_j_e_functor<T, U>(n,
                                                            x.derived(),
                                                            e.derived()));
}

// ========================================
// bessel second kind
// ========================================

template <typename T>
inline T sbessel_y_impl(const int n, const T x)
{
    UNSUPPORTED_TYPE(T);
}

template <>
inline double sbessel_y_impl(const int n, const double x)
{
    return gsl_sf_bessel_yl(n, x);
}

template <typename T>
class sbessel_y_functor
{
  public:
    typedef Array<typename T::Scalar,
                  T::RowsAtCompileTime,
                  T::ColsAtCompileTime,
                  T::Flags & RowMajorBit ? RowMajor : ColMajor,
                  T::MaxRowsAtCompileTime,
                  T::MaxColsAtCompileTime>
        ArrayType;

    sbessel_y_functor(const int n, const T &x)
        : m_n(n)
        , m_x(x)
    {
    }

    const typename T::Scalar operator()(Index i, Index j) const
    {
        return sbessel_y_impl<typename T::Scalar>(m_n, m_x(i, j));
    }

  private:
    const int m_n;
    const T &m_x;
};

template <typename T>
inline CwiseNullaryOp<sbessel_y_functor<T>,
                      typename sbessel_y_functor<T>::ArrayType>
sbessel_y(const int n, const ArrayBase<T> &x)
{
    typedef typename sbessel_y_functor<T>::ArrayType ArrayType;
    return ArrayType::NullaryExpr(x.rows(),
                                  x.cols(),
                                  sbessel_y_functor<T>(n, x.derived()));
}

template <typename T>
inline T sbessel_y_e_impl(const int n, const T x, T &e)
{
    UNSUPPORTED_TYPE(T);
}

template <>
inline double sbessel_y_e_impl(const int n, const double x, double &e)
{
    gsl_sf_result r;
    if (gsl_sf_bessel_yl_e(n, x, &r) == GSL_SUCCESS) {
        e = r.err;
        return r.val;
    }
    THROW_OR_RETURN_NAN(std::runtime_error("sbessel_y"));
}

template <typename T, typename U>
class sbessel_y_e_functor
{
  public:
    typedef Array<typename T::Scalar,
                  T::RowsAtCompileTime,
                  T::ColsAtCompileTime,
                  T::Flags & RowMajorBit ? RowMajor : ColMajor,
                  T::MaxRowsAtCompileTime,
                  T::MaxColsAtCompileTime>
        ArrayType;

    sbessel_y_e_functor(const int n, const T &x, U &e)
        : m_n(n)
        , m_x(x)
        , m_e(e)
    {
    }

    const typename T::Scalar operator()(Index i, Index j) const
    {
        return sbessel_y_e_impl<typename T::Scalar>(m_n, m_x(i, j), m_e(i, j));
    }

  private:
    const int m_n;
    const T &m_x;
    U &m_e;
};

template <typename T, typename U>
inline CwiseNullaryOp<sbessel_y_e_functor<T, U>,
                      typename sbessel_y_e_functor<T, U>::ArrayType>
sbessel_y(const int n, const ArrayBase<T> &x, ArrayBase<U> &e)
{
    typedef typename sbessel_y_e_functor<T, U>::ArrayType ArrayType;
    return ArrayType::NullaryExpr(x.rows(),
                                  x.cols(),
                                  sbessel_y_e_functor<T, U>(n,
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

#endif /* __IEXP_SPHERICAL_BESSEL__ */
