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

enum sbessel_j
{
    j0,
    j1,
    j2,
    jn,
};

enum sbessel_y
{
    y0,
    y1,
    y2,
    yn,
};

// ========================================
// sbessel first kind
// ========================================

template <sbessel_j order, typename T>
inline T sbessel_j_impl(int n, T x)
{
    UNSUPPORTED_TYPE(T);
}

template <>
inline double sbessel_j_impl<j0, double>(int n, double x)
{
    return gsl_sf_bessel_j0(x);
}

template <>
inline double sbessel_j_impl<j1, double>(int n, double x)
{
    return gsl_sf_bessel_j1(x);
}

template <>
inline double sbessel_j_impl<j2, double>(int n, double x)
{
    return gsl_sf_bessel_j2(x);
}

template <>
inline double sbessel_j_impl<jn, double>(int n, double x)
{
    return gsl_sf_bessel_jl(n, x);
}

template <sbessel_j order, typename T>
class sbessel_j_functor
{
  public:
    using Scalar = typename T::Scalar;
    using ResultType = typename dense_derive<T>::type;

    sbessel_j_functor(int n, const T &x)
        : m_n(n)
        , m_x(x)
    {
    }

    Scalar operator()(Index i, Index j) const
    {
        return sbessel_j_impl<order, Scalar>(m_n, m_x(i, j));
    }

  private:
    int m_n;
    const T &m_x;
};

template <typename T>
inline CwiseNullaryOp<sbessel_j_functor<jn, T>,
                      typename sbessel_j_functor<jn, T>::ResultType>
sbessel_j(int n, const DenseBase<T> &x)
{
    using ResultType = typename sbessel_j_functor<jn, T>::ResultType;
    return ResultType::NullaryExpr(x.rows(),
                                   x.cols(),
                                   sbessel_j_functor<jn, T>(n, x.derived()));
}

template <enum sbessel_j order, typename T>
inline CwiseNullaryOp<sbessel_j_functor<order, T>,
                      typename sbessel_j_functor<order, T>::ResultType>
sbessel_j(const DenseBase<T> &x)
{
    static_assert(order < jn, "order can only be j0/j1/j2");

    using ResultType = typename sbessel_j_functor<order, T>::ResultType;
    return ResultType::NullaryExpr(x.rows(),
                                   x.cols(),
                                   sbessel_j_functor<order, T>(order,
                                                               x.derived()));
}

template <enum sbessel_j order, typename T>
inline T sbessel_j_e_impl(int n, T x, T &e)
{
    UNSUPPORTED_TYPE(T);
}

template <>
inline double sbessel_j_e_impl<j0, double>(int n, double x, double &e)
{
    gsl_sf_result r;
    gsl_sf_bessel_j0_e(x, &r);
    e = r.err;
    return r.val;
}

template <>
inline double sbessel_j_e_impl<j1, double>(int n, double x, double &e)
{
    gsl_sf_result r;
    gsl_sf_bessel_j1_e(x, &r);
    e = r.err;
    return r.val;
}

template <>
inline double sbessel_j_e_impl<j2, double>(int n, double x, double &e)
{
    gsl_sf_result r;
    gsl_sf_bessel_j2_e(x, &r);
    e = r.err;
    return r.val;
}

template <>
inline double sbessel_j_e_impl<jn, double>(int n, double x, double &e)
{
    gsl_sf_result r;
    gsl_sf_bessel_jl_e(n, x, &r);
    e = r.err;
    return r.val;
}

template <enum sbessel_j order, typename T, typename U>
class sbessel_j_e_functor
{
  public:
    using Scalar = typename T::Scalar;
    using ResultType = typename dense_derive<T>::type;

    sbessel_j_e_functor(int n, const T &x, U &e)
        : m_n(n)
        , m_x(x)
        , m_e(e)
    {
    }

    Scalar operator()(Index i, Index j) const
    {
        return sbessel_j_e_impl<order, Scalar>(m_n, m_x(i, j), m_e(i, j));
    }

  private:
    int m_n;
    const T &m_x;
    U &m_e;
};

template <typename T, typename U>
inline CwiseNullaryOp<sbessel_j_e_functor<jn, T, U>,
                      typename sbessel_j_e_functor<jn, T, U>::ResultType>
sbessel_j(int n, const DenseBase<T> &x, DenseBase<U> &e)
{
    using ResultType = typename sbessel_j_e_functor<jn, T, U>::ResultType;
    return ResultType::NullaryExpr(x.rows(),
                                   x.cols(),
                                   sbessel_j_e_functor<jn, T, U>(n,
                                                                 x.derived(),
                                                                 e.derived()));
}

template <enum sbessel_j order, typename T, typename U>
inline CwiseNullaryOp<sbessel_j_e_functor<order, T, U>,
                      typename sbessel_j_e_functor<order, T, U>::ResultType>
sbessel_j(const DenseBase<T> &x, DenseBase<U> &e)
{
    static_assert(order < jn, "order can only be j0/j1/j2");

    using ResultType = typename sbessel_j_e_functor<order, T, U>::ResultType;
    return ResultType::
        NullaryExpr(x.rows(),
                    x.cols(),
                    sbessel_j_e_functor<order, T, U>(order,
                                                     x.derived(),
                                                     e.derived()));
}

// ========================================
// bessel second kind
// ========================================

template <enum sbessel_y y, typename T>
inline T sbessel_y_impl(int n, T x)
{
    UNSUPPORTED_TYPE(T);
}

template <>
inline double sbessel_y_impl<y0, double>(int n, double x)
{
    return gsl_sf_bessel_y0(x);
}

template <>
inline double sbessel_y_impl<y1, double>(int n, double x)
{
    return gsl_sf_bessel_y1(x);
}

template <>
inline double sbessel_y_impl<y2, double>(int n, double x)
{
    return gsl_sf_bessel_y2(x);
}

template <>
inline double sbessel_y_impl<yn, double>(int n, double x)
{
    return gsl_sf_bessel_yl(n, x);
}

template <enum sbessel_y order, typename T>
class sbessel_y_functor
{
  public:
    using Scalar = typename T::Scalar;
    using ResultType = typename dense_derive<T>::type;

    sbessel_y_functor(int n, const T &x)
        : m_n(n)
        , m_x(x)
    {
    }

    Scalar operator()(Index i, Index j) const
    {
        return sbessel_y_impl<order, Scalar>(m_n, m_x(i, j));
    }

  private:
    int m_n;
    const T &m_x;
};

template <typename T>
inline CwiseNullaryOp<sbessel_y_functor<yn, T>,
                      typename sbessel_y_functor<yn, T>::ResultType>
sbessel_y(int n, const DenseBase<T> &x)
{
    using ResultType = typename sbessel_y_functor<yn, T>::ResultType;
    return ResultType::NullaryExpr(x.rows(),
                                   x.cols(),
                                   sbessel_y_functor<yn, T>(n, x.derived()));
}

template <enum sbessel_y order, typename T>
inline CwiseNullaryOp<sbessel_y_functor<order, T>,
                      typename sbessel_y_functor<order, T>::ResultType>
sbessel_y(const DenseBase<T> &x)
{
    static_assert(order < yn, "order can only be y0/y1/y2");

    using ResultType = typename sbessel_y_functor<order, T>::ResultType;
    return ResultType::NullaryExpr(x.rows(),
                                   x.cols(),
                                   sbessel_y_functor<order, T>(order,
                                                               x.derived()));
}

template <typename T>
inline T sbessel_y_e_impl(int n, T x, T &e)
{
    UNSUPPORTED_TYPE(T);
}

template <enum sbessel_y order, typename T>
inline T sbessel_y_e_impl(int n, T x, T &e)
{
    UNSUPPORTED_TYPE(T);
}

template <>
inline double sbessel_y_e_impl<y0, double>(int n, double x, double &e)
{
    gsl_sf_result r;
    gsl_sf_bessel_y0_e(x, &r);
    e = r.err;
    return r.val;
}

template <>
inline double sbessel_y_e_impl<y1, double>(int n, double x, double &e)
{
    gsl_sf_result r;
    gsl_sf_bessel_y1_e(x, &r);
    e = r.err;
    return r.val;
}

template <>
inline double sbessel_y_e_impl<y2, double>(int n, double x, double &e)
{
    gsl_sf_result r;
    gsl_sf_bessel_y2_e(x, &r);
    e = r.err;
    return r.val;
}

template <>
inline double sbessel_y_e_impl<yn, double>(int n, double x, double &e)
{
    gsl_sf_result r;
    gsl_sf_bessel_yl_e(n, x, &r);
    e = r.err;
    return r.val;
}

template <enum sbessel_y order, typename T, typename U>
class sbessel_y_e_functor
{
  public:
    using Scalar = typename T::Scalar;
    using ResultType = typename dense_derive<T>::type;

    sbessel_y_e_functor(int n, const T &x, U &e)
        : m_n(n)
        , m_x(x)
        , m_e(e)
    {
    }

    Scalar operator()(Index i, Index j) const
    {
        return sbessel_y_e_impl<order, Scalar>(m_n, m_x(i, j), m_e(i, j));
    }

  private:
    int m_n;
    const T &m_x;
    U &m_e;
};

template <typename T, typename U>
inline CwiseNullaryOp<sbessel_y_e_functor<yn, T, U>,
                      typename sbessel_y_e_functor<yn, T, U>::ResultType>
sbessel_y(int n, const DenseBase<T> &x, DenseBase<U> &e)
{
    using ResultType = typename sbessel_y_e_functor<yn, T, U>::ResultType;
    return ResultType::NullaryExpr(x.rows(),
                                   x.cols(),
                                   sbessel_y_e_functor<yn, T, U>(n,
                                                                 x.derived(),
                                                                 e.derived()));
}

template <enum sbessel_y order, typename T, typename U>
inline CwiseNullaryOp<sbessel_y_e_functor<order, T, U>,
                      typename sbessel_y_e_functor<order, T, U>::ResultType>
sbessel_y(const DenseBase<T> &x, DenseBase<U> &e)
{
    static_assert(order < yn, "order can only be y0/y1/y2");

    using ResultType = typename sbessel_y_e_functor<order, T, U>::ResultType;
    return ResultType::
        NullaryExpr(x.rows(),
                    x.cols(),
                    sbessel_y_e_functor<order, T, U>(order,
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
