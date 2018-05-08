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
    j_n,
};

enum sbessel_y
{
    y0,
    y1,
    y2,
    y_n,
};

// ========================================
// sbessel first kind
// ========================================

template <sbessel_j order>
inline double sbessel_j_impl(int n, double x)
{
    static_assert(order == j_n, "order must be j_n");
    return gsl_sf_bessel_jl(n, x);
}

template <>
inline double sbessel_j_impl<j0>(int n, double x)
{
    return gsl_sf_bessel_j0(x);
}

template <>
inline double sbessel_j_impl<j1>(int n, double x)
{
    return gsl_sf_bessel_j1(x);
}

template <>
inline double sbessel_j_impl<j2>(int n, double x)
{
    return gsl_sf_bessel_j2(x);
}

template <sbessel_j order, typename T>
class sbessel_j_functor
    : public functor_foreach<sbessel_j_functor<order, T>, T, double>
{
  public:
    sbessel_j_functor(int n, const T &x)
        : functor_foreach<sbessel_j_functor<order, T>, T, double>(x)
        , m_n(n)
    {
    }

    double foreach_impl(double x) const
    {
        return sbessel_j_impl<order>(m_n, x);
    }

  private:
    int m_n;
};

template <typename T>
inline CwiseNullaryOp<sbessel_j_functor<j_n, T>,
                      typename sbessel_j_functor<j_n, T>::ResultType>
sbessel_j(int n, const DenseBase<T> &x)
{
    using ResultType = typename sbessel_j_functor<j_n, T>::ResultType;
    return ResultType::NullaryExpr(x.rows(),
                                   x.cols(),
                                   sbessel_j_functor<j_n, T>(n, x.derived()));
}

template <enum sbessel_j order, typename T>
inline CwiseNullaryOp<sbessel_j_functor<order, T>,
                      typename sbessel_j_functor<order, T>::ResultType>
sbessel_j(const DenseBase<T> &x)
{
    static_assert(order < j_n, "order can only be j0/j1/j2");

    using ResultType = typename sbessel_j_functor<order, T>::ResultType;
    return ResultType::NullaryExpr(x.rows(),
                                   x.cols(),
                                   sbessel_j_functor<order, T>(order,
                                                               x.derived()));
}

template <enum sbessel_j order>
inline double sbessel_j_e_impl(int n, double x, double &e)
{
    static_assert(order == j_n, "order must be j_n");

    gsl_sf_result r;
    gsl_sf_bessel_jl_e(n, x, &r);
    e = r.err;
    return r.val;
}

template <>
inline double sbessel_j_e_impl<j0>(int n, double x, double &e)
{
    gsl_sf_result r;
    gsl_sf_bessel_j0_e(x, &r);
    e = r.err;
    return r.val;
}

template <>
inline double sbessel_j_e_impl<j1>(int n, double x, double &e)
{
    gsl_sf_result r;
    gsl_sf_bessel_j1_e(x, &r);
    e = r.err;
    return r.val;
}

template <>
inline double sbessel_j_e_impl<j2>(int n, double x, double &e)
{
    gsl_sf_result r;
    gsl_sf_bessel_j2_e(x, &r);
    e = r.err;
    return r.val;
}

template <enum sbessel_j order, typename T, typename U>
class sbessel_j_e_functor
    : public functor_foreach_e<sbessel_j_e_functor<order, T, U>, T, U, double>
{
  public:
    sbessel_j_e_functor(int n, const T &x, U &e)
        : functor_foreach_e<sbessel_j_e_functor<order, T, U>, T, U, double>(x,
                                                                            e)
        , m_n(n)
    {
    }

    double foreach_e_impl(double x, double &e) const
    {
        return sbessel_j_e_impl<order>(m_n, x, e);
    }

  private:
    int m_n;
};

template <typename T, typename U>
inline CwiseNullaryOp<sbessel_j_e_functor<j_n, T, U>,
                      typename sbessel_j_e_functor<j_n, T, U>::ResultType>
sbessel_j(int n, const DenseBase<T> &x, DenseBase<U> &e)
{
    using ResultType = typename sbessel_j_e_functor<j_n, T, U>::ResultType;
    return ResultType::NullaryExpr(x.rows(),
                                   x.cols(),
                                   sbessel_j_e_functor<j_n, T, U>(n,
                                                                  x.derived(),
                                                                  e.derived()));
}

template <enum sbessel_j order, typename T, typename U>
inline CwiseNullaryOp<sbessel_j_e_functor<order, T, U>,
                      typename sbessel_j_e_functor<order, T, U>::ResultType>
sbessel_j(const DenseBase<T> &x, DenseBase<U> &e)
{
    static_assert(order < j_n, "order can only be j0/j1/j2");

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

template <enum sbessel_y order>
inline double sbessel_y_impl(int n, double x)
{
    static_assert(order == y_n, "order must be j_n");

    return gsl_sf_bessel_yl(n, x);
}

template <>
inline double sbessel_y_impl<y0>(int n, double x)
{
    return gsl_sf_bessel_y0(x);
}

template <>
inline double sbessel_y_impl<y1>(int n, double x)
{
    return gsl_sf_bessel_y1(x);
}

template <>
inline double sbessel_y_impl<y2>(int n, double x)
{
    return gsl_sf_bessel_y2(x);
}

template <enum sbessel_y order, typename T>
class sbessel_y_functor
    : public functor_foreach<sbessel_y_functor<order, T>, T, double>
{
  public:
    sbessel_y_functor(int n, const T &x)
        : functor_foreach<sbessel_y_functor<order, T>, T, double>(x)
        , m_n(n)
    {
    }

    double foreach_impl(double x) const
    {
        return sbessel_y_impl<order>(m_n, x);
    }

  private:
    int m_n;
};

template <typename T>
inline CwiseNullaryOp<sbessel_y_functor<y_n, T>,
                      typename sbessel_y_functor<y_n, T>::ResultType>
sbessel_y(int n, const DenseBase<T> &x)
{
    using ResultType = typename sbessel_y_functor<y_n, T>::ResultType;
    return ResultType::NullaryExpr(x.rows(),
                                   x.cols(),
                                   sbessel_y_functor<y_n, T>(n, x.derived()));
}

template <enum sbessel_y order, typename T>
inline CwiseNullaryOp<sbessel_y_functor<order, T>,
                      typename sbessel_y_functor<order, T>::ResultType>
sbessel_y(const DenseBase<T> &x)
{
    static_assert(order < y_n, "order can only be y0/y1/y2");

    using ResultType = typename sbessel_y_functor<order, T>::ResultType;
    return ResultType::NullaryExpr(x.rows(),
                                   x.cols(),
                                   sbessel_y_functor<order, T>(order,
                                                               x.derived()));
}

template <enum sbessel_y order>
inline double sbessel_y_e_impl(int n, double x, double &e)
{
    static_assert(order == y_n, "order must be y_n");

    gsl_sf_result r;
    gsl_sf_bessel_yl_e(n, x, &r);
    e = r.err;
    return r.val;
}

template <>
inline double sbessel_y_e_impl<y0>(int n, double x, double &e)
{
    gsl_sf_result r;
    gsl_sf_bessel_y0_e(x, &r);
    e = r.err;
    return r.val;
}

template <>
inline double sbessel_y_e_impl<y1>(int n, double x, double &e)
{
    gsl_sf_result r;
    gsl_sf_bessel_y1_e(x, &r);
    e = r.err;
    return r.val;
}

template <>
inline double sbessel_y_e_impl<y2>(int n, double x, double &e)
{
    gsl_sf_result r;
    gsl_sf_bessel_y2_e(x, &r);
    e = r.err;
    return r.val;
}

template <enum sbessel_y order, typename T, typename U>
class sbessel_y_e_functor
    : public functor_foreach_e<sbessel_y_e_functor<order, T, U>, T, U, double>
{
  public:
    sbessel_y_e_functor(int n, const T &x, U &e)
        : functor_foreach_e<sbessel_y_e_functor<order, T, U>, T, U, double>(x,
                                                                            e)
        , m_n(n)
    {
    }

    double foreach_e_impl(double x, double &e) const
    {
        return sbessel_y_e_impl<order>(m_n, x, e);
    }

  private:
    int m_n;
};

template <typename T, typename U>
inline CwiseNullaryOp<sbessel_y_e_functor<y_n, T, U>,
                      typename sbessel_y_e_functor<y_n, T, U>::ResultType>
sbessel_y(int n, const DenseBase<T> &x, DenseBase<U> &e)
{
    using ResultType = typename sbessel_y_e_functor<y_n, T, U>::ResultType;
    return ResultType::NullaryExpr(x.rows(),
                                   x.cols(),
                                   sbessel_y_e_functor<y_n, T, U>(n,
                                                                  x.derived(),
                                                                  e.derived()));
}

template <enum sbessel_y order, typename T, typename U>
inline CwiseNullaryOp<sbessel_y_e_functor<order, T, U>,
                      typename sbessel_y_e_functor<order, T, U>::ResultType>
sbessel_y(const DenseBase<T> &x, DenseBase<U> &e)
{
    static_assert(order < y_n, "order can only be y0/y1/y2");

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
