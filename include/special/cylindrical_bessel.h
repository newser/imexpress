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

#ifndef __IEXP_CYLINDRICAL_BESSEL__
#define __IEXP_CYLINDRICAL_BESSEL__

////////////////////////////////////////////////////////////
// import header files
////////////////////////////////////////////////////////////

#include <common/common.h>

#include <common/template_foreach.h>

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

enum cbessel_j
{
    J0,
    J1,
    J_n,
};

enum cbessel_y
{
    Y0,
    Y1,
    Y2,
    Y_n,
};

// ========================================
// cbessel first kind
// ========================================

template <cbessel_j order, typename N>
inline double cbessel_j_impl(N n, double x)
{
    throw std::invalid_argument("invalid template parameters");
}

template <>
inline double cbessel_j_impl<J0, int>(int n, double x)
{
    return gsl_sf_bessel_J0(x);
}

template <>
inline double cbessel_j_impl<J1, int>(int n, double x)
{
    return gsl_sf_bessel_J1(x);
}

template <>
inline double cbessel_j_impl<J_n, int>(int n, double x)
{
    return gsl_sf_bessel_Jn(n, x);
}

template <>
inline double cbessel_j_impl<J_n, double>(double n, double x)
{
    return gsl_sf_bessel_Jnu(n, x);
}

template <cbessel_j order, typename T, typename N>
class cbessel_j_functor
    : public functor_foreach<cbessel_j_functor<order, T, N>, T, double>
{
  public:
    cbessel_j_functor(N n, const T &x)
        : functor_foreach<cbessel_j_functor<order, T, N>, T, double>(x)
        , m_n(n)
    {
    }

    double foreach_impl(double x) const
    {
        return cbessel_j_impl<order, N>(m_n, x);
    }

  private:
    N m_n;
};

template <typename T>
inline CwiseNullaryOp<cbessel_j_functor<J_n, T, int>,
                      typename cbessel_j_functor<J_n, T, int>::ResultType>
cbessel_j(int n, const DenseBase<T> &x)
{
    using ResultType = typename cbessel_j_functor<J_n, T, int>::ResultType;
    return ResultType::NullaryExpr(x.rows(),
                                   x.cols(),
                                   cbessel_j_functor<J_n, T, int>(n,
                                                                  x.derived()));
}

template <typename T>
inline CwiseNullaryOp<cbessel_j_functor<J_n, T, double>,
                      typename cbessel_j_functor<J_n, T, double>::ResultType>
cbessel_j(double n, const DenseBase<T> &x)
{
    using ResultType = typename cbessel_j_functor<J_n, T, double>::ResultType;
    return ResultType::
        NullaryExpr(x.rows(),
                    x.cols(),
                    cbessel_j_functor<J_n, T, double>(n, x.derived()));
}

template <enum cbessel_j order, typename T>
inline CwiseNullaryOp<cbessel_j_functor<order, T, int>,
                      typename cbessel_j_functor<order, T, int>::ResultType>
cbessel_j(const DenseBase<T> &x)
{
    static_assert(order < J_n, "order can only be J0/J1/J2");

    using ResultType = typename cbessel_j_functor<order, T, int>::ResultType;
    return ResultType::
        NullaryExpr(x.rows(),
                    x.cols(),
                    cbessel_j_functor<order, T, int>(order, x.derived()));
}

template <enum cbessel_j order, typename N>
inline double cbessel_j_e_impl(N n, double x, double &e)
{
    throw std::invalid_argument("invalid template parameters");
}

template <>
inline double cbessel_j_e_impl<J0, int>(int n, double x, double &e)
{
    gsl_sf_result r;
    gsl_sf_bessel_J0_e(x, &r);
    e = r.err;
    return r.val;
}

template <>
inline double cbessel_j_e_impl<J1, int>(int n, double x, double &e)
{
    gsl_sf_result r;
    gsl_sf_bessel_J1_e(x, &r);
    e = r.err;
    return r.val;
}

template <>
inline double cbessel_j_e_impl<J_n, int>(int n, double x, double &e)
{
    gsl_sf_result r;
    gsl_sf_bessel_Jn_e(n, x, &r);
    e = r.err;
    return r.val;
}

template <>
inline double cbessel_j_e_impl<J_n, double>(double n, double x, double &e)
{
    gsl_sf_result r;
    gsl_sf_bessel_Jnu_e(n, x, &r);
    e = r.err;
    return r.val;
}

template <enum cbessel_j order, typename T, typename U, typename N>
class cbessel_j_e_functor
    : public functor_foreach_e<cbessel_j_e_functor<order, T, U, N>,
                               T,
                               U,
                               double>
{
  public:
    cbessel_j_e_functor(N n, const T &x, U &e)
        : functor_foreach_e<cbessel_j_e_functor<order, T, U, N>,
                            T,
                            U,
                            double>(x, e)
        , m_n(n)
    {
    }

    double foreach_e_impl(double x, double &e) const
    {
        return cbessel_j_e_impl<order, N>(m_n, x, e);
    }

  private:
    N m_n;
};

template <typename T, typename U>
inline CwiseNullaryOp<cbessel_j_e_functor<J_n, T, U, int>,
                      typename cbessel_j_e_functor<J_n, T, U, int>::ResultType>
cbessel_j(int n, const DenseBase<T> &x, DenseBase<U> &e)
{
    using ResultType = typename cbessel_j_e_functor<J_n, T, U, int>::ResultType;
    return ResultType::
        NullaryExpr(x.rows(),
                    x.cols(),
                    cbessel_j_e_functor<J_n, T, U, int>(n,
                                                        x.derived(),
                                                        e.derived()));
}

template <typename T, typename U>
inline CwiseNullaryOp<
    cbessel_j_e_functor<J_n, T, U, double>,
    typename cbessel_j_e_functor<J_n, T, U, double>::ResultType>
cbessel_j(double n, const DenseBase<T> &x, DenseBase<U> &e)
{
    using ResultType =
        typename cbessel_j_e_functor<J_n, T, U, double>::ResultType;
    return ResultType::
        NullaryExpr(x.rows(),
                    x.cols(),
                    cbessel_j_e_functor<J_n, T, U, double>(n,
                                                           x.derived(),
                                                           e.derived()));
}

template <enum cbessel_j order, typename T, typename U>
inline CwiseNullaryOp<
    cbessel_j_e_functor<order, T, U, int>,
    typename cbessel_j_e_functor<order, T, U, int>::ResultType>
cbessel_j(const DenseBase<T> &x, DenseBase<U> &e)
{
    static_assert(order < J_n, "order can only be J0/J1/J2");

    using ResultType =
        typename cbessel_j_e_functor<order, T, U, int>::ResultType;
    return ResultType::
        NullaryExpr(x.rows(),
                    x.cols(),
                    cbessel_j_e_functor<order, T, U, int>(order,
                                                          x.derived(),
                                                          e.derived()));
}

// ========================================
// bessel second kind
// ========================================

template <enum cbessel_y order, typename N>
inline double cbessel_y_impl(N n, double x)
{
    throw std::invalid_argument("invalid template parameters");
}

template <>
inline double cbessel_y_impl<Y0, int>(int n, double x)
{
    return gsl_sf_bessel_Y0(x);
}

template <>
inline double cbessel_y_impl<Y1, int>(int n, double x)
{
    return gsl_sf_bessel_Y1(x);
}

template <>
inline double cbessel_y_impl<Y_n, int>(int n, double x)
{
    return gsl_sf_bessel_Yn(n, x);
}

template <>
inline double cbessel_y_impl<Y_n, double>(double n, double x)
{
    return gsl_sf_bessel_Ynu(n, x);
}

template <enum cbessel_y order, typename T, typename N>
class cbessel_y_functor
    : public functor_foreach<cbessel_y_functor<order, T, N>, T, double>
{
  public:
    cbessel_y_functor(N n, const T &x)
        : functor_foreach<cbessel_y_functor<order, T, N>, T, double>(x)
        , m_n(n)
    {
    }

    double foreach_impl(double x) const
    {
        return cbessel_y_impl<order, N>(m_n, x);
    }

  private:
    N m_n;
};

template <typename T>
inline CwiseNullaryOp<cbessel_y_functor<Y_n, T, int>,
                      typename cbessel_y_functor<Y_n, T, int>::ResultType>
cbessel_y(int n, const DenseBase<T> &x)
{
    using ResultType = typename cbessel_y_functor<Y_n, T, int>::ResultType;
    return ResultType::NullaryExpr(x.rows(),
                                   x.cols(),
                                   cbessel_y_functor<Y_n, T, int>(n,
                                                                  x.derived()));
}

template <typename T>
inline CwiseNullaryOp<cbessel_y_functor<Y_n, T, double>,
                      typename cbessel_y_functor<Y_n, T, double>::ResultType>
cbessel_y(double n, const DenseBase<T> &x)
{
    using ResultType = typename cbessel_y_functor<Y_n, T, double>::ResultType;
    return ResultType::
        NullaryExpr(x.rows(),
                    x.cols(),
                    cbessel_y_functor<Y_n, T, double>(n, x.derived()));
}

template <enum cbessel_y order, typename T>
inline CwiseNullaryOp<cbessel_y_functor<order, T, int>,
                      typename cbessel_y_functor<order, T, int>::ResultType>
cbessel_y(const DenseBase<T> &x)
{
    static_assert(order < Y_n, "order can only be Y0/Y1/Y2");

    using ResultType = typename cbessel_y_functor<order, T, int>::ResultType;
    return ResultType::
        NullaryExpr(x.rows(),
                    x.cols(),
                    cbessel_y_functor<order, T, int>(order, x.derived()));
}

template <enum cbessel_y order, typename N>
inline double cbessel_y_e_impl(N n, double x, double &e)
{
    throw std::invalid_argument("invalid template parameters");
}

template <>
inline double cbessel_y_e_impl<Y0, int>(int n, double x, double &e)
{
    gsl_sf_result r;
    gsl_sf_bessel_Y0_e(x, &r);
    e = r.err;
    return r.val;
}

template <>
inline double cbessel_y_e_impl<Y1, int>(int n, double x, double &e)
{
    gsl_sf_result r;
    gsl_sf_bessel_Y1_e(x, &r);
    e = r.err;
    return r.val;
}

template <>
inline double cbessel_y_e_impl<Y_n, int>(int n, double x, double &e)
{
    gsl_sf_result r;
    gsl_sf_bessel_Yn_e(n, x, &r);
    e = r.err;
    return r.val;
}

template <>
inline double cbessel_y_e_impl<Y_n, double>(double n, double x, double &e)
{
    gsl_sf_result r;
    gsl_sf_bessel_Ynu_e(n, x, &r);
    e = r.err;
    return r.val;
}

template <enum cbessel_y order, typename T, typename U, typename N>
class cbessel_y_e_functor
    : public functor_foreach_e<cbessel_y_e_functor<order, T, U, N>,
                               T,
                               U,
                               double>
{
  public:
    cbessel_y_e_functor(N n, const T &x, U &e)
        : functor_foreach_e<cbessel_y_e_functor<order, T, U, N>,
                            T,
                            U,
                            double>(x, e)
        , m_n(n)
    {
    }

    double foreach_e_impl(double x, double &e) const
    {
        return cbessel_y_e_impl<order, N>(m_n, x, e);
    }

  private:
    N m_n;
};

template <typename T, typename U>
inline CwiseNullaryOp<cbessel_y_e_functor<Y_n, T, U, int>,
                      typename cbessel_y_e_functor<Y_n, T, U, int>::ResultType>
cbessel_y(int n, const DenseBase<T> &x, DenseBase<U> &e)
{
    using ResultType = typename cbessel_y_e_functor<Y_n, T, U, int>::ResultType;
    return ResultType::
        NullaryExpr(x.rows(),
                    x.cols(),
                    cbessel_y_e_functor<Y_n, T, U, int>(n,
                                                        x.derived(),
                                                        e.derived()));
}

template <typename T, typename U>
inline CwiseNullaryOp<
    cbessel_y_e_functor<Y_n, T, U, double>,
    typename cbessel_y_e_functor<Y_n, T, U, double>::ResultType>
cbessel_y(double n, const DenseBase<T> &x, DenseBase<U> &e)
{
    using ResultType =
        typename cbessel_y_e_functor<Y_n, T, U, double>::ResultType;
    return ResultType::
        NullaryExpr(x.rows(),
                    x.cols(),
                    cbessel_y_e_functor<Y_n, T, U, double>(n,
                                                           x.derived(),
                                                           e.derived()));
}

template <enum cbessel_y order, typename T, typename U>
inline CwiseNullaryOp<
    cbessel_y_e_functor<order, T, U, int>,
    typename cbessel_y_e_functor<order, T, U, int>::ResultType>
cbessel_y(const DenseBase<T> &x, DenseBase<U> &e)
{
    static_assert(order < Y_n, "order can only be Y0/Y1/Y2");

    using ResultType =
        typename cbessel_y_e_functor<order, T, U, int>::ResultType;
    return ResultType::
        NullaryExpr(x.rows(),
                    x.cols(),
                    cbessel_y_e_functor<order, T, U, int>(order,
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

#endif /* __IEXP_CYLINDRICAL_BESSEL__ */
