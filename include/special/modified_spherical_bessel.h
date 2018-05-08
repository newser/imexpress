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

#ifndef __IEXP_MODIFIED_SPHERICAL_BESSEL__
#define __IEXP_MODIFIED_SPHERICAL_BESSEL__

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

enum sbessel_i
{
    i0,
    i1,
    i2,
    i_n,
};

enum sbessel_k
{
    k0,
    k1,
    k2,
    k_n,
};

// ========================================
// modified bessel first kind
// ========================================

template <enum sbessel_i order, bool scaled>
inline double sbessel_i_impl(int n, double x)
{
    throw std::invalid_argument("must set template param 'scaled' to true");
}

template <>
inline double sbessel_i_impl<i0, true>(int n, double x)
{
    return gsl_sf_bessel_i0_scaled(x);
}

template <>
inline double sbessel_i_impl<i1, true>(int n, double x)
{
    return gsl_sf_bessel_i1_scaled(x);
}

template <>
inline double sbessel_i_impl<i2, true>(int n, double x)
{
    return gsl_sf_bessel_i2_scaled(x);
}

template <>
inline double sbessel_i_impl<i_n, true>(int n, double x)
{
    return gsl_sf_bessel_il_scaled(n, x);
}

template <enum sbessel_i order, bool scaled, typename T>
class sbessel_i_functor
    : public functor_foreach<sbessel_i_functor<order, scaled, T>, T, double>
{
  public:
    sbessel_i_functor(int n, const T &x)
        : functor_foreach<sbessel_i_functor<order, scaled, T>, T, double>(x)
        , m_n(n)
    {
    }

    double foreach_impl(double x) const
    {
        return sbessel_i_impl<order, scaled>(m_n, x);
    }

  private:
    int m_n;
};

template <bool scaled = true, typename T = void>
inline CwiseNullaryOp<sbessel_i_functor<i_n, scaled, T>,
                      typename sbessel_i_functor<i_n, scaled, T>::ResultType>
sbessel_i(int n, const DenseBase<T> &x)
{
    using ResultType = typename sbessel_i_functor<i_n, scaled, T>::ResultType;
    return ResultType::
        NullaryExpr(x.rows(),
                    x.cols(),
                    sbessel_i_functor<i_n, scaled, T>(n, x.derived()));
}

template <enum sbessel_i order, bool scaled = true, typename T = void>
inline CwiseNullaryOp<sbessel_i_functor<order, scaled, T>,
                      typename sbessel_i_functor<order, scaled, T>::ResultType>
sbessel_i(const DenseBase<T> &x)
{
    using ResultType = typename sbessel_i_functor<order, scaled, T>::ResultType;
    return ResultType::
        NullaryExpr(x.rows(),
                    x.cols(),
                    sbessel_i_functor<order, scaled, T>(order, x.derived()));
}

template <enum sbessel_i order, bool scaled>
inline double sbessel_i_e_impl(int n, double x, double &e)
{
    throw std::invalid_argument("must set template param 'scaled' to true");
}

template <>
inline double sbessel_i_e_impl<i0, true>(int n, double x, double &e)
{
    gsl_sf_result r;
    gsl_sf_bessel_i0_scaled_e(x, &r);
    e = r.err;
    return r.val;
}

template <>
inline double sbessel_i_e_impl<i1, true>(int n, double x, double &e)
{
    gsl_sf_result r;
    gsl_sf_bessel_i1_scaled_e(x, &r);
    e = r.err;
    return r.val;
}

template <>
inline double sbessel_i_e_impl<i2, true>(int n, double x, double &e)
{
    gsl_sf_result r;
    gsl_sf_bessel_i2_scaled_e(x, &r);
    e = r.err;
    return r.val;
}

template <>
inline double sbessel_i_e_impl<i_n, true>(int n, double x, double &e)
{
    gsl_sf_result r;
    gsl_sf_bessel_il_scaled_e(n, x, &r);
    e = r.err;
    return r.val;
}

template <enum sbessel_i order, bool scaled, typename T, typename U>
class sbessel_i_e_functor
    : public functor_foreach_e<sbessel_i_e_functor<order, scaled, T, U>,
                               T,
                               U,
                               double>
{
  public:
    sbessel_i_e_functor(int n, const T &x, U &e)
        : functor_foreach_e<sbessel_i_e_functor<order, scaled, T, U>,
                            T,
                            U,
                            double>(x, e)
        , m_n(n)
    {
    }

    double foreach_e_impl(double x, double &e) const
    {
        return sbessel_i_e_impl<order, scaled>(m_n, x, e);
    }

  private:
    int m_n;
};

template <bool scaled = true, typename T = void, typename U = void>
inline CwiseNullaryOp<
    sbessel_i_e_functor<i_n, scaled, T, U>,
    typename sbessel_i_e_functor<i_n, scaled, T, U>::ResultType>
sbessel_i(int n, const DenseBase<T> &x, DenseBase<U> &e)
{
    using ResultType =
        typename sbessel_i_e_functor<i_n, scaled, T, U>::ResultType;
    return ResultType::
        NullaryExpr(x.rows(),
                    x.cols(),
                    sbessel_i_e_functor<i_n, scaled, T, U>(n,
                                                           x.derived(),
                                                           e.derived()));
}

template <enum sbessel_i order,
          bool scaled = true,
          typename T = void,
          typename U = void>
inline CwiseNullaryOp<
    sbessel_i_e_functor<order, scaled, T, U>,
    typename sbessel_i_e_functor<order, scaled, T, U>::ResultType>
sbessel_i(const DenseBase<T> &x, DenseBase<U> &e)
{
    using ResultType =
        typename sbessel_i_e_functor<order, scaled, T, U>::ResultType;
    return ResultType::
        NullaryExpr(x.rows(),
                    x.cols(),
                    sbessel_i_e_functor<order, scaled, T, U>(order,
                                                             x.derived(),
                                                             e.derived()));
}

// ========================================
// modified bessel second kind
// ========================================

template <enum sbessel_k order, bool scaled>
inline double sbessel_k_impl(int n, double x)
{
    throw std::invalid_argument("must set template param 'scaled' to true");
}

template <>
inline double sbessel_k_impl<k0, true>(int n, double x)
{
    return gsl_sf_bessel_k0_scaled(x);
}

template <>
inline double sbessel_k_impl<k1, true>(int n, double x)
{
    return gsl_sf_bessel_k1_scaled(x);
}

template <>
inline double sbessel_k_impl<k2, true>(int n, double x)
{
    return gsl_sf_bessel_k2_scaled(x);
}

template <>
inline double sbessel_k_impl<k_n, true>(int n, double x)
{
    return gsl_sf_bessel_kl_scaled(n, x);
}

template <enum sbessel_k order, bool scaled, typename T>
class sbessel_k_functor
    : public functor_foreach<sbessel_k_functor<order, scaled, T>, T, double>
{
  public:
    sbessel_k_functor(int n, const T &x)
        : functor_foreach<sbessel_k_functor<order, scaled, T>, T, double>(x)
        , m_n(n)
    {
    }

    double foreach_impl(double x) const
    {
        return sbessel_k_impl<order, scaled>(m_n, x);
    }

  private:
    int m_n;
};

template <bool scaled = true, typename T = void>
inline CwiseNullaryOp<sbessel_k_functor<k_n, scaled, T>,
                      typename sbessel_k_functor<k_n, scaled, T>::ResultType>
sbessel_k(int n, const DenseBase<T> &x)
{
    using ResultType = typename sbessel_k_functor<k_n, scaled, T>::ResultType;
    return ResultType::
        NullaryExpr(x.rows(),
                    x.cols(),
                    sbessel_k_functor<k_n, scaled, T>(n, x.derived()));
}

template <enum sbessel_k order, bool scaled = true, typename T = void>
inline CwiseNullaryOp<sbessel_k_functor<order, scaled, T>,
                      typename sbessel_k_functor<order, scaled, T>::ResultType>
sbessel_k(const DenseBase<T> &x)
{
    using ResultType = typename sbessel_k_functor<order, scaled, T>::ResultType;
    return ResultType::
        NullaryExpr(x.rows(),
                    x.cols(),
                    sbessel_k_functor<order, scaled, T>(order, x.derived()));
}

template <enum sbessel_k order, bool scaled>
inline double sbessel_k_e_impl(int n, double x, double &e)
{
    throw std::invalid_argument("must set template param 'scaled' to true");
}

template <>
inline double sbessel_k_e_impl<k0, true>(int n, double x, double &e)
{
    gsl_sf_result r;
    gsl_sf_bessel_k0_scaled_e(x, &r);
    e = r.err;
    return r.val;
}

template <>
inline double sbessel_k_e_impl<k1, true>(int n, double x, double &e)
{
    gsl_sf_result r;
    gsl_sf_bessel_k1_scaled_e(x, &r);
    e = r.err;
    return r.val;
}

template <>
inline double sbessel_k_e_impl<k2, true>(int n, double x, double &e)
{
    gsl_sf_result r;
    gsl_sf_bessel_k2_scaled_e(x, &r);
    e = r.err;
    return r.val;
}

template <>
inline double sbessel_k_e_impl<k_n, true>(int n, double x, double &e)
{
    gsl_sf_result r;
    gsl_sf_bessel_kl_scaled_e(n, x, &r);
    e = r.err;
    return r.val;
}

template <enum sbessel_k order, bool scaled, typename T, typename U>
class sbessel_k_e_functor
    : public functor_foreach_e<sbessel_k_e_functor<order, scaled, T, U>,
                               T,
                               U,
                               double>
{
  public:
    sbessel_k_e_functor(int n, const T &x, U &e)
        : functor_foreach_e<sbessel_k_e_functor<order, scaled, T, U>,
                            T,
                            U,
                            double>(x, e)
        , m_n(n)
    {
    }

    double foreach_e_impl(double x, double &e) const
    {
        return sbessel_k_e_impl<order, scaled>(m_n, x, e);
    }

  private:
    int m_n;
};

template <bool scaled = true, typename T = void, typename U = void>
inline CwiseNullaryOp<
    sbessel_k_e_functor<k_n, scaled, T, U>,
    typename sbessel_k_e_functor<k_n, scaled, T, U>::ResultType>
sbessel_k(int n, const DenseBase<T> &x, DenseBase<U> &e)
{
    using ResultType =
        typename sbessel_k_e_functor<k_n, scaled, T, U>::ResultType;
    return ResultType::
        NullaryExpr(x.rows(),
                    x.cols(),
                    sbessel_k_e_functor<k_n, scaled, T, U>(n,
                                                           x.derived(),
                                                           e.derived()));
}

template <enum sbessel_k order,
          bool scaled = true,
          typename T = void,
          typename U = void>
inline CwiseNullaryOp<
    sbessel_k_e_functor<order, scaled, T, U>,
    typename sbessel_k_e_functor<order, scaled, T, U>::ResultType>
sbessel_k(const DenseBase<T> &x, DenseBase<U> &e)
{
    using ResultType =
        typename sbessel_k_e_functor<order, scaled, T, U>::ResultType;
    return ResultType::
        NullaryExpr(x.rows(),
                    x.cols(),
                    sbessel_k_e_functor<order, scaled, T, U>(order,
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

#endif /* __IEXP_MODIFIED_SPHERICAL_BESSEL__ */
