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

#ifndef __IEXP_MODIFIED_CYLINDRICAL_BESSEL__
#define __IEXP_MODIFIED_CYLINDRICAL_BESSEL__

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

enum cbessel_i
{
    I0,
    I1,
    I_n,
};

enum cbessel_k
{
    K0,
    K1,
    K_n,
};

// ========================================
// modified bessel first kind
// ========================================

template <enum cbessel_i order, bool scaled>
inline double cbessel_i_impl(int n, double x)
{
    throw std::invalid_argument("");
}

template <>
inline double cbessel_i_impl<I0, true>(int n, double x)
{
    return gsl_sf_bessel_I0_scaled(x);
}

template <>
inline double cbessel_i_impl<I0, false>(int n, double x)
{
    return gsl_sf_bessel_I0(x);
}

template <>
inline double cbessel_i_impl<I1, true>(int n, double x)
{
    return gsl_sf_bessel_I1_scaled(x);
}

template <>
inline double cbessel_i_impl<I1, false>(int n, double x)
{
    return gsl_sf_bessel_I1(x);
}

template <>
inline double cbessel_i_impl<I_n, true>(int n, double x)
{
    return gsl_sf_bessel_In_scaled(n, x);
}

template <>
inline double cbessel_i_impl<I_n, false>(int n, double x)
{
    return gsl_sf_bessel_In(n, x);
}

template <enum cbessel_i order, bool scaled, typename T>
class cbessel_i_functor
    : public functor_foreach<cbessel_i_functor<order, scaled, T>, T, double>
{
  public:
    cbessel_i_functor(int n, const T &x)
        : functor_foreach<cbessel_i_functor<order, scaled, T>, T, double>(x)
        , m_n(n)
    {
    }

    double foreach_impl(double x) const
    {
        return cbessel_i_impl<order, scaled>(m_n, x);
    }

  private:
    int m_n;
};

template <bool scaled = true, typename T = void>
inline CwiseNullaryOp<cbessel_i_functor<I_n, scaled, T>,
                      typename cbessel_i_functor<I_n, scaled, T>::ResultType>
cbessel_i(int n, const DenseBase<T> &x)
{
    using ResultType = typename cbessel_i_functor<I_n, scaled, T>::ResultType;
    return ResultType::
        NullaryExpr(x.rows(),
                    x.cols(),
                    cbessel_i_functor<I_n, scaled, T>(n, x.derived()));
}

template <enum cbessel_i order, bool scaled = true, typename T = void>
inline CwiseNullaryOp<cbessel_i_functor<order, scaled, T>,
                      typename cbessel_i_functor<order, scaled, T>::ResultType>
cbessel_i(const DenseBase<T> &x)
{
    using ResultType = typename cbessel_i_functor<order, scaled, T>::ResultType;
    return ResultType::
        NullaryExpr(x.rows(),
                    x.cols(),
                    cbessel_i_functor<order, scaled, T>(order, x.derived()));
}

template <enum cbessel_i order, bool scaled>
inline double cbessel_i_e_impl(int n, double x, double &e)
{
    throw std::invalid_argument("");
}

template <>
inline double cbessel_i_e_impl<I0, true>(int n, double x, double &e)
{
    gsl_sf_result r;
    gsl_sf_bessel_I0_scaled_e(x, &r);
    e = r.err;
    return r.val;
}

template <>
inline double cbessel_i_e_impl<I0, false>(int n, double x, double &e)
{
    gsl_sf_result r;
    gsl_sf_bessel_I0_e(x, &r);
    e = r.err;
    return r.val;
}

template <>
inline double cbessel_i_e_impl<I1, true>(int n, double x, double &e)
{
    gsl_sf_result r;
    gsl_sf_bessel_I1_scaled_e(x, &r);
    e = r.err;
    return r.val;
}

template <>
inline double cbessel_i_e_impl<I1, false>(int n, double x, double &e)
{
    gsl_sf_result r;
    gsl_sf_bessel_I1_e(x, &r);
    e = r.err;
    return r.val;
}

template <>
inline double cbessel_i_e_impl<I_n, true>(int n, double x, double &e)
{
    gsl_sf_result r;
    gsl_sf_bessel_In_scaled_e(n, x, &r);
    e = r.err;
    return r.val;
}

template <>
inline double cbessel_i_e_impl<I_n, false>(int n, double x, double &e)
{
    gsl_sf_result r;
    gsl_sf_bessel_In_e(n, x, &r);
    e = r.err;
    return r.val;
}

template <enum cbessel_i order, bool scaled, typename T, typename U>
class cbessel_i_e_functor
    : public functor_foreach_e<cbessel_i_e_functor<order, scaled, T, U>,
                               T,
                               U,
                               double>
{
  public:
    cbessel_i_e_functor(int n, const T &x, U &e)
        : functor_foreach_e<cbessel_i_e_functor<order, scaled, T, U>,
                            T,
                            U,
                            double>(x, e)
        , m_n(n)
    {
    }

    double foreach_e_impl(double x, double &e) const
    {
        return cbessel_i_e_impl<order, scaled>(m_n, x, e);
    }

  private:
    int m_n;
};

template <bool scaled = true, typename T = void, typename U = void>
inline CwiseNullaryOp<
    cbessel_i_e_functor<I_n, scaled, T, U>,
    typename cbessel_i_e_functor<I_n, scaled, T, U>::ResultType>
cbessel_i(int n, const DenseBase<T> &x, DenseBase<U> &e)
{
    using ResultType =
        typename cbessel_i_e_functor<I_n, scaled, T, U>::ResultType;
    return ResultType::
        NullaryExpr(x.rows(),
                    x.cols(),
                    cbessel_i_e_functor<I_n, scaled, T, U>(n,
                                                           x.derived(),
                                                           e.derived()));
}

template <enum cbessel_i order,
          bool scaled = true,
          typename T = void,
          typename U = void>
inline CwiseNullaryOp<
    cbessel_i_e_functor<order, scaled, T, U>,
    typename cbessel_i_e_functor<order, scaled, T, U>::ResultType>
cbessel_i(const DenseBase<T> &x, DenseBase<U> &e)
{
    using ResultType =
        typename cbessel_i_e_functor<order, scaled, T, U>::ResultType;
    return ResultType::
        NullaryExpr(x.rows(),
                    x.cols(),
                    cbessel_i_e_functor<order, scaled, T, U>(order,
                                                             x.derived(),
                                                             e.derived()));
}

// ========================================
// modified bessel second kind
// ========================================

template <enum cbessel_k order, bool scaled>
inline double cbessel_k_impl(int n, double x)
{
    throw std::invalid_argument("");
}

template <>
inline double cbessel_k_impl<K0, true>(int n, double x)
{
    return gsl_sf_bessel_K0_scaled(x);
}

template <>
inline double cbessel_k_impl<K0, false>(int n, double x)
{
    return gsl_sf_bessel_K0(x);
}

template <>
inline double cbessel_k_impl<K1, true>(int n, double x)
{
    return gsl_sf_bessel_K1_scaled(x);
}

template <>
inline double cbessel_k_impl<K1, false>(int n, double x)
{
    return gsl_sf_bessel_K1(x);
}

template <>
inline double cbessel_k_impl<K_n, true>(int n, double x)
{
    return gsl_sf_bessel_Kn_scaled(n, x);
}

template <>
inline double cbessel_k_impl<K_n, false>(int n, double x)
{
    return gsl_sf_bessel_Kn(n, x);
}

template <enum cbessel_k order, bool scaled, typename T>
class cbessel_k_functor
    : public functor_foreach<cbessel_k_functor<order, scaled, T>, T, double>
{
  public:
    cbessel_k_functor(int n, const T &x)
        : functor_foreach<cbessel_k_functor<order, scaled, T>, T, double>(x)
        , m_n(n)
    {
    }

    double foreach_impl(double x) const
    {
        return cbessel_k_impl<order, scaled>(m_n, x);
    }

  private:
    int m_n;
};

template <bool scaled = true, typename T = void>
inline CwiseNullaryOp<cbessel_k_functor<K_n, scaled, T>,
                      typename cbessel_k_functor<K_n, scaled, T>::ResultType>
cbessel_k(int n, const DenseBase<T> &x)
{
    using ResultType = typename cbessel_k_functor<K_n, scaled, T>::ResultType;
    return ResultType::
        NullaryExpr(x.rows(),
                    x.cols(),
                    cbessel_k_functor<K_n, scaled, T>(n, x.derived()));
}

template <enum cbessel_k order, bool scaled = true, typename T = void>
inline CwiseNullaryOp<cbessel_k_functor<order, scaled, T>,
                      typename cbessel_k_functor<order, scaled, T>::ResultType>
cbessel_k(const DenseBase<T> &x)
{
    using ResultType = typename cbessel_k_functor<order, scaled, T>::ResultType;
    return ResultType::
        NullaryExpr(x.rows(),
                    x.cols(),
                    cbessel_k_functor<order, scaled, T>(order, x.derived()));
}

template <enum cbessel_k order, bool scaled>
inline double cbessel_k_e_impl(int n, double x, double &e)
{
    throw std::invalid_argument("must set template param 'scaled' to true");
}

template <>
inline double cbessel_k_e_impl<K0, true>(int n, double x, double &e)
{
    gsl_sf_result r;
    gsl_sf_bessel_K0_scaled_e(x, &r);
    e = r.err;
    return r.val;
}

template <>
inline double cbessel_k_e_impl<K0, false>(int n, double x, double &e)
{
    gsl_sf_result r;
    gsl_sf_bessel_K0_e(x, &r);
    e = r.err;
    return r.val;
}

template <>
inline double cbessel_k_e_impl<K1, true>(int n, double x, double &e)
{
    gsl_sf_result r;
    gsl_sf_bessel_K1_scaled_e(x, &r);
    e = r.err;
    return r.val;
}

template <>
inline double cbessel_k_e_impl<K1, false>(int n, double x, double &e)
{
    gsl_sf_result r;
    gsl_sf_bessel_K1_e(x, &r);
    e = r.err;
    return r.val;
}

template <>
inline double cbessel_k_e_impl<K_n, true>(int n, double x, double &e)
{
    gsl_sf_result r;
    gsl_sf_bessel_Kn_scaled_e(n, x, &r);
    e = r.err;
    return r.val;
}

template <>
inline double cbessel_k_e_impl<K_n, false>(int n, double x, double &e)
{
    gsl_sf_result r;
    gsl_sf_bessel_Kn_e(n, x, &r);
    e = r.err;
    return r.val;
}

template <enum cbessel_k order, bool scaled, typename T, typename U>
class cbessel_k_e_functor
    : public functor_foreach_e<cbessel_k_e_functor<order, scaled, T, U>,
                               T,
                               U,
                               double>
{
  public:
    cbessel_k_e_functor(int n, const T &x, U &e)
        : functor_foreach_e<cbessel_k_e_functor<order, scaled, T, U>,
                            T,
                            U,
                            double>(x, e)
        , m_n(n)
    {
    }

    double foreach_e_impl(double x, double &e) const
    {
        return cbessel_k_e_impl<order, scaled>(m_n, x, e);
    }

  private:
    int m_n;
};

template <bool scaled = true, typename T = void, typename U = void>
inline CwiseNullaryOp<
    cbessel_k_e_functor<K_n, scaled, T, U>,
    typename cbessel_k_e_functor<K_n, scaled, T, U>::ResultType>
cbessel_k(int n, const DenseBase<T> &x, DenseBase<U> &e)
{
    using ResultType =
        typename cbessel_k_e_functor<K_n, scaled, T, U>::ResultType;
    return ResultType::
        NullaryExpr(x.rows(),
                    x.cols(),
                    cbessel_k_e_functor<K_n, scaled, T, U>(n,
                                                           x.derived(),
                                                           e.derived()));
}

template <enum cbessel_k order,
          bool scaled = true,
          typename T = void,
          typename U = void>
inline CwiseNullaryOp<
    cbessel_k_e_functor<order, scaled, T, U>,
    typename cbessel_k_e_functor<order, scaled, T, U>::ResultType>
cbessel_k(const DenseBase<T> &x, DenseBase<U> &e)
{
    using ResultType =
        typename cbessel_k_e_functor<order, scaled, T, U>::ResultType;
    return ResultType::
        NullaryExpr(x.rows(),
                    x.cols(),
                    cbessel_k_e_functor<order, scaled, T, U>(order,
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

#endif /* __IEXP_MODIFIED_CYLINDRICAL_BESSEL__ */
