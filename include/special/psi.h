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

#ifndef __IEXP_SF_PSI__
#define __IEXP_SF_PSI__

////////////////////////////////////////////////////////////
// import header files
////////////////////////////////////////////////////////////

#include <common/common.h>

#include <common/template_foreach.h>

#include <gsl/gsl_errno.h>
#include <gsl/gsl_sf_psi.h>

IEXP_NS_BEGIN

namespace sf {

////////////////////////////////////////////////////////////
// macro definition
////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////
// type definition
////////////////////////////////////////////////////////////

// ========================================
// digamma
// ========================================

template <typename T>
class psi_functor : public functor_foreach<psi_functor<T>, T, double>
{
  public:
    psi_functor(const T &x)
        : functor_foreach<psi_functor<T>, T, double>(x)
    {
    }

    double foreach_impl(double x) const
    {
        return gsl_sf_psi(x);
    }

    double foreach_impl(int x) const
    {
        return gsl_sf_psi_int(x);
    }
};

template <typename T>
inline CwiseNullaryOp<psi_functor<T>, typename psi_functor<T>::ResultType> psi(
    const DenseBase<T> &x)
{
    using ResultType = typename psi_functor<T>::ResultType;
    return ResultType::NullaryExpr(x.rows(),
                                   x.cols(),
                                   psi_functor<T>(x.derived()));
}

template <typename T, typename U>
class psi_e_functor
    : public functor_foreach_e<psi_e_functor<T, U>, T, U, double>
{
  public:
    psi_e_functor(const T &x, U &e)
        : functor_foreach_e<psi_e_functor<T, U>, T, U, double>(x, e)
    {
    }

    double foreach_e_impl(double x, double &e) const
    {
        gsl_sf_result r;
        gsl_sf_psi_e(x, &r);
        e = r.err;
        return r.val;
    }

    double foreach_e_impl(int x, double &e) const
    {
        gsl_sf_result r;
        gsl_sf_psi_int_e(x, &r);
        e = r.err;
        return r.val;
    }
};

template <typename T, typename U>
inline CwiseNullaryOp<psi_e_functor<T, U>,
                      typename psi_e_functor<T, U>::ResultType>
psi(const DenseBase<T> &x, DenseBase<U> &e)
{
    using ResultType = typename psi_e_functor<T, U>::ResultType;
    return ResultType::NullaryExpr(x.rows(),
                                   x.cols(),
                                   psi_e_functor<T, U>(x.derived(),
                                                       e.derived()));
}

// ========================================
// digamma(1+ix)
// ========================================

DEFINE_TEMPLATE_FOREACH(psi_1pix, double, double, gsl_sf_psi_1piy)

DEFINE_TEMPLATE_FOREACH_E(psi_1pix, double, double, gsl_sf_psi_1piy_e)

// ========================================
// trigamma
// ========================================

template <typename T>
class psi_d1_functor : public functor_foreach<psi_d1_functor<T>, T, double>
{
  public:
    psi_d1_functor(const T &x)
        : functor_foreach<psi_d1_functor<T>, T, double>(x)
    {
    }

    double foreach_impl(double x) const
    {
        return gsl_sf_psi_1(x);
    }

    double foreach_impl(int x) const
    {
        return gsl_sf_psi_1_int(x);
    }
};

template <typename T>
inline CwiseNullaryOp<psi_d1_functor<T>, typename psi_d1_functor<T>::ResultType>
psi_d1(const DenseBase<T> &x)
{
    using ResultType = typename psi_d1_functor<T>::ResultType;
    return ResultType::NullaryExpr(x.rows(),
                                   x.cols(),
                                   psi_d1_functor<T>(x.derived()));
}

template <typename T, typename U>
class psi_d1_e_functor
    : public functor_foreach_e<psi_d1_e_functor<T, U>, T, U, double>
{
  public:
    psi_d1_e_functor(const T &x, U &e)
        : functor_foreach_e<psi_d1_e_functor<T, U>, T, U, double>(x, e)
    {
    }

    double foreach_e_impl(double x, double &e) const
    {
        gsl_sf_result r;
        gsl_sf_psi_1_e(x, &r);
        e = r.err;
        return r.val;
    }

    double foreach_e_impl(int x, double &e) const
    {
        gsl_sf_result r;
        gsl_sf_psi_1_int_e(x, &r);
        e = r.err;
        return r.val;
    }
};

template <typename T, typename U>
inline CwiseNullaryOp<psi_d1_e_functor<T, U>,
                      typename psi_d1_e_functor<T, U>::ResultType>
psi_d1(const DenseBase<T> &x, DenseBase<U> &e)
{
    using ResultType = typename psi_d1_e_functor<T, U>::ResultType;
    return ResultType::NullaryExpr(x.rows(),
                                   x.cols(),
                                   psi_d1_e_functor<T, U>(x.derived(),
                                                          e.derived()));
}

// ========================================
// polygamma
// ========================================

template <typename T>
class psi_dn_functor : public functor_foreach<psi_dn_functor<T>, T, double>
{
  public:
    psi_dn_functor(int n, const T &x)
        : functor_foreach<psi_dn_functor<T>, T, double>(x)
        , m_n(n)
    {
    }

    double foreach_impl(double x) const
    {
        return gsl_sf_psi_n(m_n, x);
    }

  private:
    int m_n;
};

template <typename T>
inline CwiseNullaryOp<psi_dn_functor<T>, typename psi_dn_functor<T>::ResultType>
psi_dn(int n, const DenseBase<T> &x)
{
    using ResultType = typename psi_dn_functor<T>::ResultType;
    return ResultType::NullaryExpr(x.rows(),
                                   x.cols(),
                                   psi_dn_functor<T>(n, x.derived()));
}

template <typename T, typename U>
class psi_dn_e_functor
    : public functor_foreach_e<psi_dn_e_functor<T, U>, T, U, double>
{
  public:
    psi_dn_e_functor(int n, const T &x, U &e)
        : functor_foreach_e<psi_dn_e_functor<T, U>, T, U, double>(x, e)
        , m_n(n)
    {
    }

    double foreach_e_impl(double x, double &e) const
    {
        gsl_sf_result r;
        gsl_sf_psi_n_e(m_n, x, &r);
        e = r.err;
        return r.val;
    }

  private:
    int m_n;
};

template <typename T, typename U>
inline CwiseNullaryOp<psi_dn_e_functor<T, U>,
                      typename psi_dn_e_functor<T, U>::ResultType>
psi_dn(int n, const DenseBase<T> &x, DenseBase<U> &e)
{
    using ResultType = typename psi_dn_e_functor<T, U>::ResultType;
    return ResultType::NullaryExpr(x.rows(),
                                   x.cols(),
                                   psi_dn_e_functor<T, U>(n,
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

#endif /* __IEXP_SF_PSI__ */
