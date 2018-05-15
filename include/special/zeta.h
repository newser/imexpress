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

#ifndef __IEXP_ZETA__
#define __IEXP_ZETA__

////////////////////////////////////////////////////////////
// import header files
////////////////////////////////////////////////////////////

#include <common/common.h>

#include <gsl/gsl_errno.h>
#include <gsl/gsl_sf_zeta.h>

IEXP_NS_BEGIN

namespace sf {

////////////////////////////////////////////////////////////
// macro definition
////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////
// type definition
////////////////////////////////////////////////////////////

// ========================================
// Riemann Zeta
// ========================================

template <typename T>
class zeta_functor : public functor_foreach<zeta_functor<T>, T, double>
{
  public:
    zeta_functor(const T &x)
        : functor_foreach<zeta_functor<T>, T, double>(x)
    {
    }

    double foreach_impl(double x) const
    {
        return gsl_sf_zeta(x);
    }

    double foreach_impl(int x) const
    {
        return gsl_sf_zeta_int(x);
    }
};

template <typename T>
inline CwiseNullaryOp<zeta_functor<T>, typename zeta_functor<T>::ResultType>
zeta(const DenseBase<T> &x)
{
    using ResultType = typename zeta_functor<T>::ResultType;
    return ResultType::NullaryExpr(x.rows(),
                                   x.cols(),
                                   zeta_functor<T>(x.derived()));
}

template <typename T, typename U>
class zeta_e_functor
    : public functor_foreach_e<zeta_e_functor<T, U>, T, U, double>
{
  public:
    zeta_e_functor(const T &x, U &e)
        : functor_foreach_e<zeta_e_functor<T, U>, T, U, double>(x, e)
    {
    }

    double foreach_e_impl(double x, double &e) const
    {
        gsl_sf_result r;
        gsl_sf_zeta_e(x, &r);
        e = r.err;
        return r.val;
    }

    double foreach_e_impl(int x, double &e) const
    {
        gsl_sf_result r;
        gsl_sf_zeta_int_e(x, &r);
        e = r.err;
        return r.val;
    }
};

template <typename T, typename U>
inline CwiseNullaryOp<zeta_e_functor<T, U>,
                      typename zeta_e_functor<T, U>::ResultType>
zeta(const DenseBase<T> &x, DenseBase<U> &e)
{
    using ResultType = typename zeta_e_functor<T, U>::ResultType;
    return ResultType::NullaryExpr(x.rows(),
                                   x.cols(),
                                   zeta_e_functor<T, U>(x.derived(),
                                                        e.derived()));
}

// ========================================
// Riemann Zeta Minus One
// ========================================

template <typename T>
class zetam1_functor : public functor_foreach<zetam1_functor<T>, T, double>
{
  public:
    zetam1_functor(const T &x)
        : functor_foreach<zetam1_functor<T>, T, double>(x)
    {
    }

    double foreach_impl(double x) const
    {
        return gsl_sf_zetam1(x);
    }

    double foreach_impl(int x) const
    {
        return gsl_sf_zetam1_int(x);
    }
};

template <typename T>
inline CwiseNullaryOp<zetam1_functor<T>, typename zetam1_functor<T>::ResultType>
zetam1(const DenseBase<T> &x)
{
    using ResultType = typename zetam1_functor<T>::ResultType;
    return ResultType::NullaryExpr(x.rows(),
                                   x.cols(),
                                   zetam1_functor<T>(x.derived()));
}

template <typename T, typename U>
class zetam1_e_functor
    : public functor_foreach_e<zetam1_e_functor<T, U>, T, U, double>
{
  public:
    zetam1_e_functor(const T &x, U &e)
        : functor_foreach_e<zetam1_e_functor<T, U>, T, U, double>(x, e)
    {
    }

    double foreach_e_impl(double x, double &e) const
    {
        gsl_sf_result r;
        gsl_sf_zetam1_e(x, &r);
        e = r.err;
        return r.val;
    }

    double foreach_e_impl(int x, double &e) const
    {
        gsl_sf_result r;
        gsl_sf_zetam1_int_e(x, &r);
        e = r.err;
        return r.val;
    }
};

template <typename T, typename U>
inline CwiseNullaryOp<zetam1_e_functor<T, U>,
                      typename zetam1_e_functor<T, U>::ResultType>
zetam1(const DenseBase<T> &x, DenseBase<U> &e)
{
    using ResultType = typename zetam1_e_functor<T, U>::ResultType;
    return ResultType::NullaryExpr(x.rows(),
                                   x.cols(),
                                   zetam1_e_functor<T, U>(x.derived(),
                                                          e.derived()));
}

// ========================================
// Eta
// ========================================

template <typename T>
class eta_functor : public functor_foreach<eta_functor<T>, T, double>
{
  public:
    eta_functor(const T &x)
        : functor_foreach<eta_functor<T>, T, double>(x)
    {
    }

    double foreach_impl(double x) const
    {
        return gsl_sf_eta(x);
    }

    double foreach_impl(int x) const
    {
        return gsl_sf_eta_int(x);
    }
};

template <typename T>
inline CwiseNullaryOp<eta_functor<T>, typename eta_functor<T>::ResultType> eta(
    const DenseBase<T> &x)
{
    using ResultType = typename eta_functor<T>::ResultType;
    return ResultType::NullaryExpr(x.rows(),
                                   x.cols(),
                                   eta_functor<T>(x.derived()));
}

template <typename T, typename U>
class eta_e_functor
    : public functor_foreach_e<eta_e_functor<T, U>, T, U, double>
{
  public:
    eta_e_functor(const T &x, U &e)
        : functor_foreach_e<eta_e_functor<T, U>, T, U, double>(x, e)
    {
    }

    double foreach_e_impl(double x, double &e) const
    {
        gsl_sf_result r;
        gsl_sf_eta_e(x, &r);
        e = r.err;
        return r.val;
    }

    double foreach_e_impl(int x, double &e) const
    {
        gsl_sf_result r;
        gsl_sf_eta_int_e(x, &r);
        e = r.err;
        return r.val;
    }
};

template <typename T, typename U>
inline CwiseNullaryOp<eta_e_functor<T, U>,
                      typename eta_e_functor<T, U>::ResultType>
eta(const DenseBase<T> &x, DenseBase<U> &e)
{
    using ResultType = typename eta_e_functor<T, U>::ResultType;
    return ResultType::NullaryExpr(x.rows(),
                                   x.cols(),
                                   eta_e_functor<T, U>(x.derived(),
                                                       e.derived()));
}

// ========================================
// Hurwitz Zeta
// ========================================

template <typename T>
class hzeta_functor : public functor_m2vnum_2d<hzeta_functor<T>, T, double>
{
  public:
    hzeta_functor(const T &s_q)
        : functor_m2vnum_2d<hzeta_functor<T>, T, double>(s_q)
    {
    }

    double m2vnum_impl(double s, double q) const
    {
        gsl_sf_result r;
        gsl_sf_hzeta_e(s, q, &r);
        return r.val;
    }
};

template <typename T>
inline CwiseNullaryOp<hzeta_functor<T>, typename hzeta_functor<T>::ResultType>
hzeta(const DenseBase<T> &s_q)
{
    using ResultType = typename hzeta_functor<T>::ResultType;
    return ResultType::NullaryExpr(M2VNUM_ROW(T, s_q),
                                   M2VNUM_COL(T, s_q),
                                   hzeta_functor<T>(s_q.derived()));
}

template <typename T, typename U>
class hzeta_e_functor
    : public functor_m2vnum_2d_e<hzeta_e_functor<T, U>, T, U, double>
{
  public:
    hzeta_e_functor(const T &s_q, U &e)
        : functor_m2vnum_2d_e<hzeta_e_functor<T, U>, T, U, double>(s_q, e)
    {
    }

    double m2vnum_e_impl(double s, double q, double &e) const
    {
        gsl_sf_result r;
        gsl_sf_hzeta_e(s, q, &r);
        e = r.err;
        return r.val;
    }
};

template <typename T, typename U>
inline CwiseNullaryOp<hzeta_e_functor<T, U>,
                      typename hzeta_e_functor<T, U>::ResultType>
hzeta(const DenseBase<T> &s_q, DenseBase<U> &e)
{
    using ResultType = typename hzeta_e_functor<T, U>::ResultType;
    return ResultType::NullaryExpr(M2VNUM_ROW(T, s_q),
                                   M2VNUM_COL(T, s_q),
                                   hzeta_e_functor<T, U>(s_q.derived(),
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

#endif /* __IEXP_ZETA__ */
