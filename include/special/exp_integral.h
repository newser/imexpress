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

#ifndef __IEXP_EXP_INTEGRAL__
#define __IEXP_EXP_INTEGRAL__

////////////////////////////////////////////////////////////
// import header files
////////////////////////////////////////////////////////////

#include <common/common.h>

#include <common/template_foreach.h>

#include <gsl/gsl_errno.h>
#include <gsl/gsl_sf_expint.h>

IEXP_NS_BEGIN

namespace sf {

////////////////////////////////////////////////////////////
// macro definition
////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////
// type definition
////////////////////////////////////////////////////////////

// ========================================
// expint_1
// ========================================

DEFINE_TEMPLATE_FOREACH(expint_1, double, double, gsl_sf_expint_E1)

DEFINE_TEMPLATE_FOREACH_E(expint_1, double, double, gsl_sf_expint_E1_e)

// ========================================
// expint_2
// ========================================

DEFINE_TEMPLATE_FOREACH(expint_2, double, double, gsl_sf_expint_E2)

DEFINE_TEMPLATE_FOREACH_E(expint_2, double, double, gsl_sf_expint_E2_e)

// ========================================
// expint
// ========================================

template <typename T>
class expint_functor : public functor_foreach<expint_functor<T>, T, double>
{
  public:
    expint_functor(int n, const T &x)
        : functor_foreach<expint_functor<T>, T, double>(x)
        , m_n(n)
    {
    }

    double foreach_impl(double x) const
    {
        return gsl_sf_expint_En(m_n, x);
    }

  private:
    int m_n;
};

template <typename T>
inline CwiseNullaryOp<expint_functor<T>, typename expint_functor<T>::ResultType>
expint(int n, const DenseBase<T> &x)
{
    using ResultType = typename expint_functor<T>::ResultType;
    return ResultType::NullaryExpr(x.rows(),
                                   x.cols(),
                                   expint_functor<T>(n, x.derived()));
}

template <typename T, typename U>
class expint_e_functor
    : public functor_foreach_e<expint_e_functor<T, U>, T, U, double>
{
  public:
    expint_e_functor(int n, const T &x, U &e)
        : functor_foreach_e<expint_e_functor<T, U>, T, U, double>(x, e)
        , m_n(n)
    {
    }

    double foreach_e_impl(double x, double &e) const
    {
        gsl_sf_result r;
        gsl_sf_expint_En_e(m_n, x, &r);
        e = r.err;
        return r.val;
    }

  private:
    int m_n;
};

template <typename T, typename U>
inline CwiseNullaryOp<expint_e_functor<T, U>,
                      typename expint_e_functor<T, U>::ResultType>
expint(int n, const DenseBase<T> &x, DenseBase<U> &e)
{
    using ResultType = typename expint_e_functor<T, U>::ResultType;
    return ResultType::NullaryExpr(x.rows(),
                                   x.cols(),
                                   expint_e_functor<T, U>(n,
                                                          x.derived(),
                                                          e.derived()));
}

// ========================================
// expint_ei
// ========================================

DEFINE_TEMPLATE_FOREACH(expint_ei, double, double, gsl_sf_expint_Ei)

DEFINE_TEMPLATE_FOREACH_E(expint_ei, double, double, gsl_sf_expint_Ei_e)

// ========================================
// expint_ei3
// ========================================

DEFINE_TEMPLATE_FOREACH(expint_ei3, double, double, gsl_sf_expint_3)

DEFINE_TEMPLATE_FOREACH_E(expint_ei3, double, double, gsl_sf_expint_3_e)

// ========================================
// sinh integral
// ========================================

DEFINE_TEMPLATE_FOREACH(sinhint, double, double, gsl_sf_Shi)

DEFINE_TEMPLATE_FOREACH_E(sinhint, double, double, gsl_sf_Shi_e)

// ========================================
// cosh integral
// ========================================

DEFINE_TEMPLATE_FOREACH(coshint, double, double, gsl_sf_Chi)

DEFINE_TEMPLATE_FOREACH_E(coshint, double, double, gsl_sf_Chi_e)

// ========================================
// sin integral
// ========================================

DEFINE_TEMPLATE_FOREACH(sinint, double, double, gsl_sf_Si)

DEFINE_TEMPLATE_FOREACH_E(sinint, double, double, gsl_sf_Si_e)

// ========================================
// cos integral
// ========================================

DEFINE_TEMPLATE_FOREACH(cosint, double, double, gsl_sf_Ci)

DEFINE_TEMPLATE_FOREACH_E(cosint, double, double, gsl_sf_Ci_e)

// ========================================
// atan integral
// ========================================

DEFINE_TEMPLATE_FOREACH(atanint, double, double, gsl_sf_atanint)

DEFINE_TEMPLATE_FOREACH_E(atanint, double, double, gsl_sf_atanint_e)

////////////////////////////////////////////////////////////
// global variants
////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////
// interface declaration
////////////////////////////////////////////////////////////
}

IEXP_NS_END

#endif /* __IEXP_EXP_INTEGRAL__ */
