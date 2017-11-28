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

#ifndef __IEXP_POLY_EVALUATION__
#define __IEXP_POLY_EVALUATION__

////////////////////////////////////////////////////////////
// import header files
////////////////////////////////////////////////////////////

#include <common/common.h>

#include <gsl/gsl_poly.h>

IEXP_NS_BEGIN

namespace poly {

////////////////////////////////////////////////////////////
// macro definition
////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////
// type definition
////////////////////////////////////////////////////////////

// ========================================
// evaluation
// ========================================

template <typename T>
inline T eval_impl(const T *c, const int len, const T x)
{
    UNSUPPORTED_TYPE(T);
}

template <>
inline double eval_impl(const double *c, const int len, const double x)
{
    return gsl_poly_eval(c, len, x);
}

template <>
inline std::complex<double> eval_impl(const std::complex<double> *c,
                                      const int len,
                                      const std::complex<double> x)
{
    gsl_complex ans = gsl_complex_poly_complex_eval((gsl_complex *)c,
                                                    len,
                                                    *(gsl_complex *)&x);
    return std::complex<double>(ans.dat[0], ans.dat[1]);
}

template <typename T>
inline typename T::Scalar eval(const ArrayBase<T> &c,
                               const typename T::Scalar x)
{
    eigen_assert(IS_VEC(c));

    typename type_eval<T>::type m_c(c.eval());
    return eval_impl(m_c.data(), m_c.size(), x);
}

// ========================================
// derivative evaluation
// ========================================

template <typename T>
inline void eval_deriv_impl(
    const T *c, const int len, const T x, T *res, const int res_len)
{
    UNSUPPORTED_TYPE(T);
}

template <>
inline void eval_deriv_impl(const double *c,
                            const int len,
                            const double x,
                            double *res,
                            const int res_len)
{
    gsl_poly_eval_derivs(c, len, x, res, res_len);
}

template <typename T>
class eval_deriv_functor
{
  public:
    typedef Array<typename T::Scalar,
                  T::SizeAtCompileTime,
                  1,
                  ColMajor,
                  T::SizeAtCompileTime,
                  1>
        ArrayType;

    eval_deriv_functor(const T &c, typename T::Scalar x, int order)
        : m_result(order, 1)
    {
        typename type_eval<T>::type m_c(c.eval());
        eval_deriv_impl(m_c.data(), m_c.size(), x, m_result.data(), order);
    }

    const typename T::Scalar &operator()(Index i) const
    {
        return m_result(i);
    }

  private:
    ArrayType m_result;
};

template <typename T>
inline CwiseNullaryOp<eval_deriv_functor<T>,
                      typename eval_deriv_functor<T>::ArrayType>
eval_deriv(const ArrayBase<T> &c, typename T::Scalar x, int order)
{
    eigen_assert(IS_VEC(c));

    // order 0: [f(x)]
    // order 1: [f(x), f'(x)]
    // ...
    // order n generate (n + 1) evaluations
    eigen_assert(order >= 0);
    ++order;

    typedef typename eval_deriv_functor<T>::ArrayType ArrayType;
    return ArrayType::NullaryExpr(order,
                                  eval_deriv_functor<T>(c.derived(), x, order));
}

////////////////////////////////////////////////////////////
// global variants
////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////
// interface declaration
////////////////////////////////////////////////////////////
}

IEXP_NS_END

#endif /* __IEXP_POLY_EVALUATION__ */
