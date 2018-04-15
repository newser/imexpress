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

template <typename T, typename U>
inline U eval_impl(const T *c, int len, U x)
{
    UNSUPPORTED_TYPE(T);
}

template <>
inline double eval_impl(const double *c, int len, double x)
{
    return gsl_poly_eval(c, len, x);
}

template <>
inline std::complex<double> eval_impl(const double *c,
                                      int len,
                                      std::complex<double> x)
{
    gsl_complex ans = gsl_poly_complex_eval(c, len, *(gsl_complex *)&x);
    return std::complex<double>(ans.dat[0], ans.dat[1]);
}

template <>
inline std::complex<double> eval_impl(const std::complex<double> *c,
                                      int len,
                                      std::complex<double> x)
{
    gsl_complex ans = gsl_complex_poly_complex_eval((gsl_complex *)c,
                                                    len,
                                                    *(gsl_complex *)&x);
    return std::complex<double>(ans.dat[0], ans.dat[1]);
}

template <typename T>
inline typename T::Scalar eval(const DenseBase<T> &c, typename T::Scalar x)
{
    typename type_eval<T>::type m_c(c.eval());
    return eval_impl(m_c.data(), m_c.size(), x);
}

template <typename T, typename U>
class eval_functor
{
  public:
    // note scalar should be same with U(input)
    using Scalar = typename U::Scalar;
    using ResultType = typename dense_derive<U>::type;

    eval_functor(const T &c, const U &x)
        : m_result(new Scalar[x.size()])
    {
        // doing evaluating m_c each time calling operator() would be expensive
        // when c is an composed expression(x+y) but not a matrix or an array
        typename type_eval<T>::type m_c(c.eval());
        ASSIGN_IJ(TP4(U) == RowMajor,
                  m_result.get(),
                  x,
                  eval_impl(m_c.data(), m_c.size(), x(i, j)));
    }

    Scalar operator()(Index i) const
    {
        return m_result.get()[i];
    }

  private:
    std::shared_ptr<Scalar> m_result;
};

template <typename T, typename U>
inline CwiseNullaryOp<eval_functor<T, U>,
                      typename eval_functor<T, U>::ResultType>
eval(const DenseBase<T> &c, const DenseBase<U> &x)
{
    using ResultType = typename eval_functor<T, U>::ResultType;
    return ResultType::NullaryExpr(x.rows(),
                                   x.cols(),
                                   eval_functor<T, U>(c.derived(),
                                                      x.derived()));
}

// ========================================
// derivative evaluation
// ========================================

template <typename T>
inline void eval_deriv_impl(const T *c, int len, T x, T *res, int res_len)
{
    UNSUPPORTED_TYPE(T);
}

template <>
inline void eval_deriv_impl(
    const double *c, int len, double x, double *res, int res_len)
{
    gsl_poly_eval_derivs(c, len, x, res, res_len);
}

template <bool row_form, typename T>
class eval_deriv_functor
{
  public:
    using Scalar = typename T::Scalar;
    using ResultType = typename dense_derive<T,
                                             Scalar,
                                             row_form ? 1 : Dynamic,
                                             row_form ? Dynamic : 1>::type;

    eval_deriv_functor(const T &c, Scalar x, int order)
        : m_result(new Scalar[order])
    {
        typename type_eval<T>::type m_c(c.eval());
        eval_deriv_impl(m_c.data(), m_c.size(), x, m_result.get(), order);
    }

    Scalar operator()(Index i) const
    {
        return m_result.get()[i];
    }

  private:
    std::shared_ptr<Scalar> m_result;
};

template <bool row_form = false, typename T = void>
inline CwiseNullaryOp<eval_deriv_functor<row_form, T>,
                      typename eval_deriv_functor<row_form, T>::ResultType>
eval_deriv(const DenseBase<T> &c, typename T::Scalar x, int order)
{
    // order 0: [f(x)]
    // order 1: [f(x), f'(x)]
    // ...
    // order n generate (n + 1) evaluations
    eigen_assert(order >= 0);
    ++order;

    using ResultType = typename eval_deriv_functor<row_form, T>::ResultType;
    return ResultType::NullaryExpr(row_form ? 1 : order,
                                   row_form ? order : 1,
                                   eval_deriv_functor<row_form, T>(c.derived(),
                                                                   x,
                                                                   order));
}

template <bool row_form, typename T, typename U>
class eval2_deriv_functor
{
  public:
    // note scalar should be same with U(input)
    using Scalar = typename U::Scalar;
    using ResultType = typename dense_derive<T,
                                             Scalar,
                                             Dynamic,
                                             Dynamic,
                                             row_form ? RowMajor : 0,
                                             Dynamic,
                                             Dynamic>::type;

    eval2_deriv_functor(const T &c, const U &x, int order)
        : m_result(row_form ? x.size() : order, row_form ? order : x.size())
    {
        typename type_eval<T>::type m_c(c.eval());
        for (int i = 0; i < x.size(); ++i) {
            eval_deriv_impl(m_c.data(),
                            m_c.size(),
                            x.data()[i],
                            &m_result.data()[i * order],
                            order);
        }
    }

    Scalar operator()(Index i, Index j) const
    {
        return m_result(i, j);
    }

  private:
    buf_rc<Scalar, row_form> m_result;
};

template <bool row_form = false, typename T = void, typename U = void>
inline CwiseNullaryOp<eval2_deriv_functor<row_form, T, U>,
                      typename eval2_deriv_functor<row_form, T, U>::ResultType>
eval_deriv(const DenseBase<T> &c, const DenseBase<U> &x, int order)
{
    eigen_assert(order >= 0);
    ++order;

    using ResultType = typename eval2_deriv_functor<row_form, T, U>::ResultType;
    return ResultType::
        NullaryExpr(row_form ? x.size() : order,
                    row_form ? order : x.size(),
                    eval2_deriv_functor<row_form, T, U>(c.derived(), x, order));
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
