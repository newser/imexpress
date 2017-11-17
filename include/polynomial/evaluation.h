/* Copyright (C) 2017 haniu (niuhao.cn@gmail.com)
 *
 * Licensed to the Apache Software Foundation (ASF) under one or more
 * contributor license agreements.  See the NOTICE file distributed with
 * this work for additional information regarding copyright ownership.
 * The ASF licenses this file to You under the Apache License, Version 2.0
 * (the "License"); you may not use this file except in compliance with
 * the License.  You may obtain a copy of the License at
 *
 *     http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */

#ifndef __IEXP_POLY_EVALUATION__
#define __IEXP_POLY_EVALUATION__

////////////////////////////////////////////////////////////
// import header files
////////////////////////////////////////////////////////////

#include <common/common.h>

#include <config.h>
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
    throw std::invalid_argument("todo");
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
inline typename T::Scalar eval(const Eigen::ArrayBase<T> &c,
                               const typename T::Scalar x)
{
    eigen_assert((c.derived().cols() == 1) || (c.derived().rows() == 1));

    typename Eigen::internal::eval<T>::type m_c(c.eval());
    return eval_impl<typename T::Scalar>(m_c.data(), m_c.size(), x);
}

// ========================================
// derivative evaluation
// ========================================

template <typename T>
inline void eval_deriv_impl(
    const T *c, const int len, const T x, T *res, const int res_len)
{
    throw std::invalid_argument("todo");
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
                  T::MaxSizeAtCompileTime,
                  1>
        ArrayType;

    eval_deriv_functor(const T &c, typename T::Scalar x, int order)
        : m_c(c.eval())
        , m_result(order)
    {
        eval_deriv_impl<typename T::Scalar>(m_c.data(),
                                            m_c.size(),
                                            x,
                                            m_result.data(),
                                            order);
    }

    const typename T::Scalar &operator()(Index i) const
    {
        return m_result(i);
    }

  private:
    typename Eigen::internal::eval<T>::type m_c;
    ArrayType m_result;
};

template <typename T>
CwiseNullaryOp<eval_deriv_functor<T>, typename eval_deriv_functor<T>::ArrayType>
eval_deriv(const Eigen::ArrayBase<T> &c, typename T::Scalar x, int order)
{
    eigen_assert((c.derived().cols() == 1) || (c.derived().rows() == 1));

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
