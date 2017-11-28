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

#ifndef __IEXP_POLY_CUBIC__
#define __IEXP_POLY_CUBIC__

////////////////////////////////////////////////////////////
// import header files
////////////////////////////////////////////////////////////

#include <common/common.h>
#include <math/inf_nan.h>

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
// solve in real field
// ========================================

template <typename T>
inline void solve_cubic_impl(
    const T a, const T b, const T c, double *x0, double *x1, double *x2)
{
    UNSUPPORTED_TYPE(T);
}

template <>
inline void solve_cubic_impl(const double a,
                             const double b,
                             const double c,
                             double *x0,
                             double *x1,
                             double *x2)
{
    gsl_poly_solve_cubic(a, b, c, x0, x1, x2);
}

template <typename T>
class solve_cubic_functor
{
  public:
    typedef Array<double, 3, 1, ColMajor, 3, 1> ArrayType;

    solve_cubic_functor(const T &c)
        : m_result(IEXP_NAN, IEXP_NAN, IEXP_NAN)
    {
        typename type_eval<T>::type m_c(c.eval());
        solve_cubic_impl(m_c[0],
                         m_c[1],
                         m_c[2],
                         &m_result[0],
                         &m_result[1],
                         &m_result[2]);
    }

    const double &operator()(Index i) const
    {
        return m_result(i);
    }

  private:
    ArrayType m_result;
};

template <typename T>
inline CwiseNullaryOp<solve_cubic_functor<T>,
                      typename solve_cubic_functor<T>::ArrayType>
solve_cubic(const ArrayBase<T> &c)
{
    eigen_assert(IS_VEC(c) && (c.size() == 3));

    typedef typename solve_cubic_functor<T>::ArrayType ArrayType;
    return ArrayType::NullaryExpr(3, solve_cubic_functor<T>(c.derived()));
}

// ========================================
// solve in complex field
// ========================================

template <typename T>
inline void complex_solve_cubic_impl(const T a,
                                     const T b,
                                     const T c,
                                     std::complex<double> *x0,
                                     std::complex<double> *x1,
                                     std::complex<double> *x2)
{
    UNSUPPORTED_TYPE(T);
}

template <>
inline void complex_solve_cubic_impl(const double a,
                                     const double b,
                                     const double c,
                                     std::complex<double> *x0,
                                     std::complex<double> *x1,
                                     std::complex<double> *x2)
{
    gsl_poly_complex_solve_cubic(a,
                                 b,
                                 c,
                                 (gsl_complex *)x0,
                                 (gsl_complex *)x1,
                                 (gsl_complex *)x2);
}

template <typename T>
class complex_solve_cubic_functor
{
  public:
    typedef Array<std::complex<double>, 3, 1, ColMajor, 3, 1> ArrayType;

    complex_solve_cubic_functor(const T &c)
        : m_result(IEXP_NAN, IEXP_NAN, IEXP_NAN)
    {
        typename type_eval<T>::type m_c(c.eval());
        complex_solve_cubic_impl(c[0],
                                 c[1],
                                 c[2],
                                 &m_result[0],
                                 &m_result[1],
                                 &m_result[2]);
    }

    const typename std::complex<double> &operator()(Index i) const
    {
        return m_result(i);
    }

  private:
    ArrayType m_result;
};

template <typename T>
inline CwiseNullaryOp<complex_solve_cubic_functor<T>,
                      typename complex_solve_cubic_functor<T>::ArrayType>
complex_solve_cubic(const ArrayBase<T> &c)
{
    eigen_assert(IS_VEC(c) && (c.size() == 3));

    typedef typename complex_solve_cubic_functor<T>::ArrayType ArrayType;
    return ArrayType::NullaryExpr(3,
                                  complex_solve_cubic_functor<T>(c.derived()));
}

////////////////////////////////////////////////////////////
// global variants
////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////
// interface declaration
////////////////////////////////////////////////////////////
}

IEXP_NS_END

#endif /* __IEXP_POLY_CUBIC__ */
