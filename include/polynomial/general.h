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

#ifndef __IEXP_POLY_GENERAL__
#define __IEXP_POLY_GENERAL__

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

template <typename T>
inline void complex_solve_impl(const T *a,
                               const int len,
                               std::complex<double> *z)
{
    UNSUPPORTED_TYPE(T);
}

template <>
inline void complex_solve_impl(const double *a,
                               const int len,
                               std::complex<double> *z)
{
    gsl_poly_complex_workspace *w = gsl_poly_complex_workspace_alloc(len);
    IEXP_NOT_NULLPTR(w);

    gsl_poly_complex_solve(a, len, w, (gsl_complex_packed_ptr)z);

    gsl_poly_complex_workspace_free(w);
}

template <typename T>
class complex_solve_functor
{
  public:
    using ArrayType = Array<std::complex<double>,
                            T::SizeAtCompileTime,
                            1,
                            ColMajor,
                            T::SizeAtCompileTime,
                            1>;

    complex_solve_functor(const T &a)
        : m_result(a.size() - 1, 1)
    {
        typename type_eval<T>::type m_a(a.eval());
        complex_solve_impl(m_a.data(), m_a.size(), m_result.data());
    }

    const typename std::complex<double> &operator()(Index i) const
    {
        return m_result(i);
    }

  private:
    ArrayType m_result;
};

template <typename T>
inline CwiseNullaryOp<complex_solve_functor<T>,
                      typename complex_solve_functor<T>::ArrayType>
complex_solve(const ArrayBase<T> &a)
{
    eigen_assert(IS_VEC(a));

    using ArrayType = typename complex_solve_functor<T>::ArrayType;
    return ArrayType::NullaryExpr(a.derived().size() - 1,
                                  complex_solve_functor<T>(a.derived()));
}

////////////////////////////////////////////////////////////
// global variants
////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////
// interface declaration
////////////////////////////////////////////////////////////
}

IEXP_NS_END

#endif /* __IEXP_POLY_GENERAL__ */
