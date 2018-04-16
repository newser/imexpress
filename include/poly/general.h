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

inline void solve(const double *a, int len, std::complex<double> *z)
{
    gsl_poly_complex_workspace *w = gsl_poly_complex_workspace_alloc(len);
    IEXP_NOT_NULLPTR(w);

    gsl_poly_complex_solve(a, len, w, (gsl_complex_packed_ptr)z);

    gsl_poly_complex_workspace_free(w);
}

inline void solve(const std::initializer_list<double> &a,
                  std::complex<double> *z)
{
    gsl_poly_complex_workspace *w = gsl_poly_complex_workspace_alloc(a.size());
    IEXP_NOT_NULLPTR(w);

    gsl_poly_complex_solve(a.begin(), a.size(), w, (gsl_complex_packed_ptr)z);

    gsl_poly_complex_workspace_free(w);
}

template <bool row_form, typename T>
class solve_functor
{
  public:
    using Scalar = typename TYPE_CHOOSE(IS_COMPLEX(typename T::Scalar),
                                        typename T::Scalar,
                                        std::complex<typename T::Scalar>);
    using ResultType = typename dense_derive<T,
                                             Scalar,
                                             row_form ? 1 : Dynamic,
                                             row_form ? Dynamic : 1>::type;

    solve_functor(const T &a)
        : m_result(new Scalar[a.size() - 1])
    {
        typename type_eval<T>::type m_a(a.eval());
        solve(m_a.data(), m_a.size(), m_result.get());
    }

    Scalar operator()(Index i) const
    {
        return m_result.get()[i];
    }

  private:
    std::shared_ptr<Scalar> m_result;
};

template <bool row_form = false, typename T = void>
inline CwiseNullaryOp<solve_functor<row_form, T>,
                      typename solve_functor<row_form, T>::ResultType>
solve(const DenseBase<T> &a)
{
    using ResultType = typename solve_functor<row_form, T>::ResultType;
    return ResultType::NullaryExpr(row_form ? 1 : (a.size() - 1),
                                   row_form ? (a.size() - 1) : 1,
                                   solve_functor<row_form, T>(a.derived()));
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
