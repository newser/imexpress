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

#ifndef __IEXP_RANDIST_LAPLACE__
#define __IEXP_RANDIST_LAPLACE__

////////////////////////////////////////////////////////////
// import header files
////////////////////////////////////////////////////////////

#include <common/common.h>

#include <gsl/gsl_cdf.h>
#include <gsl/gsl_randist.h>

IEXP_NS_BEGIN

namespace rdist {

////////////////////////////////////////////////////////////
// macro definition
////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////
// type definition
////////////////////////////////////////////////////////////

class laplace
{
  public:
    laplace(double a)
        : m_a(a)
    {
    }

    double pdf(double x) const
    {
        return gsl_ran_laplace_pdf(x, m_a);
    }

    double p(double x) const
    {
        return gsl_cdf_laplace_P(x, m_a);
    }

    double invp(double x) const
    {
        return gsl_cdf_laplace_Pinv(x, m_a);
    }

    double q(double x) const
    {
        return gsl_cdf_laplace_Q(x, m_a);
    }

    double invq(double x) const
    {
        return gsl_cdf_laplace_Qinv(x, m_a);
    }

  private:
    double m_a;
};

template <typename T>
class laplace_pdf_functor
{
  public:
    using ArrayType = Array<typename T::Scalar,
                            T::RowsAtCompileTime,
                            T::ColsAtCompileTime,
                            T::Flags & RowMajorBit ? RowMajor : ColMajor,
                            T::MaxRowsAtCompileTime,
                            T::MaxColsAtCompileTime>;

    laplace_pdf_functor(const T &x, typename T::Scalar a)
        : m_x(x)
        , m_laplace(a)
    {
    }

    typename T::Scalar operator()(Index i, Index j) const
    {
        return m_laplace.pdf(m_x(i, j));
    }

  private:
    const T &m_x;
    laplace m_laplace;
};

template <typename T>
inline CwiseNullaryOp<laplace_pdf_functor<T>,
                      typename laplace_pdf_functor<T>::ArrayType>
laplace_pdf(const ArrayBase<T> &x, typename T::Scalar a)
{
    static_assert(TYPE_IS(typename T::Scalar, double),
                  "scalar can only be double");

    using ArrayType = typename laplace_pdf_functor<T>::ArrayType;
    return ArrayType::NullaryExpr(x.rows(),
                                  x.cols(),
                                  laplace_pdf_functor<T>(x.derived(), a));
}

////////////////////////////////////////////////////////////
// global variants
////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////
// indexerface declaration
////////////////////////////////////////////////////////////
}

IEXP_NS_END

#endif /* __IEXP_RANDIST_LAPLACE__ */
