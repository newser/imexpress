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

#ifndef __IEXP_RANDIST_GUMBEL1__
#define __IEXP_RANDIST_GUMBEL1__

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

class gbl1
{
  public:
    gbl1(double a, double b)
        : m_a(a)
        , m_b(b)
    {
    }

    double pdf(double x) const
    {
        return gsl_ran_gumbel1_pdf(x, m_a, m_b);
    }

    double p(double x) const
    {
        return gsl_cdf_gumbel1_P(x, m_a, m_b);
    }

    double invp(double x) const
    {
        return gsl_cdf_gumbel1_Pinv(x, m_a, m_b);
    }

    double q(double x) const
    {
        return gsl_cdf_gumbel1_Q(x, m_a, m_b);
    }

    double invq(double x) const
    {
        return gsl_cdf_gumbel1_Qinv(x, m_a, m_b);
    }

  private:
    double m_a, m_b;
};

template <typename T>
class gbl1_pdf_functor
{
  public:
    using ArrayType = Array<typename T::Scalar,
                            T::RowsAtCompileTime,
                            T::ColsAtCompileTime,
                            T::Flags & RowMajorBit ? RowMajor : ColMajor,
                            T::MaxRowsAtCompileTime,
                            T::MaxColsAtCompileTime>;

    gbl1_pdf_functor(const T &x, typename T::Scalar a, typename T::Scalar b)
        : m_x(x)
        , m_gbl1(a, b)
    {
    }

    typename T::Scalar operator()(Index i, Index j) const
    {
        return m_gbl1.pdf(m_x(i, j));
    }

  private:
    const T &m_x;
    gbl1 m_gbl1;
};

template <typename T>
inline CwiseNullaryOp<gbl1_pdf_functor<T>,
                      typename gbl1_pdf_functor<T>::ArrayType>
gbl1_pdf(const ArrayBase<T> &x, typename T::Scalar a, typename T::Scalar b)
{
    static_assert(TYPE_IS(typename T::Scalar, double),
                  "scalar can only be double");

    using ArrayType = typename gbl1_pdf_functor<T>::ArrayType;
    return ArrayType::NullaryExpr(x.rows(),
                                  x.cols(),
                                  gbl1_pdf_functor<T>(x.derived(), a, b));
}

////////////////////////////////////////////////////////////
// global variants
////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////
// indexerface declaration
////////////////////////////////////////////////////////////
}

IEXP_NS_END

#endif /* __IEXP_RANDIST_GUMBEL1__ */
