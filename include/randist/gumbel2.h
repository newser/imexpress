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

#ifndef __IEXP_RANDIST_GUMBEL2__
#define __IEXP_RANDIST_GUMBEL2__

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

class gbl2
{
  public:
    gbl2(double a, double b)
        : m_a(a)
        , m_b(b)
    {
    }

    double pdf(double x) const
    {
        return gsl_ran_gumbel2_pdf(x, m_a, m_b);
    }

    double p(double x) const
    {
        return gsl_cdf_gumbel2_P(x, m_a, m_b);
    }

    double invp(double x) const
    {
        return gsl_cdf_gumbel2_Pinv(x, m_a, m_b);
    }

    double q(double x) const
    {
        return gsl_cdf_gumbel2_Q(x, m_a, m_b);
    }

    double invq(double x) const
    {
        return gsl_cdf_gumbel2_Qinv(x, m_a, m_b);
    }

  private:
    double m_a, m_b;
};

template <typename T>
class gbl2_pdf_functor
{
  public:
    using ArrayType = Array<typename T::Scalar,
                            T::RowsAtCompileTime,
                            T::ColsAtCompileTime,
                            T::Flags & RowMajorBit ? RowMajor : ColMajor,
                            T::MaxRowsAtCompileTime,
                            T::MaxColsAtCompileTime>;

    gbl2_pdf_functor(const T &x, typename T::Scalar a, typename T::Scalar b)
        : m_x(x)
        , m_gbl2(a, b)
    {
    }

    typename T::Scalar operator()(Index i, Index j) const
    {
        return m_gbl2.pdf(m_x(i, j));
    }

  private:
    const T &m_x;
    gbl2 m_gbl2;
};

template <typename T>
inline CwiseNullaryOp<gbl2_pdf_functor<T>,
                      typename gbl2_pdf_functor<T>::ArrayType>
gbl2_pdf(const ArrayBase<T> &x, typename T::Scalar a, typename T::Scalar b)
{
    static_assert(TYPE_IS(typename T::Scalar, double),
                  "scalar can only be double");

    using ArrayType = typename gbl2_pdf_functor<T>::ArrayType;
    return ArrayType::NullaryExpr(x.rows(),
                                  x.cols(),
                                  gbl2_pdf_functor<T>(x.derived(), a, b));
}

////////////////////////////////////////////////////////////
// global variants
////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////
// indexerface declaration
////////////////////////////////////////////////////////////
}

IEXP_NS_END

#endif /* __IEXP_RANDIST_GUMBEL2__ */
