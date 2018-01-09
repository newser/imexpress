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

#ifndef __IEXP_RANDIST_NEG_BINOMIAL__
#define __IEXP_RANDIST_NEG_BINOMIAL__

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

class nbnom
{
  public:
    nbnom(double p, unsigned int n)
        : m_p(p)
        , m_n(n)
    {
    }

    double pdf(double x) const
    {
        return gsl_ran_negative_binomial_pdf(x, m_p, m_n);
    }

    double p(double x) const
    {
        return gsl_cdf_negative_binomial_P(x, m_p, m_n);
    }

    double q(double x) const
    {
        return gsl_cdf_negative_binomial_Q(x, m_p, m_n);
    }

  private:
    double m_p;
    unsigned int m_n;
};

template <typename T>
class nbnom_pdf_functor
{
  public:
    using ArrayType = Array<double,
                            T::RowsAtCompileTime,
                            T::ColsAtCompileTime,
                            T::Flags & RowMajorBit ? RowMajor : ColMajor,
                            T::MaxRowsAtCompileTime,
                            T::MaxColsAtCompileTime>;

    nbnom_pdf_functor(const T &x, double p, unsigned int n)
        : m_x(x)
        , m_nbnom(p, n)
    {
    }

    double operator()(Index i, Index j) const
    {
        return m_nbnom.pdf(m_x(i, j));
    }

  private:
    const T &m_x;
    nbnom m_nbnom;
};

template <typename T>
inline CwiseNullaryOp<nbnom_pdf_functor<T>,
                      typename nbnom_pdf_functor<T>::ArrayType>
nbnom_pdf(const ArrayBase<T> &x, double p, unsigned int n)
{
    static_assert(TYPE_IS(typename T::Scalar, int) ||
                      TYPE_IS(typename T::Scalar, unsigned int),
                  "scalar can only be int or unsigned int");

    using ArrayType = typename nbnom_pdf_functor<T>::ArrayType;
    return ArrayType::NullaryExpr(x.rows(),
                                  x.cols(),
                                  nbnom_pdf_functor<T>(x.derived(), p, n));
}

////////////////////////////////////////////////////////////
// global variants
////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////
// indexerface declaration
////////////////////////////////////////////////////////////
}

IEXP_NS_END

#endif /* __IEXP_RANDIST_NEG_BINOMIAL__ */
